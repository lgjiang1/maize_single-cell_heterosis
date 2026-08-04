#!/bin/bash
#SBATCH --job-name=extract_mappable_regions
#SBATCH --output=logs/create_control_regions_%A_%a.out
#SBATCH --error=logs/create_control_regions_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=120gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-3 

ml Bioinformatics bwa

ref_dir="../../../reference_genomes/zea"
fq_dir="2_uniquely_mappable_regions/sim_reads"

mkdir -p ${fq_dir}

genomes=(
  "${ref_dir}/B73/Zm-B73-REFERENCE-NAM-5.0-chrs-mt-pt.fa"
  "${ref_dir}/Ki3/Zm-Ki3-REFERENCE-NAM-1.0-chrs-mt-pt.fa"
  "${ref_dir}/Oh43/Zm-Oh43-REFERENCE-NAM-1.0-chrs-mt-pt.fa"
)

genome=${genomes[$((SLURM_ARRAY_TASK_ID - 1))]}
prefix=$(basename "$genome" | sed 's/-REFERENCE.*//')


#======================================================================
##  Run wgsim to create synthetic reads -------------------------------
#======================================================================
wgsim -N 400000000 -1 150 -2 150 -d 300 $genome ${fq_dir}/${prefix}_R1.fq ${fq_dir}/${prefix}_R2.fq
pigz -p 24 ${fq_dir}/${prefix}_R1.fq
pigz -p 24 ${fq_dir}/${prefix}_R2.fq
echo "Finished wgsim for $genome"


#======================================================================
##  Map synthetic reads -----------------------------------------------
#======================================================================
indexs=(
  "${ref_dir}/B73/B73v5_bwa"
  "${ref_dir}/Ki3/Ki3_bwa"
  "${ref_dir}/Oh43/Oh43_bwa"
)

index=${indexs[$((SLURM_ARRAY_TASK_ID - 1))]}
bam_dir="2_uniquely_mappable_regions/bamfiles"
bed_dir="2_uniquely_mappable_regions/bedfiles"
mkdir -p $bam_dir
mkdir -p $bed_dir

raw_bam=${bam_dir}/${prefix}.raw.bam
mq10_uni_bam=${bam_dir}/${prefix}.mq10.unique.bam

# align synthetic reads
echo "BWA mapping to index: $index"
bwa mem -M -t 24 $index ${fq_dir}/${prefix}_R1.fq.gz ${fq_dir}/${prefix}_R2.fq.gz \
	| samtools view -bS - \
	| samtools sort -@ 24 - > $raw_bam
echo "Finished bwa mapping for $prefix"


# filter synthetic bam alignments (uniquely map: q > 10; no XA tag)
echo " filtering alignments ..."
samtools view -h -q 10 -f 3 $raw_bam \
       | grep -v -E -e '\bXA:Z:' \
       | samtools view -bSh - > $mq10_uni_bam
echo "Finished filtering for $prefix"


#======================================================================
##  Extract mappable genomic regions ----------------------------------
#======================================================================
echo " extracting framents ..."
samtools sort -@ 24 -n $mq10_uni_bam \
| bedtools bamtobed -bedpe -i - \
| awk 'BEGIN{OFS="\t"}
       $1==$4 {                        
         s = ($2<$5 ? $2 : $5);         
         e = ($3>$6 ? $3 : $6);        
         if (e > s && s >= 0) print $1, s, e
       }' \
| sort -k1,1V -k2,2n \
| bedtools merge -i - \
> ${bed_dir}/${prefix}.mappable_genomic_regions.bed

echo "Finished extracting mappable genomic regions for $prefix"



