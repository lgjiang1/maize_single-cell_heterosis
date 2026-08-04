#!/bin/bash
#SBATCH --job-name=wasp
#SBATCH --output=logs/wasp_%A_%a.out
#SBATCH --error=logs/wasp_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=amarand1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --array=1-14

ml Bioinformatics bwa

WASP="../../../software/WASP"
bam_files=($(ls 1_WASP_analysis/0_bam_file/*.mq10.BC.rmdup.mm.bam))
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID-1]}
basename=$(basename $bam_file .mq10.BC.rmdup.mm.bam)
mkdir -p 1_WASP_analysis/1_WASP_bam

##========================== identify reads that may have mapping biases ===================================
#A minor modification was made to the original script (find_intersecting_snps.py), as it could not properly handle BAM files with incomplete pairing
python $WASP/mapping/modified_find_intersecting_snps.py \
    --is_paired_end \
    --is_sorted \
    --output_dir 1_WASP_analysis/1_WASP_bam \
    --snp_tab 1_WASP_analysis/0_hdf5_file/25NAM_tab.h5 \
    --snp_index 1_WASP_analysis/0_hdf5_file/25NAM_index.h5 \
    --haplotype 1_WASP_analysis/0_hdf5_file/haplotypes_25NAM.h5 \
    ${bam_file}


##=============================== remap reads using bwa ====================================================
reference="../../..//reference_genomes/zea/B73/B73v5_bwa"
fq1="1_WASP_analysis/1_WASP_bam/${basename}.mq10.BC.rmdup.mm.remap.fq1.gz"
fq2="1_WASP_analysis/1_WASP_bam/${basename}.mq10.BC.rmdup.mm.remap.fq2.gz"

mkdir -p 1_WASP_analysis/1_WASP_bam/remap

bwa mem -t 8 $reference $fq1 $fq2 | samtools view -@ 8 -b -q 10 - | samtools sort -@ 8 -o 1_WASP_analysis/1_WASP_bam/remap/${basename}.remap.sorted.bam
samtools index -@ 8 1_WASP_analysis/1_WASP_bam/remap/${basename}.remap.sorted.bam


##======================= filter out reads where the allelic versions of the reads fail to map back to the same location ===================================
python $WASP/mapping/filter_remapped_reads.py \
       1_WASP_analysis/1_WASP_bam/${basename}.mq10.BC.rmdup.mm.to.remap.bam \
       1_WASP_analysis/1_WASP_bam/remap/${basename}.remap.sorted.bam \
       1_WASP_analysis/1_WASP_bam/remap/${basename}.remap.sorted_keep.bam


##============================== merge ${SAMPLE_NAME}.keep.bam and ${SAMPLE_NAME}.remap.sorted_keep.bam ==================================
samtools merge -@ 8 \
	      1_WASP_analysis/1_WASP_bam/remap/${basename}.keep.merge.bam \
              1_WASP_analysis/1_WASP_bam/${basename}.mq10.BC.rmdup.mm.keep.bam  \
              1_WASP_analysis/1_WASP_bam/remap/${basename}.remap.sorted_keep.bam


samtools sort -@ 8 -o  1_WASP_analysis/1_WASP_bam/remap/${basename}.keep.merge.sorted.bam \
              1_WASP_analysis/1_WASP_bam/remap/${basename}.keep.merge.bam

samtools index -@ 8 1_WASP_analysis/1_WASP_bam/remap/${basename}.keep.merge.sorted.bam


##==================================== filter duplicate reads =================================================
mkdir -p 1_WASP_analysis/1_WASP_bam/final_bam

python $WASP/mapping/rmdup_pe.py \
       1_WASP_analysis/1_WASP_bam/remap/${basename}.keep.merge.sorted.bam \
       1_WASP_analysis/1_WASP_bam/final_bam/${basename}.keep.merge.sorted.rmdup.bam

##======================================= sort and index =========================================================
samtools sort -@ 8 -o /1_WASP_analysis/1_WASP_bam/final_bam/${basename}.non-ref_biased.bam \
	      1_WASP_analysis/1_WASP_bam/final_bam/${basename}.keep.merge.sorted.rmdup.bam

samtools index 1_WASP_analysis/1_WASP_bam/final_bam/${basename}.non-ref_biased.bam
rm 1_WASP_analysis/1_WASP_bam/final_bam/${basename}.keep.merge.sorted.rmdup.bam

echo 'mappability filtering finished'
