#!/bin/bash
#SBATCH --job-name=vartrix
#SBATCH --output=logs/vartrix_%A_%a.out
#SBATCH --error=logs/vartrix_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-18

input_dir="4_genotyping_QC/5_Cell_validation"
bams=($(ls ${input_dir}/*_CB.bam))
bam=${bams[$SLURM_ARRAY_TASK_ID-1]}

filename=$(basename "$bam")
sample=$(basename "$filename" _CB.bam)

barcode_file="${input_dir}/${sample}_barcode.txt"
ref="../../../reference_genomes/zea/B73/raw_B73v5_genome_file/Zm-B73-REFERENCE-NAM-5.0.fa"
vcf="${input_dir}/Final_Ref_scATAC_merged_SNP.vcf.gz"
outdir="${input_dir}/${sample}_vartrix"

mkdir -p $outdir

### Note: vartrix output variant is 0-based coordinates, while VCF is 1-based.
# ======================= run vartrix ================================
vartrix \
  --bam $bam \
  --bam-tag CB \
  --cell-barcodes ${barcode_file} \
  --fasta $ref \
  --vcf $vcf \
  --out-matrix "$outdir/alt.mtx" \
  --ref-matrix "$outdir/ref.mtx" \
  --out-variants "$outdir/variants.tsv" \
  --scoring-method coverage \
  --mapq 30 \
  --threads 10  

