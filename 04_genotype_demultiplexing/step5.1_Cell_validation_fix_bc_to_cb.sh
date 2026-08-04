#!/bin/bash
#SBATCH --job-name=BC_to_CB
#SBATCH --output=logs/bctocb_%A_%a.out
#SBATCH --error=logs/bctocb_%A_%a.err
#SBATCH --time=00:20:00
#SBATCH --mem=15gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-18

subset_bam_dir="4_genotyping_QC/1_subset_bam_file"
bams=($(ls ${subset_bam_dir}/*_subset.rg.bam))
bam=${bams[$SLURM_ARRAY_TASK_ID-1]}

filename=$(basename "$bam")
base=$(basename "$filename" _subset.rg.bam)

out_dir="4_genotyping_QC/5_Cell_validation"

mkdir -p ${out_dir}

# ======================= fix BC to CB tag ================================
samtools view -h "$bam" | \
  awk '{if($0 ~ /^@/) {print $0} else {gsub("BC:Z:", "CB:Z:", $0); print $0}}' | \
  samtools view -b -o ${out_dir}/"${base}_CB.bam"

samtools index ${out_dir}/"${base}_CB.bam"
echo "Done with $bam"
