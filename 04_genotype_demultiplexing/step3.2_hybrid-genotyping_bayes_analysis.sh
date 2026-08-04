#!/bin/bash
#SBATCH --job-name=bayes_analysis
#SBATCH --output=logs/bayes_%A_%a.out
#SBATCH --error=logs/bayes_%A_%a.err
#SBATCH --time=00:01:00
#SBATCH --mem=1gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-12

barcode_dir=$(realpath 3_hybrid_genotyping_based_on_read_counts/1_read_count_file)
count_files=($(ls ${barcode_dir}/*_count.txt))
count_file=${count_files[$((SLURM_ARRAY_TASK_ID-1))]}
prefix=$(basename $count_file _count.txt)

outdir=$(realpath 3_hybrid_genotyping_based_on_read_counts/2_bayes_results)
mkdir -p "$outdir"

cd "$outdir"

Rscript ../../bayes_hybrid.R $count_file  $prefix  0.9  "${barcode_dir}/${prefix}.final_filtered_metadata.txt"

