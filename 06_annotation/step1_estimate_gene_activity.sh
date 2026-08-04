#!/bin/bash
#SBATCH --job-name=estimate_gene_activity
#SBATCH --output=logs/estimate_gene_activity.out
#SBATCH --error=logs/estimate_gene_activity.err
#SBATCH --time=06:00:00
#SBATCH --mem=260G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=largemem

input_dir="1_estimate_gene_activity/input_file"
output_dir="1_estimate_gene_activity"

# ============= Run estimate_gene_activity R script =======================
Rscript ./estimate_gene_activity.R \
	${input_dir}/all_high-quality_genotyped_clustered_tn5.bed.gz \
	${input_dir}/Zm-B73-REFERENCE-NAM-5.0-chrs-mt-pt.gff3.gz \
	${input_dir}/merged.clustering.txt \
	${input_dir}/all_high-quality_genotyped_clustered_peaks.narrowPeak \
	${output_dir}
