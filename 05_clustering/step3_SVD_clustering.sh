#!/bin/bash
#SBATCH --job-name=SVD_clustering
#SBATCH --output=logs/SVD-clustering_%j.out
#SBATCH --error=logs/SVD-clustering_%j.err
#SBATCH --time=01:30:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

rds="1_processed_Socrates_rds/merged.soc.rds"       # Folder containing merged.rds files
feature_rate=0.4
pcs=8
knear=30
res=1.0
out_dir="2_clustering_results"        # Where to save output results

mkdir -p ${out_dir}

# =================== Run the SVD clustering R script ====================
Rscript ./SVD_clustering_v2.R ${rds} ${feature_rate} ${pcs} ${knear} ${res} ${out_dir}

