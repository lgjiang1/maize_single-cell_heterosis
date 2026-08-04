#!/bin/bash
#SBATCH --job-name=rds_merge
#SBATCH --output=logs/rds_merge_%j.out
#SBATCH --error=logs/rds_merge_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

rds_dir="1_processed_Socrates_rds"       # Folder containing *.filtered.soc.rds files
out_dir="1_processed_Socrates_rds"        # Where to save merged rds file

# ===== Run the merge R script =====
Rscript ./merge_rds.R ${rds_dir} ${out_dir}

