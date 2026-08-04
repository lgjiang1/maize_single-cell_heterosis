#!/bin/bash
#SBATCH --job-name=step2
#SBATCH --output=logs/QC-step2_%A_%a.out
#SBATCH --error=logs/QC-step2_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=30gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-14

file_dir="Socrates_results"
files=(${file_dir}/*.raw.soc.rds)
rdsfile=${files[$SLURM_ARRAY_TASK_ID - 1]}
prefix=${rdsfile%.raw.soc.rds}

Rscript ./filter_low_qual_cells.R $rdsfile ${file_dir}/$prefix
