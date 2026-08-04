#!/bin/bash
#SBATCH --job-name=cis_trans
#SBATCH --output=logs/cis_trans_%A_%a.out
#SBATCH --error=logs/cis_trans_%A_%a.err
#SBATCH --mail-user=lgjiang@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --time=02:00:00
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=amarand1
#SBATCH --partition=standard
#SBATCH --array=1-14  # cell type numbers

Rscript fit_cis_trans_model.R ${SLURM_ARRAY_TASK_ID}
