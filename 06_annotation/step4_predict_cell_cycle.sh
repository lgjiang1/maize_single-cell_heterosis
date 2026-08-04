#!/bin/bash
#SBATCH --job-name=cycle_prediction
#SBATCH --output=logs/cycle_predict-328.out
#SBATCH --error=logs/cycle_predict-328.err
#SBATCH --time=01:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

geneact="3_predict_cell_cycle/merged.scaled_gene_activity.rds"
marker="3_predict_cell_cycle/maize_cell_cycle_markers.txt"
output="3_predict_cell_cycle/cell_cycle_prediction_result.txt"

# ============= Run cell cycle phase prediction R script =======================
Rscript ./cell_cycle_prediction.R $geneact $marker $output
