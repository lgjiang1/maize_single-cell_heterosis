#!/bin/bash
#SBATCH --job-name=plot_markers
#SBATCH --output=logs/plot_markers.out
#SBATCH --error=logs/plot_markers.err
#SBATCH --time=00:05:00
#SBATCH --mem=7G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

rds="2_plot_marker_accessibility_scores/merged_marker_clusters.rds"
marker="2_plot_marker_accessibility_scores/input_file//maize_marker_genes_v5_update_2026.3.10.bed"
gene="2_plot_marker_accessibility_scores/all_markers.txt"
outdir="2_plot_marker_accessibility_scores"

# ============= Run plot marker gene R script =======================
# Note: $gene can be either (1) a file containing one gene name per line (e.g. all_markers.txt), 
# or (2) a comma-separated string of gene names (e.g. "tb1" or "tb1,dct2").
#Rscript ./visualize_marker_accessibility.R $rds $marker $gene $outdir markerActivityScores 

Rscript ./visualize_marker_accessibility.R $rds $marker $gene $outdir seedling_with_genename
#Rscript ./visualize_marker_accessibility.R $rds $marker "BK1" $outdir BK1
