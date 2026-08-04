#!/bin/bash
#SBATCH --job-name=plot_markers
#SBATCH --output=logs/plot_markers_%j.out
#SBATCH --error=logs/plot_markers_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=180G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

input_dir="2_plot_marker_accessibility_scores/input_file"
threads=25  
prefix="merged"
output_dir="2_plot_marker_accessibility_scores"

mkdir -p  ${output_dir}

# ============= Run plot marker gene R script =======================
Rscript ./get_marker_accessibility.R \
	${input_dir}/merged.clustering.txt \
	${input_dir}/merged.scaled_gene_activity.rds \
	${input_dir}/merged.reduced_dimensions.txt \
	${input_dir}/maize_marker_genes_v5_update_2026.3.10.bed \
        $threads \
        $prefix \
	${output_dir}


## maize_marker_genes_v5_update_2026.3.10.bed file format:

#chr     start   end     geneID  name    type
#9       25172386        25177273        Zm00001eb378150 NST1    abaxial_bundle_sheath_cell
#1       1025166 1028066 Zm00001eb000240 STP3    abaxial_bundle_sheath_cell
#1       281720679       281723151       Zm00001eb056900 UMI12   abaxial_bundle_sheath_cell
#7       49104779        49117584        Zm00001eb307070 FRA1    abaxial_bundle_sheath_cell
#2       59030105        59032522        Zm00001eb083840 ACA1    bundle_sheath
#1       185941269       185953288       Zm00001eb033390 DCT2    bundle_sheath



