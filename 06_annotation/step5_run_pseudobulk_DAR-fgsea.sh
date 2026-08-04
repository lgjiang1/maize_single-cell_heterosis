#!/bin/bash
#SBATCH --job-name=pseudobulk_DEG_fgsea
#SBATCH --output=logs/pseudobulk_DEG_fgsea.out
#SBATCH --error=logs/pseudobulk_DEG_fgsea.err
#SBATCH --time=00:10:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

raw_geneact="4_cell-type_validation/02_pseudobulk_DEG_fgsea/input/merged.raw_gene_activity.rds"
meta="4_cell-type_validation/02_pseudobulk_DEG_fgsea/input/merged.clustering.txt"
gmt="../../../maize_GO/Zm-B73-REFERENCE-NAM-5.0-GOTerms.gmt"
marker="4_cell-type_validation/02_pseudobulk_DEG_fgsea/input/maize_marker_genes_v5_update_2025.8.1.bed"
outdir="4_cell-type_validation/02_pseudobulk_DEG_fgsea"

# ============= Run fgsea R script =======================
Rscript ./pseudobulk_DAR-fgsea.R \
  --counts ${raw_geneact} \
  --meta $meta \
  --gmt $gmt \
  --minSize 10 \
  --maxSize 600 \
  --nperm 10000 \
  --markers $marker \
  --out $outdir
