#!/bin/bash
#SBATCH --job-name=rds_filtering
#SBATCH --output=logs/rds_filtering_%A_%a.out
#SBATCH --error=logs/rds_filtering_%A_%a.err
#SBATCH --time=00:05:00
#SBATCH --mem=7G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-14

# === Set and export environment variables ===
export RDS_DIR="./0_Socrates_raw_rds"
export META_DIR="./1_processed_Socrates_rds/extra_meta_info"
export OUT_DIR="./1_processed_Socrates_rds"

SAMPLES=($(ls ${RDS_DIR}/*.rds | sed 's|.*/||' | sed 's/.raw.soc.rds//' | sort))
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing sample: $SAMPLE_NAME"
echo "Using RDS:    ${RDS_DIR}/${SAMPLE_NAME}.raw.soc.rds"
echo "Using Meta:   ${META_DIR}/${SAMPLE_NAME}_barcode.txt"
echo "Saving to:    ${OUT_DIR}/"

# ================ Run the filter R script ===============
Rscript ./filter_rds.R $SAMPLE_NAME
