#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=04d_cov
#SBATCH --output=_logs/04d_count_acr_coverage_%j.log
#SBATCH --error=_logs/04d_count_acr_coverage_%j.log

# 04d — Count fragments per ACR per condition (168 samples, sequential).
#
# Reads fragment files + region BEDs, outputs per-ACR coverage table
# used for depth correction: δ_hat = z/√n
#
# Usage:
#   sbatch 04d_count_acr_coverage.sh


PROJ="./TF_footprinting_analysis"
cd $PROJ

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

python -u 04d_count_acr_coverage.py

echo ""
echo "End: $(date)"
