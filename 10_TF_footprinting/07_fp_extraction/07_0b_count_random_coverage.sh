#!/bin/bash
#SBATCH --account=amarand1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=07_0b_rand_cov
#SBATCH --output=_logs/07_0b_random_coverage_%j.log
#SBATCH --error=_logs/07_0b_random_coverage_%j.log

# 07_0b — Count fragments at random non-ACR 2kb windows per sample.
# Per-sample P95(random) is the "above-background" threshold used in 07_0d.
#
# Usage:
#   cd /nfs/turbo/lsa-amarand/fabio_home/Projects/15_Heterosis
#   # Full run (all 168 samples, ~4h)
#   sbatch 00_scripts/07_fp_extraction/07_0b_count_random_coverage.sh
#   # C5 quad test (~minutes)
#   sbatch --export=ALL,TASK_IDS=9,23,51,65 \
#          00_scripts/07_fp_extraction/07_0b_count_random_coverage.sh

set -euo pipefail

source ~/home_turbo/fabio_home/LocalInstall/miniconda3/etc/profile.d/conda.sh
conda activate scprinter-cpu

PROJ=/nfs/turbo/lsa-amarand/fabio_home/Projects/15_Heterosis
cd "$PROJ"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

if [[ -n "${TASK_IDS:-}" ]]; then
  python -u 00_scripts/07_fp_extraction/07_0b_count_random_coverage.py \
    --task-ids "${TASK_IDS}"
else
  python -u 00_scripts/07_fp_extraction/07_0b_count_random_coverage.py
fi

echo ""
echo "End: $(date)"
