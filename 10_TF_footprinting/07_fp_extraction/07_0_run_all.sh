#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --job-name=07_0_run_all
#SBATCH --output=_logs/07_0_run_all_%j.log
#SBATCH --error=_logs/07_0_run_all_%j.log

# 07_0 pipeline: random-coverage pre-cleaning of ACRs (4 steps end-to-end).
#
# Steps:
#   07_0a  random non-ACR 2kb windows per parent  (~10s, local-style)
#   07_0b  count frags per random window per sample  (~4-5 min/sample)
#   07_0c  aggregate per-sample percentiles  (~30s)
#   07_0d  derive usable peak set per (cross × cell type)  (~30s)
#
# Positional args (avoid SLURM --export comma parsing):
#   $1  TASK_IDS    Comma-separated task IDs for 07_0b. Empty = all 168.
#   $2  CELL_TYPES  Comma-separated cell types for 07_0d. Empty = all 14.
#
# Usage:
#   # C5 quad test (4 samples, ~25 min total)
#   sbatch 00_scripts/07_fp_extraction/07_0_run_all.sh 9,23,51,65 C5
#
#   # Full pipeline (all 168 samples, ~5h total)
#   sbatch 00_scripts/07_fp_extraction/07_0_run_all.sh
#
# After this job: submit 07a_v2 with --usable-peaks-dir to score on the
# cleaned ACR set.

PROJ="10_TF_footprinting"
cd "$PROJ"

TASK_IDS="${1:-}"
CELL_TYPES="${2:-}"

echo "================================================================"
echo " 07_0 pipeline — random-coverage pre-cleaning"
echo "================================================================"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       $(hostname)"
echo "Start:      $(date)"
echo "TASK_IDS:   ${TASK_IDS:-<all 168>}"
echo "CELL_TYPES: ${CELL_TYPES:-<all 14>}"
echo

# ----------------------------------------------------------------
# 07_0a — random non-ACR 2kb windows per parent
# ----------------------------------------------------------------
echo "----------------------------------------------------------------"
echo " 07_0a  Generate random non-ACR 2kb windows"
echo "----------------------------------------------------------------"
python -u 07_0a_make_random_regions.py
echo

# ----------------------------------------------------------------
# 07_0b — count fragments per random window per sample
# ----------------------------------------------------------------
echo "----------------------------------------------------------------"
echo " 07_0b  Count fragments per random window per sample"
echo "----------------------------------------------------------------"
if [[ -n "${TASK_IDS}" ]]; then
  python -u 07_0b_count_random_coverage.py \
    --task-ids "${TASK_IDS}"
else
  python -u 07_0b_count_random_coverage.py
fi
echo

# ----------------------------------------------------------------
# 07_0c — per-sample percentiles
# ----------------------------------------------------------------
echo "----------------------------------------------------------------"
echo " 07_0c  Aggregate per-sample random-coverage percentiles"
echo "----------------------------------------------------------------"
python -u 07_0c_aggregate_random_summary.py
echo

# ----------------------------------------------------------------
# 07_0d — usable peak sets per (cross × cell type)
# ----------------------------------------------------------------
echo "----------------------------------------------------------------"
echo " 07_0d  Derive usable peak set per (cross × cell type)"
echo "----------------------------------------------------------------"
if [[ -n "${CELL_TYPES}" ]]; then
  python -u 07_0d_build_usable_peaks.py \
    --only-cell-types "${CELL_TYPES}"
else
  python -u 07_0d_build_usable_peaks.py
fi

echo
echo "================================================================"
echo " End: $(date)"
echo "================================================================"
