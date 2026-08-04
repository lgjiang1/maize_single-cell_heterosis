#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --job-name=07a_fp_extract
#SBATCH --output=_logs/07a_fp_extract_%A_%a.log
#SBATCH --error=_logs/07a_fp_extract_%A_%a.log
#SBATCH --array=0-167

# Phase 3 Step 07a: Extract FP scores at motif hits with a coverage-
# stratified empirical null. Reads scPrinter's -log10(p) (default) or
# the legacy z-scores from 04c.
#
# Defaults:
#   --bg-percentile 20      Bottom-20% TFBS tiles as the null pool
#   --spread iqr            (P75-P25)/1.349; non-degenerate up to 75% zeros
#   --presence-z 2.0        Presence threshold (z_global > 2 = footprint)
#
# Optional env-var overrides (set on submission):
#   USABLE_PEAKS_DIR    Restrict to peaks from 07_0d (recommended)
#   OUTDIR              Output root (default 7_fp_extraction)
#   BG_PCT              Override --bg-percentile
#   SPREAD              Override --spread {iqr,mad,sd}
#
# Examples:
#   # B-K C5 quad on cleaned peak set (recommended)
#   sbatch --export=ALL,USABLE_PEAKS_DIR=regions/usable_peaks \
#          --array=9,23,51,65 \
#          07a_extract_fp_at_hits.sh
#
#   # Full run on all 168 samples
#   sbatch --export=ALL,USABLE_PEAKS_DIR=regions/usable_peaks \
#          07a_extract_fp_at_hits.sh


export HDF5_USE_FILE_LOCKING=FALSE

PROJ="10_TF_Footprinting"
cd "$PROJ"

echo "Job ID: ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

EXTRA_ARGS=()
if [[ -n "${USABLE_PEAKS_DIR:-}" ]]; then
  EXTRA_ARGS+=(--usable-peaks-dir "${USABLE_PEAKS_DIR}")
fi
if [[ -n "${OUTDIR:-}" ]]; then
  EXTRA_ARGS+=(--outdir "${OUTDIR}")
fi
if [[ -n "${BG_PCT:-}" ]]; then
  EXTRA_ARGS+=(--bg-percentile "${BG_PCT}")
fi
if [[ -n "${SPREAD:-}" ]]; then
  EXTRA_ARGS+=(--spread "${SPREAD}")
fi

python -u 07a_extract_fp_at_hits.py \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  "${EXTRA_ARGS[@]}"

echo ""
echo "End: $(date)"
