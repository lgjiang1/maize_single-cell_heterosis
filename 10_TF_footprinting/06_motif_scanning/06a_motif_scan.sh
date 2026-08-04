#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=150G
#SBATCH --time=06:00:00
#SBATCH --output=_logs/06a_scan_%x_%A_%a.log
#SBATCH --error=_logs/06a_scan_%x_%A_%a.log
#SBATCH --array=0-49

# Phase 3 Step 06a: MOODS motif scanning.
#
# Parameterized SLURM array template — submit via 06a_submit_all.sh which
# passes per-scan configuration as environment variables:
#   SCAN_FASTA   — per-parent genome FASTA
#   SCAN_BED     — per-coordinate-system region BED (2000bp)
#   SCAN_OUTBASE — output base directory for this scan
#
# Array job: splits BED into 50 chunks, scans each with MOODS.
# Output per chunk: ${SCAN_OUTBASE}/chunk_NN/motif_hits.tsv.gz

export HDF5_USE_FILE_LOCKING=FALSE

PROJ="10_TF_footprinting"
cd "$PROJ"

# -- Per-scan configuration (set by 06a_submit_all.sh via --export) --
FASTA="${SCAN_FASTA:?SCAN_FASTA not set}"
BED_FULL="${SCAN_BED:?SCAN_BED not set}"
OUTBASE="${SCAN_OUTBASE:?SCAN_OUTBASE not set}"
N_CHUNKS=50

SPLIT_DIR="${OUTBASE}/splits"
PREFIX="${SPLIT_DIR}/region_part_"

mkdir -p "$SPLIT_DIR" _logs "$OUTBASE"

# Task 0 creates BED splits; others wait
if [[ "${SLURM_ARRAY_TASK_ID}" == "0" ]]; then
  if [[ ! -f "${PREFIX}00.bed" ]]; then
    echo "[INFO] Creating BED splits..." >&2
    rm -f "${PREFIX}"*.bed
    sort -k1,1 -k2,2n "$BED_FULL" > "${SPLIT_DIR}/regions_sorted.bed"
    split -d -a 2 -n l/${N_CHUNKS} "${SPLIT_DIR}/regions_sorted.bed" "$PREFIX"
    for f in "${PREFIX}"*; do
      [[ "$f" == *.bed ]] && continue
      mv "$f" "$f.bed"
    done
    echo "[INFO] BED splits ready." >&2
  else
    echo "[INFO] BED splits already exist; skipping split." >&2
  fi
fi

# Barrier: wait until splits exist
while [[ ! -f "${PREFIX}00.bed" ]]; do
  sleep 2
done

CHUNK=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")
BED="${PREFIX}${CHUNK}.bed"
OUT="${OUTBASE}/chunk_${CHUNK}"

echo "[INFO] Scan: ${SLURM_JOB_NAME}"
echo "[INFO] Chunk ${CHUNK}: scanning ${BED}"
echo "[INFO] FASTA: ${FASTA}"
echo "Job ID: ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

python -u 00_scripts/06_motif_scanning/06a_motif_scan.py \
  --fasta "$FASTA" \
  --bed-in "$BED" \
  --outdir "$OUT" \
  --meme "5_motif_scanning/signatures/Plant_Motif_SignatureDB.meme" \
  --metadata "5_motif_scanning/signatures/signature_metadata.tsv" \
  --n-jobs 2 \
  --chunk-size 100 \
  --pvalue 5e-5

echo ""
echo "End: $(date)"
