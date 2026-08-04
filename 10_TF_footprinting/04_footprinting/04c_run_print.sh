#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --array=0-167
#SBATCH --job-name=04c_fp
#SBATCH --partition=standard
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --output=_logs/04c_fp_%A_%a.log

# 04c_run_print.sh — SLURM CPU array: scPrinter import + TFBS + NucBS + FP.
#
# Reads _logs/04b_jobs.tsv (same table as 04b).
# Requires: scprinter-gpu conda env (newest scPrinter; runs on CPU partition).
#
# Usage:
#     sbatch 04c_run_print.sh

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

# ── Read job table ───────────────────────────────────────────────────────
JOBS_TSV=_logs/04b_jobs.tsv
if [ ! -f "$JOBS_TSV" ]; then
    echo "[ERROR] $JOBS_TSV not found" >&2
    exit 1
fi

LINE=$((SLURM_ARRAY_TASK_ID + 2))
ROW=$(sed -n "${LINE}p" "$JOBS_TSV")
if [ -z "$ROW" ]; then
    echo "[ERROR] empty row at line $LINE" >&2
    exit 1
fi

CROSS=$(echo "$ROW" | cut -f1)
COORD_SYSTEM=$(echo "$ROW" | cut -f2)
SAMPLE_ID=$(echo "$ROW" | cut -f4)
FRAG_OUT=$(echo "$ROW" | cut -f6)
GENOME_OBJ=$(echo "$ROW" | cut -f7)
REGION_BED=$(echo "$ROW" | cut -f8)
FP_OUTDIR=$(echo "$ROW" | cut -f9)

# Extract cell type from sample_id (e.g. B-K_B73_B73-C5 -> C5,
# B-K_Ki3_BxK-C5.Ki3allele -> C5).
if [[ "$SAMPLE_ID" =~ -(C[0-9]+) ]]; then
    CELL_TYPE="${BASH_REMATCH[1]}"
else
    echo "[ERROR] could not parse cell type from sample_id $SAMPLE_ID" >&2
    exit 1
fi

# ── Optional pre-filter to coverage-usable peaks (from 07_0d) ───────────
# If a per-quad usable BED exists, restrict scPrinter to that subset to cut
# runtime + h5ad size. Override location via env var USABLE_PEAKS_DIR.
USABLE_PEAKS_DIR="${USABLE_PEAKS_DIR:-regions/usable_peaks}"
USABLE_BED="${USABLE_PEAKS_DIR}/${CROSS}/${CELL_TYPE}/${COORD_SYSTEM}_coord.bed"

SUBSET_ARGS=()
if [ -f "$USABLE_BED" ]; then
    SUBSET_ARGS+=(--regions-subset "$USABLE_BED")
    SUBSET_MSG="$USABLE_BED"
else
    SUBSET_MSG="(none; running on full $REGION_BED)"
fi

echo "[task ${SLURM_ARRAY_TASK_ID}] sample=$SAMPLE_ID"
echo "[task ${SLURM_ARRAY_TASK_ID}] cell_type=$CELL_TYPE  coord=$COORD_SYSTEM"
echo "[task ${SLURM_ARRAY_TASK_ID}] frag=$FRAG_OUT"
echo "[task ${SLURM_ARRAY_TASK_ID}] genome=$GENOME_OBJ"
echo "[task ${SLURM_ARRAY_TASK_ID}] regions=$REGION_BED"
echo "[task ${SLURM_ARRAY_TASK_ID}] regions-subset=$SUBSET_MSG"
echo "[task ${SLURM_ARRAY_TASK_ID}] outdir=$FP_OUTDIR"

# Verify fragment file exists
if [ ! -f "$FRAG_OUT" ]; then
    echo "[ERROR] fragment file missing: $FRAG_OUT" >&2
    exit 1
fi

# ── Run scPrinter ────────────────────────────────────────────────────────
python 00_scripts/04_footprinting/04c_run_print.py \
    --frag "$FRAG_OUT" \
    --genome-obj "$GENOME_OBJ" \
    --regions "$REGION_BED" \
    "${SUBSET_ARGS[@]}" \
    --outdir "$FP_OUTDIR" \
    --sample-id "$SAMPLE_ID"

# ── Dedup duplicate h5ads ────────────────────────────────────────────────
# scPrinter writes each output h5ad twice: as <PREFIX>__ALL.h5ad (the
# canonical name downstream 07a/09c read) AND as a sample-named copy
# <PREFIX>_<sample_id_normalized>.h5ad with identical bytes. The named
# copy is read by nothing in this pipeline. Replace it with a relative
# symlink to reclaim ~50 GB per FP file.
SUPP_DIR="$FP_OUTDIR/printer_supp"
if [ -d "$SUPP_DIR" ]; then
    for prefix in TFBS NucBS FP_neglog10p; do
        canonical="${SUPP_DIR}/${prefix}__ALL.h5ad"
        [ -f "$canonical" ] || continue
        for h5 in "${SUPP_DIR}/${prefix}_"*.h5ad; do
            [ -f "$h5" ] || continue
            [ -L "$h5" ] && continue
            base="$(basename "$h5")"
            [ "$base" = "${prefix}__ALL.h5ad" ] && continue
            rm -f "$h5"
            ln -s "${prefix}__ALL.h5ad" "$h5"
            echo "[dedup] $base -> ${prefix}__ALL.h5ad"
        done
    done
fi
