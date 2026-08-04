#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-167
#SBATCH --job-name=01b_split
#SBATCH --partition=standard
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --output=_logs/01b_split_%A_%a.log

# 01b_split_and_rename.sh — SLURM array wrapper for 01b_split_and_rename.py
#
# Reads the per-task (input, target_parent, output) tuple from
# `_logs/01b_jobs.tsv`, which must have been generated first by:
#
#     python 01b_make_jobs.py
#
# Usage (from project root):
#     python 01b_make_jobs.py        # one-time, generates _logs/01b_jobs.tsv
#     sbatch 01b_split_and_rename.sh # 168-task array


# Tame native threads (matches reference Arabidopsis pipeline)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

JOBS_TSV=_logs/01b_jobs.tsv
if [ ! -f "$JOBS_TSV" ]; then
    echo "[ERROR] $JOBS_TSV not found. Run 01b_make_jobs.py first." >&2
    exit 1
fi

# 1-based line index for sed (SLURM_ARRAY_TASK_ID is 0-based)
LINE=$((SLURM_ARRAY_TASK_ID + 1))
ROW=$(sed -n "${LINE}p" "$JOBS_TSV")
if [ -z "$ROW" ]; then
    echo "[ERROR] empty row at line $LINE in $JOBS_TSV" >&2
    exit 1
fi

INPUT_BAM=$(echo "$ROW" | cut -f1)
TARGET_PARENT=$(echo "$ROW" | cut -f2)
OUTPUT_BAM=$(echo "$ROW" | cut -f3)

echo "[task ${SLURM_ARRAY_TASK_ID}] $INPUT_BAM"
echo "[task ${SLURM_ARRAY_TASK_ID}]   -> $OUTPUT_BAM (parent=$TARGET_PARENT)"

python 01b_split_and_rename.py \
    --input  "$INPUT_BAM" \
    --output "$OUTPUT_BAM" \
    --target-parent "$TARGET_PARENT"
