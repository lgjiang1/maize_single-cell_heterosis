#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --job-name=03b_idx
#SBATCH --partition=standard
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --output=_logs/03b_index_%j.log

# 03b_index_equalized.sh — Index any equalized BAMs that lack a .bai.
#
# Symlinked BAMs already have symlinked indices from 03a. This step catches
# any that were missed and verifies completeness.
#
# Usage:
#     sbatch --dependency=afterok:$JOB_03A 03b_index_equalized.sh


EQ_DIR=2_equalized
MISSING=0
INDEXED=0
TOTAL=0

for bam in $(find "$EQ_DIR" -name "*.bam" -type f -o -name "*.bam" -type l | sort); do
    TOTAL=$((TOTAL + 1))
    if [ ! -f "${bam}.bai" ] && [ ! -L "${bam}.bai" ]; then
        echo "[INDEX] $bam"
        samtools index "$bam"
        INDEXED=$((INDEXED + 1))
    fi
done

echo "[OK] total=$TOTAL indexed=$INDEXED already_indexed=$((TOTAL - INDEXED))"

# Verify all BAMs have indices
for bam in $(find "$EQ_DIR" -name "*.bam" -type f -o -name "*.bam" -type l | sort); do
    if [ ! -f "${bam}.bai" ] && [ ! -L "${bam}.bai" ]; then
        echo "[ERROR] missing index: ${bam}.bai" >&2
        MISSING=$((MISSING + 1))
    fi
done

if [ "$MISSING" -gt 0 ]; then
    echo "[ERROR] $MISSING BAMs still lack indices" >&2
    exit 1
fi

echo "[OK] all $TOTAL BAMs indexed"
