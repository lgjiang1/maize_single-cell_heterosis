#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name=01c_index
#SBATCH --partition=standard
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --output=_logs/01c_index_%A.log

# 01c_index_bams.sh — index all split BAMs produced by 01b_split_and_rename.
#
# Walks 1_split/{B-K,B-O}/{B73,Ki3,Oh43}_coords/*.bam and runs `samtools index`
# in parallel via xargs. Single SLURM job, ~hour wall time for 168 BAMs.
#
# Usage (from project root, after 01b finishes):
#     sbatch 01_split/01c_index_bams.sh

PARALLEL=8        # number of concurrent samtools index calls
THREADS_EACH=2    # threads per call → 8 * 2 = 16 cores total (matches --cpus-per-task)

echo "[INFO] indexing all BAMs under 1_split/"
echo "[INFO] parallel = $PARALLEL, threads-per-call = $THREADS_EACH"

# Collect BAM list
BAM_LIST=$(find 1_split/B-K 1_split/B-O -type f -name '*.bam')
N_BAMS=$(echo "$BAM_LIST" | wc -l)
echo "[INFO] $N_BAMS BAMs to index"

# Index in parallel
echo "$BAM_LIST" | xargs -n1 -P $PARALLEL samtools index -@ $THREADS_EACH

# Verify every BAM has a sibling .bai
MISSING=0
for bam in $BAM_LIST; do
    if [ ! -f "${bam}.bai" ]; then
        echo "[ERROR] missing index: ${bam}.bai" >&2
        MISSING=$((MISSING + 1))
    fi
done

if [ "$MISSING" -gt 0 ]; then
    echo "[ERROR] $MISSING BAMs missing indices" >&2
    exit 1
fi

echo "[OK] all $N_BAMS BAMs indexed"
