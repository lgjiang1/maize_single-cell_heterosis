#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --array=0-2
#SBATCH --job-name=02a_bias
#SBATCH --partition=standard
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --output=_logs/02a_bias_%A_%a.log

# 02a_make_bias.sh — SLURM array (3 tasks) for the per-parent bias models.
#
# Each task runs 02a_make_bias.py for one parent (B73, Ki3, Oh43). Tasks
# are independent and can run in parallel. Each task uses the pretrained
# scPrinter Tn5 bias CNN (no training data needed beyond the per-parent
# fasta produced by 01a_split_genomes.py).
#
# Resource budget mirrors `bias_model/make_maize_bias.sh` (the original
# scaffold) and `13_Arbidopsis_protoplast/5_TF_FP/v1/3_00_Create_bias_and_genome_obj.sh`.
#
# Usage (from project root, after 01a_split_genomes.sh finishes):
#     sbatch 02a_make_bias.sh

# Tame native threads & HDF5 locking (matches reference setup)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

# Map task index → parent name
PARENTS=(B73 Ki3 Oh43)
PARENT=${PARENTS[$SLURM_ARRAY_TASK_ID]}

echo "[task ${SLURM_ARRAY_TASK_ID}] building bias model for ${PARENT}"

python 02a_make_bias.py --parent "${PARENT}"
