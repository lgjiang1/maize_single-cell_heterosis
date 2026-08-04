#!/bin/bash
#SBATCH --job-name=01a_split_genomes
#SBATCH --output=logs/01a_split_genomes_%j.out
#SBATCH --error=logs/01a_split_genomes_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

# 01a_split_genomes.sh — split per-cross merged fastas into per-parent fastas.
# See 00_scripts/01_split/01a_split_genomes.py for behavior.
#
# Usage (from project root):
#   sbatch 00_scripts/01_split/01a_split_genomes.sh
#
# Output:
#   1_split/genomes/{B73,Ki3,Oh43}.fa[.fai]

# Tame native threads (matches reference Arabidopsis pipeline)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

python 01a_split_genomes.py
