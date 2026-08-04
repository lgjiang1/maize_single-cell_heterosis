#!/bin/bash
#SBATCH --job-name=motif_cluster
#SBATCH --account=amarand1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/02a_motif_clustering_%j.log
#SBATCH --error=_logs/02a_motif_clustering_%j.log

export PROJ_ROOT="./10_TF_footprinting"
cd "$PROJ_ROOT"

echo "[$(date)] Starting motif clustering"

Rscript 02a_motif_clustering.R

