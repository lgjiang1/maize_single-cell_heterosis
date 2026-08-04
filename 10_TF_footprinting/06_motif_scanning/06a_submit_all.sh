#!/bin/bash
# Submit all 4 motif scans as separate SLURM array jobs.
#
# Each scan targets one (cross x parent-coordinate-system) combination,
# using the same per-parent FASTA and region BED used by 04c scPrinter scoring.
#
# Usage:
#   cd 10_TF_footprinting
#   bash 06a_submit_all.sh


SCRIPT="06a_motif_scan.sh"

echo "Submitting 4 motif scans..."

# B-K / B73 coordinates
sbatch --export=ALL,SCAN_FASTA=1_split/genomes/B73.fa,SCAN_BED=regions/B-K/B73_coords/regions_2000bp.bed,SCAN_OUTBASE=5_motif_scanning/hits/B-K/B73_coords \
  --job-name=06a_BK_B73 "$SCRIPT"

# B-K / Ki3 coordinates
sbatch --export=ALL,SCAN_FASTA=1_split/genomes/Ki3.fa,SCAN_BED=regions/B-K/Ki3_coords/regions_2000bp.bed,SCAN_OUTBASE=5_motif_scanning/hits/B-K/Ki3_coords \
  --job-name=06a_BK_Ki3 "$SCRIPT"

# B-O / B73 coordinates
sbatch --export=ALL,SCAN_FASTA=1_split/genomes/B73.fa,SCAN_BED=regions/B-O/B73_coords/regions_2000bp.bed,SCAN_OUTBASE=5_motif_scanning/hits/B-O/B73_coords \
  --job-name=06a_BO_B73 "$SCRIPT"

# B-O / Oh43 coordinates
sbatch --export=ALL,SCAN_FASTA=1_split/genomes/Oh43.fa,SCAN_BED=regions/B-O/Oh43_coords/regions_2000bp.bed,SCAN_OUTBASE=5_motif_scanning/hits/B-O/Oh43_coords \
  --job-name=06a_BO_Oh43 "$SCRIPT"

echo "All 4 scans submitted."
