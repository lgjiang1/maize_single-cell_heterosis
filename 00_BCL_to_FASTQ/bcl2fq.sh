#!/bin/bash
#SBATCH --job-name=bcl2fastq
#SBATCH --output=logs/bcl2fastq.out
#SBATCH --error=logs/bcl2fastq.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=30gb
#SBATCH --time=20:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

ml Bioinformatics  bcl2fastq2/2.20.0.422-oxq6lf3

# vars
bcldir="RawData/BCLs/admera_24350_01_20241224_LH00150_0650_A22HYMNLT4"
rawfastq="01_raw_fastqs"

mkdir -p $rawfastq

##==================== Run bcl2fastq ====================
bcl2fastq --use-bases-mask=Y150n*,I8n*,Y16n*,Y150n* \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          --no-lane-splitting \
          --barcode-mismatches 1 \
          -r 15 -w 15 -p 15 \
          -R $bcldir \
          --output-dir=$rawfastq \
          --interop-dir=${bcldir}/InterOp \
          --sample-sheet=./SampleSheet.csv

