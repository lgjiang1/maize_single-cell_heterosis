#!/bin/bash
#SBATCH --job-name=step1
#SBATCH --output=logs/scATAC-QC_step1_%A_%a.out
#SBATCH --error=logs/scATAC-QC_step1_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=150gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-14

ml Bioinformatics py-macs2

bed_files=(03_mapping/3_final_bam_and_tn5bed/*.tn5.bed.gz)
bed=${bed_files[$SLURM_ARRAY_TASK_ID - 1]}
out=$(basename $bed .tn5.bed.gz)
ann="../../../reference_genomes/zea/B73/Zm-B73-REFERENCE-NAM-5.0-chrs-mt-pt.gff3.gz"
ref="./maize_chr"

Rscript ./QC_scATAC_data.v2.R $bed $out $ann $ref



