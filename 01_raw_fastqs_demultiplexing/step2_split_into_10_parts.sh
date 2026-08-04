#!/bin/bash
#SBATCH --job-name=fastq-splitting
#SBATCH --output=logs/astq-splitting_%A_%a.out
#SBATCH --error=logs/fastq-splitting_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=5gb
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-2   # Adjust the range based on the number of paired files

#This script splits the *bc1.bc2.fastq.gz file from step1 into 10 parts because the file is usually very big, and skipping this step would make step3 take too long.
#The split files are saved in the same folder (input_dir) as the original *bc1.bc2.fastq.gz file.

input_dir="02_clean_reads"
prefixes=($(ls ${input_dir}/*_R1_001.bc1.bc2.fastq.gz | sed 's/_R1_001.bc1.bc2.fastq.gz//' | xargs -n1 basename))
prefix=${prefixes[$SLURM_ARRAY_TASK_ID-1]}

r1=${input_dir}/${prefix}_R1_001.bc1.bc2.fastq.gz
r2=${input_dir}/${prefix}_R3_001.bc1.bc2.fastq.gz

#========================================= splitting ===========================================
ml Bioinformatics seqkit

seqkit split2 -p 10 -j 10 -1 $r1 -2 $r2 -O $input_dir
