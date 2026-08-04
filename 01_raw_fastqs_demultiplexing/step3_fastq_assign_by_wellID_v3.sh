#!/bin/bash
#SBATCH --job-name=fastq_assign
#SBATCH --output=logs/fastq_assign_%A_%a.out
#SBATCH --error=logs/fastq_assign_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=50mb
#SBATCH --time=2-00:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-40  ##Adjust the range based on the number of files

input_dir="02_clean_reads"  ### specify the input path for the clean_reads files (generated from step2)
outdir="02_clean_reads/Assign"   ### specify the output path for the split FASTQ files
layout_file="./96well_Tn5_bc_layout.txt"
sample_well_map="./96well_sample_layout.txt"
py_script="./fastq_sample_assign_and_fixTn5_v3.py"

mkdir -p $outdir

files=($(ls ${input_dir}/*part*.fastq.gz))
input_file=${files[$SLURM_ARRAY_TASK_ID-1]}

#============================== Run Python script =======================================
python ${py_script} ${layout_file} ${input_file} ${sample_well_map} --output ${outdir}
