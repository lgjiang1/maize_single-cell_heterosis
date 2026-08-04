#!/bin/bash
#SBATCH --job-name=allele_ratio
#SBATCH --output=logs/allele_ratio_%A_%a.out
#SBATCH --error=logs/allele_ratio_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=140gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-18

vartrix_dir="4_genotyping_QC/5_Cell_validation/"
samples=($(find "${vartrix_dir}" -maxdepth 1 -type d -name "*_vartrix" | sed 's|.*/||' | sed 's/_vartrix$//' | sort))
sample=${samples[$((SLURM_ARRAY_TASK_ID - 1))]}
input_dir="${vartrix_dir}/${sample}_vartrix"

# =========== Run the Python script to compute allele ratio matrix ====================
python ./build_allele_ratio_matrix.py --input_dir $input_dir --sample_name $sample
