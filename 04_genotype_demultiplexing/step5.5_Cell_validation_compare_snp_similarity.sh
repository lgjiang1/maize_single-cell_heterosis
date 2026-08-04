#!/bin/bash
#SBATCH --job-name=compare_snp_similarity
#SBATCH --output=logs/snp_similarity_%A_%a.out
#SBATCH --error=logs/snp_similarity_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-68

wk_dir="4_genotyping_QC/5_Cell_validation"
chunk_id=${SLURM_ARRAY_TASK_ID}
input_file="${wk_dir}/chunks/output_part_${chunk_id}.txt"
output_file="${wk_dir}/chunks/similarity_result_${chunk_id}.tsv"

# ======================= Run compare_snp_similarity python script ================================
python ./compare_snp_similarity.py \
  --vcf_file ${wk_dir}/Final_Ref_SNP.vcf \
  --tsv_file ${input_file} \
  --output_file ${output_file} \
  --min_snps 100 \
  --chunk_size 500

