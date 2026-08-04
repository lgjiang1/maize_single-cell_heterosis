#!/bin/bash
#SBATCH --job-name=popscle_pileup
#SBATCH --output=logs/popscle_pileup_%A_%a.out
#SBATCH --error=logs/popscle_pileup_%A_%a.err
#SBATCH --time=6-00:00:00    # it takes a lot of time
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-2

ml singularity

bam_files=(2_parent_genotyping_using_demuxlet/0_bam_file/*.non-ref_biased.bam)
bam_file=${bam_files[$((SLURM_ARRAY_TASK_ID-1))]}
sample_name=$(basename "$bam_file" .non-ref_biased.bam)
out_dir="2_parent_genotyping_using_demuxlet/1_demuxlet_results/${sample_name}"
mkdir -p $out_dir

## Before running popscle, check if "3parents_relative_to_B73v5_maf0_nomissing.vcf" vcf file includes AF, AC, and AN information in INFO field,
## If they’re missing, populate them with "bcftools +fill-tags" command: "bcftools +fill-tags 3parents_relative_to_B73v5_maf0_nomissing.vcf.gz   -Oz -o 3parents_relative_to_B73v5_maf0_nomissing.vcf.with_AF_AC_AN.vcf.gz  -- -t AF,AC,AN"

#======================================= Run pileup using popscle =============================================
singularity exec ../../../software/Demuxafy.sif popscle dsc-pileup \
        --sam ${bam_file} \
        --vcf 2_parent_genotyping_using_demuxlet/0_vcf_file/3parents_relative_to_B73v5_maf0_nomissing.vcf \
        --tag-group BC \
        --out ${out_dir}/pileup

