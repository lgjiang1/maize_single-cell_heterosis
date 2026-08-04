#!/bin/bash
#SBATCH --job-name=run_demuxlet
#SBATCH --output=logs/run_demuxlet_%A_%a.out
#SBATCH --error=logs/run_demuxlet_%A_%a.err
#SBATCH --time=00:15:00
#SBATCH --mem=80gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-2

ml singularity

bam_files=(2_parent_genotyping_using_demuxlet/0_bam_file/*.non-ref_biased.bam)
bam_file=${bam_files[$((SLURM_ARRAY_TASK_ID-1))]}
sample_name=$(basename "$bam_file" .non-ref_biased.bam)

singularity exec ../../../software/Demuxafy.sif popscle demuxlet \
        --plp 2_parent_genotyping_using_demuxlet/1_demuxlet_results/${sample_name}/pileup \
        --vcf 2_parent_genotyping_using_demuxlet/0_vcf_file/3parents_relative_to_B73v5_maf0_nomissing.vcf \
        --field GT \
        --out 2_parent_genotyping_using_demuxlet/1_demuxlet_results/${sample_name}/${sample_name}_demuxlet
