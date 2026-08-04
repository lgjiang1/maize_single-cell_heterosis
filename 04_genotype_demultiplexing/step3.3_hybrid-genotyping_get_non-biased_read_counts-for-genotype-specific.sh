#!/bin/bash
#SBATCH --job-name=counting_reads
#SBATCH --output=logs/counting_reads_%A_%a.out
#SBATCH --error=logs/counting_reads_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=25gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=amarand1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --array=1-12


#===================variables=================================
bam_dir="3_hybrid_genotyping_based_on_read_counts/0_bam_file"
vcf_dir="3_hybrid_genotyping_based_on_read_counts/0_vcf_file"

bam_files=(
    ${bam_dir}/LSC2_BxK.non-ref_biased.bam
    ${bam_dir}/LSC3_BxK.non-ref_biased.bam
    ${bam_dir}/LSC2_KxB.non-ref_biased.bam
    ${bam_dir}/LSC3_KxB.non-ref_biased.bam
    ${bam_dir}/LSC2_BxO.non-ref_biased.bam
    ${bam_dir}/LSC3_BxO.non-ref_biased.bam
    ${bam_dir}/LSC2_OxB.non-ref_biased.bam
    ${bam_dir}/LSC3_OxB.non-ref_biased.bam
    ${bam_dir}/LSC2_KxO.non-ref_biased.bam
    ${bam_dir}/LSC3_KxO.non-ref_biased.bam
    ${bam_dir}/LSC2_OxK.non-ref_biased.bam
    ${bam_dir}/LSC3_OxK.non-ref_biased.bam
)

vcf_files=(
    ${vcf_dir}/Oh43_specific.vcf
    ${vcf_dir}/Oh43_specific.vcf
    ${vcf_dir}/Oh43_specific.vcf
    ${vcf_dir}/Oh43_specific.vcf
    ${vcf_dir}/Ki3_specific.vcf
    ${vcf_dir}/Ki3_specific.vcf
    ${vcf_dir}/Ki3_specific.vcf
    ${vcf_dir}/Ki3_specific.vcf
    ${vcf_dir}/B73_specific.vcf
    ${vcf_dir}/B73_specific.vcf
    ${vcf_dir}/B73_specific.vcf
    ${vcf_dir}/B73_specific.vcf
)

bam_file=${bam_files[$((SLURM_ARRAY_TASK_ID-1))]}
vcf_file=${vcf_files[$((SLURM_ARRAY_TASK_ID-1))]}

bam_filename=$(basename $bam_file)
prefix=${bam_filename%.non-ref_biased.bam}
mkdir -p 3_hybrid_genotyping_based_on_read_counts/3_genotype-specific-read_count_file
count="../3_hybrid_genotyping_based_on_read_counts/3_genotype-specific-read_count_file/${prefix}_count.txt"

echo "Running count_reads.py for:"
echo "  BAM file: $bam_file"
echo "  VCF file: $vcf_file"
echo "  Output:   $count"

#================= get counts ==========================
python ./genotype_specific_reads_counting.py $bam_file $vcf_file $count




