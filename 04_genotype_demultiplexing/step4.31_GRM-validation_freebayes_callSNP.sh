#!/bin/bash
#SBATCH --job-name=freebayes
#SBATCH --output=logs/freebayes_%j.out
#SBATCH --error=logs/freebayes_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80gb
#SBATCH --time=04:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

ref="../../../reference_genomes/zea/B73/raw_B73v5_genome_file/Zm-B73-REFERENCE-NAM-5.0.fa"

mkdir -p 4_genotyping_QC/3_GRM_validation
cd 4_genotyping_QC/1_subset_bam_file
ls *.rg.bam > bam.list

#=================================== Run FreeBayes in parallel across 24 cores ====================================
freebayes-parallel \
    <(fasta_generate_regions.py ${ref}.fai 5000000) 24 -f $ref -L bam.list \
    --limit-coverage 10000 -q 10  -n 4 \
    --strict-vcf > ../3_GRM_validation/high-quality_and_genotyped_cells.vcf


