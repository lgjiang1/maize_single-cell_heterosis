#!/bin/bash
#SBATCH --job-name=snp2h5
#SBATCH --output=logs/snp2h5.out
#SBATCH --error=logs/snp2h5.err
#SBATCH --time=01:30:00
#SBATCH --mem=4gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

WASP="../../../software/WASP"
mkdir -p 1_WASP_analysis/0_hdf5_file


#=======================Run snp2h5 to convert VCF to HDF5 =======================
$WASP/snp2h5/snp2h5 --chrom ./maize_chr \
    --format vcf \
    --snp_index 1_WASP_analysis/0_hdf5_file/25NAM_index.h5 \
    --snp_tab 1_WASP_analysis/0_hdf5_file/25NAM_tab.h5 \
    --haplotype 1_WASP_analysis/0_hdf5_file/haplotypes_25NAM.h5 \
    1_WASP_analysis/0_vcf_file/25NAM_inbreds_relative_to_B73v5_chr*.vcf.gz
