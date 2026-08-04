#!/bin/bash
#SBATCH --job-name=SNP_processing
#SBATCH --output=logs/SNP_processing_%j.out
#SBATCH --error=logs/SNP_processing_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=03:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

ml Bioinformatics bcftools
cd 4_genotyping_QC/3_GRM_validation

#======================================================== SNP processing =====================================================================
#------------------------Step1: Only retain bi-allelic SNP----------------------------------
vcftools --vcf high-quality_and_genotyped_cells.vcf \
         --min-alleles 2 --max-alleles 2 \
         --remove-indels --recode --recode-INFO-all \
         --out high-quality_and_genotyped_cells_bi-allelic-SNP


#-------------Step2: Simplify VCF: set ID as chr_pos, clear INFO field, keep only GT, and set genotype to ./. if read depth < 2 -----------------------
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    $3 = $1 "_" $2;
    $8 = ".";
    $9 = "GT";
    for (i = 10; i <= NF; i++) {
        split($i, a, ":");
        split(a[2], dp, ",");
        if (dp[1] < 2) a[1] = "./.";
        $i = a[1];
    }
    print;
}' high-quality_and_genotyped_cells_bi-allelic-SNP.recode.vcf > high-quality_and_genotyped_cells_bi-allelic-SNP.vcf

rm high-quality_and_genotyped_cells_bi-allelic-SNP.recode.vcf


#------------Step3: Identifies and merges SNPs shared between scATAC-seq and reference data within accessible chromatin regions---------------------

# 1.Define output files
MERGED_BED=merged.peak
scATAC_VCF=high-quality_and_genotyped_cells_bi-allelic-SNP.vcf     # Genotypes called by freebayes using BAMs from high-quality cells (output from Step 2)
ref_VCF=3parents_hybrids.vcf              # Reference genotypes of B73, Ki3, and Oh43 inbred lines and hybrids  

# 2.Merge the ATAC peak files and combine overlapping regions
cat ./ATAC-peak/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > $MERGED_BED

# 3.Only extract SNPs located in ACRs
bedtools intersect -a $scATAC_VCF -b $MERGED_BED -header > scATAC_ACR-localized.vcf
bedtools intersect -a $ref_VCF -b $MERGED_BED -header > ref_ACR-localized.vcf

# 4.Compress and index the filtered VCFs
bgzip scATAC_ACR-localized.vcf
bgzip ref_ACR-localized.vcf
tabix -p vcf scATAC_ACR-localized.vcf.gz
tabix -p vcf ref_ACR-localized.vcf.gz

# 5.Find SNPs shared between scATAC and reference VCFs in ACRs
bcftools isec -n=2 -w1 -Oz -o shared_between_scATAC_and_ref.vcf.gz scATAC_ACR-localized.vcf.gz ref_ACR-localized.vcf.gz

# 6.Extract the shared SNPs from both filtered VCFs
bcftools view -R shared_between_scATAC_and_ref.vcf.gz -Oz -o scATAC_intersect_with_ref_ACR-localized.vcf.gz scATAC_ACR-localized.vcf.gz
bcftools view -R shared_between_scATAC_and_ref.vcf.gz -Oz -o ref_intersect_with_scATAC_ACR-localized.vcf.gz ref_ACR-localized.vcf.gz
tabix  -p vcf scATAC_intersect_with_ref_ACR-localized.vcf.gz
tabix  -p vcf ref_intersect_with_scATAC_ACR-localized.vcf.gz

# 7.Set the INFO field to '.'
bcftools annotate -x INFO ref_intersect_with_scATAC_ACR-localized.vcf.gz -Ov -o Final_ref_intersect_with_scATAC_ACR-localized.vcf.gz
bcftools annotate -x INFO scATAC_intersect_with_ref_ACR-localized.vcf.gz -Ov -o Final_scATAC_intersect_with_ref_ACR-localized.vcf.gz
tabix  -p vcf Final_ref_intersect_with_scATAC_ACR-localized.vcf.gz
tabix  -p vcf Final_scATAC_intersect_with_ref_ACR-localized.vcf.gz

# 7.Merge the reference and scATAC SNPs at shared ACR-localized sites
bcftools merge -Oz -o Final_Ref_scATAC_merged_SNP.vcf.gz Final_ref_intersect_with_scATAC_ACR-localized.vcf.gz Final_scATAC_intersect_with_ref_ACR-localized.vcf.gz

# 8.Remove intermediate files
rm merged.peak ref* scATAC* shared* *intersect* high-quality_and_genotyped_cells_bi-allelic-SNP.vcf

echo "SNP processing completed. Proceed to compute the genetic relationship matrix (GRM)."
