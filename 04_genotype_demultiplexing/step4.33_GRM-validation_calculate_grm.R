setwd("4_genotyping_QC/3_GRM_validation")

library(SNPRelate)

# Step 1: Convert the VCF file to GDS format, keeping only biallelic SNPs
vcf.fn <- "Final_Ref_scATAC_merged_SNP.vcf.gz"
gds.fn <- "Final_Ref_scATAC_merged_SNP.gds"
snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")

# Step 2: Open the GDS file
genofile <- snpgdsOpen(gds.fn)

# Step 3: Remove SNPs with MAF<0.05 and missing rate>0
snp.set <- unlist(snpgdsSelectSNP(genofile, maf=0.05, missing.rate=0))

# Step 4: Calculate the Genetic Relationship Matrix (GRM) using LD-pruned SNPs
grm <- snpgdsGRM(genofile, snp.id = snp.set, method = "Corr")

# Step 5: Extract the GRM matrix and save to file
grm_mat <- grm$grm
rownames(grm_mat) <- grm$sample.id
colnames(grm_mat) <- grm$sample.id
write.table(grm_mat, file = "GRM_matrix.txt", quote = FALSE)


# ======================== PCA visualization ============================
# only keep 18 real samples, removing ref and pseudo genotypes
vcf.fn <- "Final_Ref_scATAC_merged_SNP_for_PCA.vcf.gz"
gds.fn <- "Final_Ref_scATAC_merged_SNP_for_PCA.gds"
snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")
genofile <- snpgdsOpen(gds.fn)

snp.set <- unlist(snpgdsSelectSNP(genofile, maf=0.05, missing.rate=0))

# get pca
pca <- snpgdsPCA(genofile, snp.id=snp.set, autosome.only=FALSE)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

df_pca <- data.frame(
  sample = sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2]
)

# save results
write.table(df_pca, file = "snp_pca_result.tab", sep = "\t", quote = FALSE, row.names = FALSE)

