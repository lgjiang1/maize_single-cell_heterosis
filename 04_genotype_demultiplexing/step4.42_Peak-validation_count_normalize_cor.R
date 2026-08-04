library(edgeR)
library(limma)

setwd("4_genotyping_QC/4_peak_validation/")

# Load the raw count matrix
count_mat <- read.table("all_ACR_counts.txt", header = TRUE)
count_mat <- count_mat[, 4:ncol(count_mat)]   # Remove 'chr', 'start', 'end'

# CPM normalization
dge <- DGEList(counts = count_mat)
cpm_mat <- cpm(dge, log=T, prior.count=5)

# Quantile normalization
qnorm_mat <- normalizeQuantiles(cpm_mat)

# Remove batch effect using limma
library.batch <- c(rep("published", 6), rep("LSC2", 12),rep("LSC3", 12))
qnorm_corrected <- removeBatchEffect(qnorm_mat, batch = library.batch)

# Compute Spearman correlation between samples
cor_mat <- cor(qnorm_corrected, method = "spearman")

# Save results
write.table(qnorm_corrected, file = "normalized_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cor_mat, file = "all_sample_spearman_correlation.txt", sep = "\t", quote = FALSE)

