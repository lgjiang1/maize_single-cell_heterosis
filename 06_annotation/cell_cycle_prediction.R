# ---------------------------
#    cell cycle prediction
# ---------------------------
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript cell_cycle_prediction.R <gene_activity_rds> <marker_file_txt> [output_file]")
}
gene_activity_file <- args[1]
marker_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "cell_cycle_prediction_result.txt")

# ================ 1. Load data ==========================
gene_activity <- readRDS(gene_activity_file) 
marker_df <- read.table(marker_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)


# ================ 2. Build marker lists by cell_cycle_phase =========================
marker_list_raw <- split(marker_df$Gene_v5, marker_df$Cell_cycle_phase)
# Keep only markers present in the gene activity matrix
overlap_list <- lapply(marker_list_raw, function(v) intersect(v, rownames(gene_activity)))

cat("Number of overlapping markers per phase:\n")
print(sapply(overlap_list, length))

# Create the union of all overlapping markers (cell cycle gene pool for following permutation)
cellcycle_pool <- unique(unlist(overlap_list))


# ================ 3. Compute real cell cycle scores: colMeans per phase =================
cell_names <- colnames(gene_activity)
# Compute mean gene activity per phase for each cell
real_scores_list <- lapply(names(overlap_list), function(stage) {
  genes <- overlap_list[[stage]]
  if (length(genes) == 0) {
    setNames(rep(NA_real_, length(cell_names)), cell_names)
  } else {
    Matrix::colMeans(gene_activity[genes, , drop = FALSE], na.rm = TRUE)
  }
})
names(real_scores_list) <- names(overlap_list)
real_scores <- do.call(cbind, real_scores_list)
colnames(real_scores) <- names(overlap_list)
rownames(real_scores) <- cell_names
real_scores <- as.data.frame(real_scores, check.names = FALSE)

# ---- Identify G0 (non-cycling): all phases are 0, and at least one non-NA exists ----
is_G0 <- apply(real_scores, 1, function(x) all(x == 0, na.rm = TRUE) & any(!is.na(x)))
cat("Non-cycling (G0) cells detected:", sum(is_G0), "of", length(is_G0), "cells\n")


# ================ 4. Generate permutation scores (skip G0 cells)  ===========
set.seed(2025)
n_perm <- 1000
all_genes <- rownames(gene_activity)
cells_use <- cell_names[!is_G0]  # only cycling candidates

# For each phase: randomly sample genes (same count as overlapping markers) and compute scores
perm_scores_list <- lapply(names(overlap_list), function(stage) {
  n <- length(overlap_list[[stage]])
  if (n == 0 || length(cells_use) == 0) {
    return(matrix(NA_real_, nrow = length(cells_use), ncol = n_perm,
                  dimnames = list(cells_use, NULL)))
  }
  replicate(n_perm, {
    random_genes <- sample(cellcycle_pool, n)
    cm <- Matrix::colMeans(gene_activity[random_genes, , drop = FALSE], na.rm = TRUE)
    cm[cells_use]
  })
})
names(perm_scores_list) <- names(overlap_list)


# ================ 5. Calculate Z-scores for cycling cells =======================
z_scores <- matrix(NA_real_, nrow = length(cell_names), ncol = length(overlap_list),
                   dimnames = list(cell_names, names(overlap_list)))

if (length(cells_use) > 0) {
  for (stage in names(overlap_list)) {
    real_vec <- real_scores[cells_use, stage]
    perm_mat <- perm_scores_list[[stage]] 
    perm_means <- rowMeans(perm_mat, na.rm = TRUE)
    perm_sds   <- apply(perm_mat, 1, sd, na.rm = TRUE)
    z_vec <- (real_vec - perm_means) / perm_sds
    z_vec[!is.finite(z_vec)] <- NA
    z_scores[cells_use, stage] <- z_vec
  }
}

colnames(z_scores) <- paste0("Z_", colnames(z_scores))
z_scores <- as.data.frame(z_scores, check.names = FALSE)


# =============== 6. Assign stage: max Z for cycling; G0 stays G0 ================
predicted_stage <- rep(NA_character_, length(cell_names))
names(predicted_stage) <- cell_names

# For cycling cells: assign the stage with the highest Z-score
if (length(cells_use) > 0) {
  cycling_max <- apply(z_scores[cells_use, , drop = FALSE], 1, function(x) {
    if (all(is.na(x))) return(NA_character_)
    colnames(z_scores)[which.max(x)]
  })
  # Remove the "Z_" prefix to get the stage name
  predicted_stage[cells_use] <- sub("^Z_", "", cycling_max)
}
# Assign G0 and Unknown
predicted_stage[is_G0] <- "G0"
predicted_stage[is.na(predicted_stage)] <- "Unknown"


# =============== 7. Save final results ======================
result <- data.frame(
  cell = cell_names,
  real_scores,
  z_scores,
  predicted_stage = predicted_stage,
  is_G0 = is_G0,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(result, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Cell cycle prediction done! Results saved to:", output_file, "\n")


