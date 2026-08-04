#!/home/lgjiang/miniconda3/envs/myenv/bin/Rscript
# This script is executed using the R environment from the conda environment "myenv"
# USAGE: /home/lgjiang/miniconda3/envs/myenv/bin/Rscript step02_dACR_parent_and_hybrid_DESeq2.R
# Purpose:
#   Identify dACRs (DESeq2) in three contrasts, per cell type:
#     1) P1 vs P2
#     2) P1 vs Hybrid (Hybrid = BxK + KxB pooled as 4 reps)
#     3) P2 vs Hybrid (Hybrid = BxK + KxB pooled as 4 reps)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(DESeq2)
library(edgeR)

# ---------------- settings ----------------
in_tsv  <- "1_filtered_ACR/Final_raw_count_for_all_nonPAV-peaks_and_all_samples.annot.noExonic.tsv"
out_dir <- "3_dACRs_between_parents_and_hybrids"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

P1 <- "B73"
P2 <- "Ki3"
HYB1 <- "BxK"
HYB2 <- "KxB"

ann_n <- 8
padj_cut <- 0.05
min_row_sum <- 4

# SPA thresholds on CPM (single-parent ACR), only for P1 vs P2
spa_on  <- 1
spa_off <- 0.1

# ---------------- Read file ----------------
df <- read_tsv(in_tsv, show_col_types = FALSE)
ann_cols    <- colnames(df)[1:ann_n]
sample_cols <- colnames(df)[(ann_n + 1):ncol(df)]

# =============================================================
# ------------- CPM normalization (all samples) --------------
# =============================================================
count_mat_all <- as.matrix(df[, (ann_n + 1):ncol(df)])
storage.mode(count_mat_all) <- "integer"
cpm_mat_all <- edgeR::cpm(count_mat_all, log = FALSE)

df_cpm  <- bind_cols(df[, 1:ann_n], as_tibble(cpm_mat_all))
out_cpm <- file.path(out_dir, "Final_CPM_for_all_nonPAV-peaks_and_all_samples.annot.noExonic.tsv")
write_tsv(df_cpm, out_cpm)

# ==========================================================
# ------------- Annotation table (ensure unique peak_id) ----
# ==========================================================
ann_df_all <- df %>%
  select(all_of(ann_cols)) %>%
  distinct(peak_id, .keep_all = TRUE)

# ============================================
# ------------- Sample metadata --------------
# ============================================
meta <- tibble(sample = sample_cols) %>%
  mutate(
    cell_type   = str_extract(sample, "(?<=-)C\\d+"),
    prefix      = str_extract(sample, "^[^-]+"),
    prefix_core = str_extract(prefix, "^[^_]+"),
    group = case_when(
      prefix_core == P1 ~ P1,
      prefix_core == P2 ~ P2,
      prefix_core %in% c(HYB1, HYB2) ~ "Hybrid",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cell_type), !is.na(group))

cell_types <- meta %>% distinct(cell_type) %>% arrange(cell_type) %>% pull(cell_type)

# ================================================================================
# ------------- Define function: run DESeq2 contrast per cell type --------------
# ================================================================================
run_one_ct_contrast_wide <- function(ct, g1, g2, tag) {

  cols_ct <- meta %>%
    filter(cell_type == ct, group %in% c(g1, g2)) %>%
    arrange(group, sample)

  rep_n <- dplyr::count(cols_ct, group)
  if (nrow(rep_n) < 2 || any(rep_n$n < 2)) {
    message(sprintf("[SKIP] %s | %s : insufficient replicates (%s)",
                    ct, tag, paste(rep_n$group, rep_n$n, sep="=", collapse=", ")))
    return(NULL)
  }

  # Subset counts
  mat <- df %>% select(all_of(ann_cols), all_of(cols_ct$sample))
  count_mat <- as.matrix(mat %>% select(all_of(cols_ct$sample)))
  storage.mode(count_mat) <- "integer"
  rownames(count_mat) <- mat$peak_id

  # Prefilter
  keep <- rowSums(count_mat) >= min_row_sum
  count_mat <- count_mat[keep, , drop = FALSE]

  # colData
  coldata <- data.frame(
    group = factor(cols_ct$group[match(colnames(count_mat), cols_ct$sample)],
                   levels = c(g1, g2))
  )
  rownames(coldata) <- colnames(count_mat)

  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", g1, g2))

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("peak_id") %>%
    mutate(cell_type = ct)

  # mean CPM per group
  cpm_mat <- edgeR::cpm(count_mat, log = FALSE)
  g1_cols <- rownames(coldata)[coldata$group == g1]
  g2_cols <- rownames(coldata)[coldata$group == g2]
  meanCPM_G1 <- rowMeans(cpm_mat[, g1_cols, drop = FALSE])
  meanCPM_G2 <- rowMeans(cpm_mat[, g2_cols, drop = FALSE])

  # SPA only for parent-parent
  do_SPA <- (g1 %in% c(P1, P2) && g2 %in% c(P1, P2))

  out <- res_df %>%
    mutate(
      log2FC = log2FoldChange,
      FC_abs = ifelse(is.na(log2FC), NA_real_, 2^abs(log2FC)),

      meanCPM_G1 = meanCPM_G1,
      meanCPM_G2 = meanCPM_G2,

      is_SPA_G1 = if (do_SPA) (meanCPM_G1 >= spa_on & meanCPM_G2 < spa_off) else FALSE,
      is_SPA_G2 = if (do_SPA) (meanCPM_G2 >= spa_on & meanCPM_G1 < spa_off) else FALSE,

      direction = case_when(
        is.na(log2FC) ~ NA_character_,
        log2FC > 0 ~ paste0(g1, "_higher"),
        log2FC < 0 ~ paste0(g2, "_higher"),
        TRUE ~ "tie"
      ),

      dACR_class = case_when(
        is.na(padj) | padj >= padj_cut ~ "non-dACR",
        is_SPA_G1 ~ paste0("SPA_", g1),
        is_SPA_G2 ~ paste0("SPA_", g2),
        FC_abs >= 4 ~ "dACR_FC>4",
        FC_abs >= 2 ~ "dACR_FC2-4",
        TRUE ~ "non-dACR"
      )
    ) %>%
    select(
      peak_id, cell_type,
      log2FC, pvalue, padj, direction, dACR_class,
      meanCPM_G1, meanCPM_G2
    )

  # ---- rename CPM columns without := (compatible with older dplyr/rlang) ----
  colnames(out)[colnames(out) == "meanCPM_G1"] <- paste0("meanCPM_", g1)
  colnames(out)[colnames(out) == "meanCPM_G2"] <- paste0("meanCPM_", g2)

  # ---- add tag suffix to all non-key columns ----
  stat_cols_now <- setdiff(colnames(out), c("peak_id", "cell_type"))
  colnames(out)[match(stat_cols_now, colnames(out))] <- paste0(stat_cols_now, "__", tag)

  return(out)
}

# =========================================================================
# ------------- For each cell type: UNION across 3 contrasts --------------
# =========================================================================
merge_three_contrasts_one_ct <- function(ct) {
  tag_pp  <- paste0(P1, "_vs_", P2)
  tag_p1h <- paste0(P1, "_vs_Hybrid")
  tag_p2h <- paste0(P2, "_vs_Hybrid")

  r_pp  <- run_one_ct_contrast_wide(ct, P1, P2,       tag = tag_pp)
  r_p1h <- run_one_ct_contrast_wide(ct, P1, "Hybrid", tag = tag_p1h)
  r_p2h <- run_one_ct_contrast_wide(ct, P2, "Hybrid", tag = tag_p2h)

  key <- c("peak_id", "cell_type")
  out <- NULL

  if (!is.null(r_pp))  out <- r_pp
  if (!is.null(r_p1h)) out <- if (is.null(out)) r_p1h else full_join(out, r_p1h, by = key)
  if (!is.null(r_p2h)) out <- if (is.null(out)) r_p2h else full_join(out, r_p2h, by = key)

  return(out)
}

# ==========================================================
# ------------- Output file per cell type ----
# ==========================================================
tag_pp  <- paste0(P1, "_vs_", P2)
tag_p1h <- paste0(P1, "_vs_Hybrid")
tag_p2h <- paste0(P2, "_vs_Hybrid")

stat_cols_template <- c(
  paste0("log2FC__", tag_pp),    paste0("pvalue__", tag_pp),    paste0("padj__", tag_pp),
  paste0("direction__", tag_pp), paste0("dACR_class__", tag_pp),
  paste0("meanCPM_", P1, "__", tag_pp), paste0("meanCPM_", P2, "__", tag_pp),

  paste0("log2FC__", tag_p1h),    paste0("pvalue__", tag_p1h),    paste0("padj__", tag_p1h),
  paste0("direction__", tag_p1h), paste0("dACR_class__", tag_p1h),
  paste0("meanCPM_", P1, "__", tag_p1h), paste0("meanCPM_Hybrid__", tag_p1h),

  paste0("log2FC__", tag_p2h),    paste0("pvalue__", tag_p2h),    paste0("padj__", tag_p2h),
  paste0("direction__", tag_p2h), paste0("dACR_class__", tag_p2h),
  paste0("meanCPM_", P2, "__", tag_p2h), paste0("meanCPM_Hybrid__", tag_p2h)
)

for (ct in cell_types) {
  message(sprintf("[RUN] Processing %s ...", ct))

  ct_base <- ann_df_all %>%
    mutate(cell_type = ct) %>%
    select(all_of(ann_cols), cell_type)

  ct_core <- merge_three_contrasts_one_ct(ct)

  if (!is.null(ct_core) && nrow(ct_core) > 0) {
    ct_wide <- ct_base %>%
      left_join(ct_core, by = c("peak_id", "cell_type"))
  } else {
    ct_wide <- ct_base
  }

  missing_cols <- setdiff(stat_cols_template, colnames(ct_wide))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) ct_wide[[mc]] <- NA
  }

  ct_wide <- ct_wide %>%
    select(all_of(ann_cols), cell_type, all_of(stat_cols_template))

  out_file <- file.path(out_dir, paste0("dACR_DESeq2_", ct, ".tsv"))
  write_tsv(ct_wide, out_file)

  message(sprintf("[DONE] %s -> %s (rows=%d, cols=%d)",
                  ct, out_file, nrow(ct_wide), ncol(ct_wide)))
}


