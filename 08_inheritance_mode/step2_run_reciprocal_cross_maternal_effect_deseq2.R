# This script is used to perform DESeq2-based analysis to assess maternal effects in reciprocal cross.

library(DESeq2)
library(readr)
library(dplyr)
library(stringr)
library(tibble)

in_tsv  <- "1_filtered_ACR/Final_raw_count_for_all_nonPAV-peaks_and_all_samples.annot.noExonic.tsv"
out_dir <- "2_maternal_effects_evaluation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_total <- 4                 # Minimum total read count across all 4 samples required for an ACR to be tested
fdr_cut   <- 0.05              # FDR threshold for significance 
lfc_cut   <- 1                 # Absolute log2 fold-change threshold for calling a maternal effect

# load input file
df <- read_tsv(in_tsv, show_col_types = FALSE) %>%
  mutate(row_id = row_number())

ann_cols  <- colnames(df)[1:8]             # Annotation columns
count_cols <- colnames(df)[9:ncol(df)]     # Raw count columns


# ---- parse count column names: <geno>_rep<rep>-<celltype> ----
meta <- tibble(col = count_cols) %>%
  mutate(
    geno = str_match(col, "^(.+)_rep[0-9]+-C[0-9]+$")[,2],
    rep  = str_match(col, "_rep([0-9]+)-C")[,2],
    cell = str_match(col, "-(C[0-9]+)$")[,2]
  ) %>%
  filter(!is.na(geno) & !is.na(rep) & !is.na(cell)) %>%
  mutate(rep = as.integer(rep))

# only keep rep1/rep2 (because we only have 2 replicates for each genotype and cell type)
meta <- meta %>% filter(rep %in% c(1,2))

# ---- find reciprocal genotype pairs: "A x B" vs "B x A" ----
# For a genotype string like "BxK", treat it as "B" and "K" around 'x'
geno_pairs <- meta %>%
  distinct(geno) %>%
  mutate(
    A = str_match(geno, "^(.+)x(.+)$")[,2],
    B = str_match(geno, "^(.+)x(.+)$")[,3]
  ) %>%
  filter(!is.na(A) & !is.na(B)) %>%
  mutate(key = ifelse(A < B, paste0(A, "|", B), paste0(B, "|", A))) %>%
  group_by(key) %>%
  summarise(genos = list(sort(geno)), .groups = "drop") %>%
  filter(lengths(genos) == 2)

cells <- sort(unique(meta$cell))


# ==========================================================================================
# ---- function: run DESeq2 for one pair + one cell type (full rows, NA for low counts) ----
# ==========================================================================================
run_one <- function(g1, g2, ct) {

  cols_ct <- meta %>%
    filter(geno %in% c(g1, g2), cell == ct, rep %in% c(1,2)) %>%
    arrange(geno, rep)

  # make sure ordering is g1 rep1/2 then g2 rep1/2
  cols_order <- c(
    cols_ct %>% filter(geno == g1) %>% arrange(rep) %>% pull(col),
    cols_ct %>% filter(geno == g2) %>% arrange(rep) %>% pull(col)
  )
  if (length(cols_order) != 4) return(NULL)

  counts_full <- df %>% select(all_of(cols_order)) %>% as.data.frame()
  counts_full[] <- lapply(counts_full, function(x) as.integer(round(x)))
  total <- rowSums(counts_full)
  testable <- total >= min_total

  res_all <- tibble(
    row_id = df$row_id,
    baseMean = NA_real_, log2FoldChange = NA_real_, lfcSE = NA_real_,
    stat = NA_real_, pvalue = NA_real_, padj = NA_real_
  )

  if (sum(testable) > 0) {
    counts_sub <- counts_full[testable, , drop = FALSE]
    rownames(counts_sub) <- df$row_id[testable]
    colnames(counts_sub) <- c("g1_r1","g1_r2","g2_r1","g2_r2")

    coldata <- data.frame(
      cross = factor(c(g1, g1, g2, g2), levels = c(g1, g2))
    )
    rownames(coldata) <- colnames(counts_sub)

    dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = coldata, design = ~ cross)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("cross", g2, g1))  # log2FC = g2 / g1

    res_df <- as.data.frame(res) %>%
      rownames_to_column("row_id") %>%
      mutate(row_id = as.integer(row_id)) %>%
      select(row_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

    res_all <- res_all %>%
      select(row_id) %>%
      left_join(res_df, by = "row_id")
  }

  out_df <- df %>%
    select(row_id, all_of(ann_cols)) %>%
    left_join(res_all, by = "row_id") %>%
    bind_cols(df %>% select(all_of(cols_order))) %>%
    mutate(
      pair = paste0(g1, "_vs_", g2),
      cell_type = ct,
      tested = ifelse(total >= min_total, "yes", "no"),
      maternal_effect = ifelse(!is.na(padj) & padj < fdr_cut & abs(log2FoldChange) > lfc_cut, "yes", "no"),
      direction = case_when(
        maternal_effect == "yes" & log2FoldChange > 0 ~ paste0(g2, "_higher"),
        maternal_effect == "yes" & log2FoldChange < 0 ~ paste0(g1, "_higher"),
        TRUE ~ NA_character_
      )
    )

  out_file <- file.path(out_dir, paste0("DESeq2_", g1, "_vs_", g2, "_", ct, ".tsv"))
  write_tsv(out_df, out_file)

  tibble(
    pair = paste0(g1, "_vs_", g2),
    cell_type = ct,
    n_peaks = nrow(out_df),
    n_tested = sum(out_df$tested == "yes"),
    n_sig = sum(out_df$maternal_effect == "yes", na.rm = TRUE)
  )
}


# ===================================================================
# ---- main loop over each reciprocal cross pair and cell type ----
# ===================================================================
summary_list <- list()

for (i in seq_len(nrow(geno_pairs))) {
  gset <- geno_pairs$genos[[i]]
  g1 <- gset[1]
  g2 <- gset[2]

  for (ct in cells) {
    s <- run_one(g1, g2, ct)
    if (!is.null(s)) summary_list[[length(summary_list)+1]] <- s
  }
}

summary_df <- bind_rows(summary_list) %>%
  mutate(frac_sig = ifelse(n_tested > 0, n_sig / n_tested, NA_real_))

write_tsv(summary_df, file.path(out_dir, "summary_by_pair_and_celltype.tsv"))

