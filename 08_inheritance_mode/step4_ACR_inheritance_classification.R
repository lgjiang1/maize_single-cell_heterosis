library(readr)
library(dplyr)
library(purrr)
library(stringr)

# ---------------- User settings ----------------
cpm_file  <- "3_dACRs_between_parents_and_hybrids/Final_CPM_for_all_nonPAV-peaks_and_all_samples.annot.noExonic.tsv"
deseq_dir <- "3_dACRs_between_parents_and_hybrids"   # dACR_DESeq2_C*.tsv
out_dir   <- "4_inheritance_mode"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
summary_file <- file.path(out_dir, "B-K_inheritance_mode_summary.tsv")  # change B-K if needed

# Step1 thresholds (DESeq2-based gating)
padj_cut <- 0.05
lfc_cut  <- 1   # |log2FC| < 1  ==> FC < 2

# Step2 thresholds (fold-based inheritance classification)
fc_cut <- 1.5   # 1.5-fold threshold (OD/UD + additive window around MPV)

# -----------------------------------------------------------------
# Define Step1 functions (identify constant and non-constant ACRs)
# -----------------------------------------------------------------
# 1. One DESeq2 comparison is "constant-like"
# If log2FC or padj is NA (often due to all-zero counts / not estimable), treat as "pass" (non-significant by default).
is_constant_comp <- function(log2fc, padj, lfc_cut, padj_cut) {
  pass_na  <- is.na(log2fc) | is.na(padj)
  # A peak is classified as constant only if all three comparisons show |log2FC| < lfc_cut and padj > padj_cut.
  pass_val <- is.finite(log2fc) & (abs(log2fc) < lfc_cut) & (padj > padj_cut)
  pass_na | pass_val
}

# 2. Discover paired (log2FC__, padj__) columns automatically
# For our data, we have 3 comparisons: P1-P2; P1-hybird; P2hybird
discover_comparisons <- function(cols) {
  lfc_cols  <- cols[str_detect(cols, "^log2FC__")]
  padj_cols <- cols[str_detect(cols, "^padj__")]
  lfc_tags  <- str_replace(lfc_cols,  "^log2FC__", "")
  padj_tags <- str_replace(padj_cols, "^padj__",  "")
  common_tags <- intersect(lfc_tags, padj_tags)
  tibble(tag = common_tags, lfc_col  = paste0("log2FC__", common_tags), padj_col = paste0("padj__",  common_tags))
}

# 3. Parse parents from CPM header:
# col1: B73_chr -> "B73"; col6: Ki3_chr -> "Ki3"
get_parents_from_header <- function(df) {
  cn <- colnames(df)
  p1 <- str_replace(cn[1], "_chr$", "")
  p2 <- str_replace(cn[6], "_chr$", "")
  list(P1 = p1, P2 = p2)
}

# 4. Subset CPM columns for the current cell type: keep first 8 annotation cols + exactly 8 sample cols ending with "-C#"
subset_cpm_for_ct <- function(cpm_df, ct) {
  ann_n <- 8
  ann_cols <- colnames(cpm_df)[1:ann_n]

  # Match columns ending with "-C1", "-C2", ...
  ct_pattern <- paste0("-", ct, "$")
  ct_cols <- colnames(cpm_df)[str_detect(colnames(cpm_df), ct_pattern)]
  cpm_df %>% select(all_of(c(ann_cols, sort(ct_cols))))
}

# 5. Sort DESeq2 files by cell index if present (C1, C2, ...)
extract_c_index <- function(x) {
  m <- str_match(basename(x), "^dACR_DESeq2_C(\\d+)\\.tsv$")
  ifelse(is.na(m[, 2]), Inf, as.integer(m[, 2]))
}


# -------------------------------------------------------------------------
# Define Step2 functions (classify inheritance modes for non-constant ACRs)
# -------------------------------------------------------------------------
# 1. Identify 2 parental replicate columns each + 4 hybrid columns inside a Step1 output table
get_sample_columns <- function(df, p1, p2) {
  cn <- colnames(df)
  ann_cols <- cn[1:8]
  sample_cols <- setdiff(cn, c(ann_cols, "constancy"))

  # Parent columns: "B73_rep1-C10" etc.
  p1_cols <- sample_cols[str_detect(sample_cols, paste0("^", p1, "_rep[0-9]+-C"))]
  p2_cols <- sample_cols[str_detect(sample_cols, paste0("^", p2, "_rep[0-9]+-C"))]
  hy_cols <- setdiff(sample_cols, c(p1_cols, p2_cols))

  list(
    p1_cols = p1_cols[order(p1_cols)],
    p2_cols = p2_cols[order(p2_cols)],
    hy_cols = hy_cols[order(hy_cols)]
  )
}

# 2. calculate average between replicates
row_mean <- function(df, cols) {
  mat <- as.matrix(df[, cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rowMeans(mat, na.rm = TRUE)
}


# 3. Mutually exclusive, strict order:
# 1) OD: F1 > 1.5 * HP
# 2) UD: F1 < LP / 1.5
# 3) Additive: MPV/1.5 <= F1 <= MPV*1.5
# 4) Dominance: else -> closer to which parent (P1-dominance / P2-dominance)

classify_inheritance <- function(P1, P2, F1, p1_name, p2_name, fc_cut = 1.5) {
  HP  <- pmax(P1, P2)
  LP  <- pmin(P1, P2)
  MPV <- (P1 + P2) / 2

  od  <- F1 > fc_cut * HP
  ud  <- F1 < (LP / fc_cut)
  add <- (!od) & (!ud) & (F1 >= (MPV / fc_cut)) & (F1 <= (MPV * fc_cut))

  closer_p1 <- abs(F1 - P1) <= abs(F1 - P2)

  ifelse(
    od,  "overdominance",
    ifelse(
      ud,  "underdominance",
      ifelse(
        add, "additive",
        ifelse(
          closer_p1,
          paste0(p1_name, "-dominance"),
          paste0(p2_name, "-dominance")
        )
      )
    )
  )
}

# -----------------------------------------------------------------
# Run Step1 -> produce CPM_C*_with_constancy.tsv
# -----------------------------------------------------------------
cpm_all <- suppressMessages(read_tsv(cpm_file, show_col_types = FALSE)) # Read CPM file

# Find all DESeq2 files
deseq_paths <- list.files(path = deseq_dir, pattern = "^dACR_DESeq2_C.*\\.tsv$", full.names = TRUE)
deseq_paths <- deseq_paths[order(extract_c_index(deseq_paths), deseq_paths)]
message("Found ", length(deseq_paths), " DESeq2 files.")

for (path in deseq_paths) {
  fname <- basename(path)
  ct <- str_match(fname, "^dACR_DESeq2_(C[^\\.]+)\\.tsv$")[, 2]
  message("Step1: processing ", fname, " ...")

  deseq <- suppressMessages(read_tsv(path, show_col_types = FALSE))
  comps <- discover_comparisons(colnames(deseq))

  const_mat <- map2_dfc(comps$lfc_col, comps$padj_col,
    ~ is_constant_comp(deseq[[.x]], deseq[[.y]], lfc_cut, padj_cut))

  flag <- tibble(peak_id = deseq$peak_id, constancy = ifelse(apply(const_mat, 1, all), "constant", "non-constant"))

  # Subset CPM for this cell type (8 sample columns only)
  cpm_ct <- subset_cpm_for_ct(cpm_all, ct)

  # Add constancy
  step1_out <- cpm_ct %>%left_join(flag, by = "peak_id") %>%mutate(constancy = ifelse(is.na(constancy), "non-constant", constancy))
  step1_file <- file.path(out_dir, paste0("CPM_", ct, "_with_constancy.tsv"))
  write_tsv(step1_out, step1_file)


# -----------------------------------------------------------------
# Run Step2 -> produce CPM_C*_with_constancy_and_inheritance.tsv
# -----------------------------------------------------------------
message("Step2: classifying inheritance for ", ct, " ...")

parents <- get_parents_from_header(step1_out)
  p1 <- parents$P1
  p2 <- parents$P2

  cols <- get_sample_columns(step1_out, p1, p2)

  P1_mean <- row_mean(step1_out, cols$p1_cols)
  P2_mean <- row_mean(step1_out, cols$p2_cols)
  F1_mean <- row_mean(step1_out, cols$hy_cols)  # pooled reciprocals (4 cols)

  inh <- rep(NA_character_, nrow(step1_out))
  inh[step1_out$constancy == "constant"] <- "constant"

  idx <- which(step1_out$constancy == "non-constant")
  if (length(idx) > 0) {
    inh[idx] <- classify_inheritance(P1_mean[idx], P2_mean[idx], F1_mean[idx], p1_name = p1, p2_name = p2,fc_cut  = fc_cut)
  }

  step2_out <- step1_out %>% mutate(inheritance_mode = inh)
  step2_file <- file.path(out_dir, paste0("CPM_", ct, "_with_constancy_and_inheritance.tsv"))
  write_tsv(step2_out, step2_file)
}


# -----------------------------------------------------------------
# Step3 -> summarize inheritance_mode counts per cell type
# -----------------------------------------------------------------
message("Step3: summarizing inheritance modes across all cell types ...")

step2_paths <- list.files(path = out_dir, pattern = "^CPM_C.*_with_constancy_and_inheritance\\.tsv$", full.names = TRUE)
step2_paths <- step2_paths[order(extract_c_index(step2_paths), step2_paths)]

# parse C1, C2, ... from filename
get_ct_from_step2 <- function(path) {
  m <- str_match(basename(path),
                 "^CPM_(C[^_]+)_with_constancy_and_inheritance\\.tsv$")
  m[, 2]
}

# read + stack
all_step2 <- map_dfr(step2_paths, function(f) {
  x <- suppressMessages(read_tsv(f, show_col_types = FALSE))
  ct <- get_ct_from_step2(f)
  tibble(cell_type = ct, inheritance_mode = x$inheritance_mode)
})

# all observed modes (keeps B73-dominance / Ki3-dominance etc.)
mode_levels <- sort(unique(na.omit(all_step2$inheritance_mode)))

# summary
summary_df <- all_step2 %>%filter(!is.na(inheritance_mode)) %>%count(cell_type, inheritance_mode, name = "n") %>%
  tidyr::complete(cell_type, inheritance_mode = mode_levels, fill = list(n = 0)) %>%group_by(cell_type) %>%
  mutate(total = sum(n), constant_n = n[inheritance_mode == "constant"], nonconstant_n = total - constant_n,
  frac_total = ifelse(total > 0, n / total, NA_real_), constant_frac = ifelse(total > 0, constant_n / total, NA_real_),
  frac_nonconstant = ifelse(inheritance_mode == "constant", NA_real_, ifelse(nonconstant_n > 0, n / nonconstant_n, 0))) %>%ungroup()
write_tsv(summary_df, summary_file)

