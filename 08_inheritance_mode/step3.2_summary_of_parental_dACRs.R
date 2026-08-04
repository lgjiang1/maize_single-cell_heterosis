library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
in_dir  <- "3_dACRs_between_parents_and_hybrids"

# (1) per-cell-type summary (all 14 cell types)
out_tsv <- file.path(in_dir, "dACR_parental_B-K_summary.tsv")

# (2) feature-level direction consistency
# Note: Cell types C13 and C14 were excluded from this analysis because they contain substantially fewer nuclei, resulting in very few
# detected ACRs and consequently very few parental dACRs (<10). Including these cell types would introduce excessive sparsity 
# and instability in cross–cell-type consistency estimates.
exclude_ct <- c("C13", "C14")
total_ct   <- 12
out_feature_tsv <- file.path(in_dir, "dACR_parental_B-K_direction_consistency_C1-12.per_feature.tsv")

# List all per-cell-type DESeq2 result files
files <- list.files(in_dir, pattern = "^dACR_DESeq2_C[0-9]+\\.tsv$", full.names = TRUE)


# ------------------------------------------
# PART 1: per-cell-type summary
# ------------------------------------------
wide_one <- function(f){
  ct <- str_extract(basename(f), "C[0-9]+")
  df <- read_tsv(f, show_col_types = FALSE)
  
  # only extract comparison between two parents for following analysis
  dir_col   <- names(df)[13]
  class_col <- names(df)[14]

  sub <- df %>%
    filter(.data[[class_col]] %in% c("dACR_FC>4","dACR_FC2-4") |
             str_detect(.data[[class_col]], "^SPA_")) %>%
    mutate(category = recode(.data[[class_col]],
        "dACR_FC>4"  = "dACR_FC_gt4",
        "dACR_FC2-4" = "dACR_FC_2to4",
        .default = as.character(.data[[class_col]])),
        direction = as.character(.data[[dir_col]]),
        key = paste0(category, "__", direction)
    )

  w1 <- sub %>%
    count(key) %>%
    pivot_wider(names_from = key, values_from = n, values_fill = 0)

  w2 <- sub %>%
    count(direction) %>%
    mutate(key = paste0("dACR_total__", direction)) %>%
    select(key, n) %>%
    pivot_wider(names_from = key, values_from = n, values_fill = 0)
  bind_cols(tibble(cell_type = ct), w1, w2)
}

wide <- map_dfr(files, wide_one) %>%
  arrange(as.integer(str_remove(cell_type, "^C")))

write_tsv(wide, out_tsv)


# -------------------------------------------------
# PART 2: feature-level direction consistency
# -------------------------------------------------
long_one <- function(f){
  ct <- str_extract(basename(f), "C[0-9]+")
  if (is.na(ct) || ct %in% exclude_ct) return(NULL)

  df <- read_tsv(f, show_col_types = FALSE)

  dir_col   <- names(df)[13]
  class_col <- names(df)[14]

  df %>%
    filter(.data[[class_col]] %in% c("dACR_FC>4","dACR_FC2-4") |
             str_detect(.data[[class_col]], "^SPA_")) %>%
    transmute(
      cell_type = ct,
      peak_id   = as.character(.data[["peak_id"]]),
      direction = as.character(.data[[dir_col]]),
      is_SPA    = str_detect(.data[[class_col]], "^SPA_")
    ) %>%
    distinct()
}

long_all <- map_dfr(files, long_one)

per_feature <- long_all %>%
  group_by(peak_id) %>%
  summarise(
    is_SPA = any(is_SPA),
    n_ct   = n_distinct(cell_type),
    n_dir  = n_distinct(direction),
    dirs   = paste(sort(unique(direction)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    consistency = if_else(n_dir == 1, "consistent", "mixed"),
    dir_consistent = if_else(n_dir == 1, dirs, NA_character_),
    bin_nct = case_when(
      n_ct >= 1  & n_ct <= 3  ~ "1~3",
      n_ct >= 4  & n_ct <= 9  ~ "4~9",
      n_ct >= 10 & n_ct <= 12 ~ "10~12"
    )
  )

write_tsv(per_feature, out_feature_tsv)

