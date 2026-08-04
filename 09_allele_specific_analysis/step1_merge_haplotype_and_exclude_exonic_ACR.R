library(data.table)
library(stringr)

hap1_name <- "B73"
hap2_name <- "Ki3"
in_dir <- "0_input"
out_dir <- "1_merged"
file1 <- file.path(in_dir, paste0("nonPAV_counts_all_samples_", hap1_name, "ref.tsv"))
file2 <- file.path(in_dir, paste0("nonPAV_counts_all_samples_", hap2_name, "ref.tsv"))
filter_file <- file.path(in_dir, "Final_raw_count_for_all_nonPAV-peaks_and_all_samples.annot.noExonic.tsv")

out_file <- file.path(out_dir, paste0("merged_", hap1_name, "_", hap2_name, "_counts_without_exonic_ACR.tsv"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read input count tables
dt1 <- fread(file1)
dt2 <- fread(file2)

dt <- merge(dt1, dt2, by = "peak_id")
filter_ids <- unique(fread(filter_file, select = "peak_id")$peak_id)
dt <- dt[peak_id %in% filter_ids]

hap1_coord <- paste0(hap1_name, c("_chr", "_start", "_end"))
hap2_coord <- paste0(hap2_name, c("_chr", "_start", "_end"))

count1 <- setdiff(names(dt1), c(hap1_coord, "peak_id"))
count2 <- setdiff(names(dt2), c(hap2_coord, "peak_id"))

get_order <- function(cols) {
  info <- data.table(
    col = cols,
    cell = as.integer(str_remove(str_extract(cols, "C\\d+"), "C")),
    rep  = as.integer(str_remove(str_extract(cols, "rep\\d+"), "rep"))
  )
  info[order(cell, rep, col), col]
}

final_cols <- c("peak_id", hap1_coord, hap2_coord, get_order(count1), get_order(count2))

fwrite(dt[, ..final_cols], out_file, sep = "\t", quote = FALSE)
