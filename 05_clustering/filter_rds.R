# This script filters Socrates objects by keeping only high-quality genotyped cells; 
# input: a raw Socrates RDS file and a corresponding metadata (only genotyped cells) file.

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]  # e.g., "LSC2_BxK"

library(Socrates)
library(Matrix)

cat("Processing sample:", sample_name, "\n")

# ========= Get paths from environment variables ======
rds_dir  <- Sys.getenv("RDS_DIR")
meta_dir <- Sys.getenv("META_DIR")
out_dir  <- Sys.getenv("OUT_DIR")

# ============== Construct full paths =================
rds_file  <- file.path(rds_dir,  paste0(sample_name, ".raw.soc.rds"))
meta_file <- file.path(meta_dir, paste0(sample_name, "_barcode.txt"))

cat("Reading:\n", rds_file, "\n", meta_file, "\n")

# ============= Load data ===============
obj <- readRDS(rds_file)
meta_info <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(meta_info) <- meta_info$cellID

# ============== Only keep high-quality and genotyped cells =============
keep_cells <- intersect(colnames(obj$counts), rownames(meta_info))
if (length(keep_cells) == 0) stop("No overlapping cells for ", sample_name)

obj$counts <- obj$counts[, keep_cells]
obj$meta   <- obj$meta[keep_cells, ]

# === Add all extra columns from meta_info (including "library" and "genotype" information) ===
for (col in setdiff(colnames(meta_info), "cellID")) {
  obj$meta[[col]] <- meta_info[keep_cells, col]
}

# ============== Double-check cellID consistency ===================
stopifnot(identical(colnames(obj$counts), rownames(obj$meta)))

# ====================== Save filtered object ===================
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(obj, file.path(out_dir, paste0(sample_name, ".filtered.soc.rds")))

cat("Finished sample:", sample_name, "\n")

