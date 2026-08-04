# -------------------------------
# Usage:
# Rscript estimate_gene_activity.R [bed] [ann] [meta] [peaks] [outdir] 
# -------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript estimate_gene_activity.R [bed] [ann] [meta] [peaks] [outdir]")
}

bed     <- args[1]
ann     <- args[2]
meta    <- args[3]
peaks   <- args[4]
outdir  <- args[5]

# -------------------------------
# Load required libraries
# -------------------------------
library(Socrates)
library(GenomicFeatures)
library(Matrix)
library(data.table)

# -------------------------------
# Define prepareGeneActivityInput() function to prepare input for gene activity estimation
# -------------------------------
prepareGeneActivityInput <- function(bed, ann, meta, peaks, dirout = "GeneActivity") {
  bed    <- normalizePath(bed)
  meta   <- normalizePath(meta)
  ann    <- normalizePath(ann)
  peaks  <- normalizePath(peaks)
  dirout <- normalizePath(dirout, mustWork = FALSE)
  if (!dir.exists(dirout)) dir.create(dirout, recursive = TRUE)

  obj <- list()
  obj$bedpath <- bed
  
  # --------------------------------------------------------------------------
  message("Loading BED file (Tn5 insertions)...")
  obj$bed <- fread(bed, sep = "\t", header = FALSE, data.table = FALSE)
  if (ncol(obj$bed) < 4) stop("ERROR: BED file must have at least 4 columns")

  # ---------------------------------------------------------------------------
  message("Loading metadata file ...")
  obj$metadata <- if (is.character(meta)) {
    read.table(meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    stop("ERROR: meta must be a file path")
  }
  stopifnot("cellID" %in% colnames(obj$metadata))
  bed_cellIDs  <- as.character(obj$bed[, 4])
  meta_cellIDs <- as.character(obj$metadata$cellID)
  valid_cells  <- meta_cellIDs %in% bed_cellIDs
  if (!any(valid_cells)) stop("ERROR: No cellIDs in metadata found in BED file.")

  obj$metadata <- obj$metadata[valid_cells, , drop = FALSE]
  obj$bed <- obj$bed[bed_cellIDs %in% obj$metadata$cellID, ]
  message("Metadata cells: ", length(meta_cellIDs))
  message("BED cells: ", length(unique(bed_cellIDs)))
  message("Cells retained for analysis: ", nrow(obj$metadata)) 

  # -----------------------------------------------------------------------------
  message("Loading genome annotation file ...")
  anntype <- if (grepl("\\.gtf$", ann, ignore.case = TRUE)) "gtf" else "gff3"
  obj$gff <- makeTxDbFromGFF(ann, format = anntype, dbxrefTag = "Parent")

  # -----------------------------------------------------------------------------  
  message("Loading peak file ...")
  obj$acr <- fread(peaks, sep = "\t", header = FALSE, data.table = FALSE)
  if (ncol(obj$acr) < 3) stop("ERROR: Peak file must have at least 3 columns")

  return(obj)
}



# =========================================
# Prepare input and estimate gene activity
# =========================================
obj <- prepareGeneActivityInput(bed, ann, meta, peaks, dirout = outdir)

message("Estimating gene activity ...")
obj <- Socrates:::estGeneActivity(obj,
                                  pRange = 500,
                                  FeatureName = 'gene',
                                  con = 5000)

# =========================================
# Save outputs
# =========================================
out_prefix <- file.path(outdir, "merged")

saveRDS(obj$gene_activity, file = paste0(out_prefix, ".raw_gene_activity.rds"))
saveRDS(obj$peak_gene_score, file = paste0(out_prefix, ".peak_gene_score.rds"))
saveRDS(obj$gba, file = paste0(out_prefix, ".gba.rds"))
saveRDS(obj$pa, file = paste0(out_prefix, ".pa.rds"))
saveRDS(obj$pgw, file = paste0(out_prefix, ".pgw.rds"))

# Gene activity scores were scaled to sum to 10,000 per nucleus
cell.names <- colnames(obj$gene_activity)
obj$scaled_gene_activity <- obj$gene_activity %*% Diagonal(x = 10000 / Matrix::colSums(obj$gene_activity))
colnames(obj$scaled_gene_activity) <- cell.names
saveRDS(obj$scaled_gene_activity, file = paste0(out_prefix, ".scaled_gene_activity.rds"))


message("All outputs saved to: ", outdir)

