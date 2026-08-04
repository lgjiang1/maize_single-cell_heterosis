## ===== pseudo-bulk DEG and fgsea analysis =====
## Rscript pseudobulk_DAR-fgsea.R --counts $raw_geneact --meta $meta --gmt $gmt --minSize 10 --maxSize 600 --nperm 10000 --markers $markers --out $outdir

library(optparse)
library(limma)
library(edgeR)
library(fgsea)
library(Matrix)
library(data.table)

## -------- Command line arguments --------
option_list <- list(
  make_option(c("--counts"),  type="character", help="RDS: counts (features x cells)"),
  make_option(c("--meta"),    type="character", help="TSV: metadata with cellID, LouvainClusters, library"),
  make_option(c("--gmt"),     type="character", help="GMT file"),
  make_option(c("--minSize"), type="integer",  default=10,    help="min geneset size [default %default]"),
  make_option(c("--maxSize"), type="integer",  default=600,   help="max geneset size [default %default]"),
  make_option(c("--nperm"),   type="integer",  default=10000, help="fgsea permutations [default %default]"),
  make_option(c("--markers"), type="character", default=NULL, help="TSV with header; column name 'geneID'"),
  make_option(c("--out"),     type="character", default="gsea_out", help="output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
counts_rds  <- opt$counts
meta_tsv    <- opt$meta
gmt_file    <- opt$gmt
min_gsize   <- opt$minSize
max_gsize   <- opt$maxSize
nperm       <- opt$nperm
markers_file<- opt$markers
outdir      <- opt$out
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
set.seed(123)


## -------- Step1: Load data and align cells --------
cat("1) Loading data and aligning cells...\n")
X <- readRDS(counts_rds)                         # genes x cells gene activity raw counts
meta <- read.table(meta_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
meta <- as.data.frame(meta)
rownames(meta) <- as.character(meta[[1]])

req_cols <- c("LouvainClusters","library")       
if (!all(req_cols %in% colnames(meta))) {
    stop(paste0("Metadata must contain columns: ", paste(req_cols, collapse=", ")))
}
lib_col <- "library"                              

cells <- intersect(colnames(X), rownames(meta))
X <- X[, cells, drop=FALSE]
meta <- meta[cells, , drop=FALSE]
clu <- factor(meta$LouvainClusters); names(clu) <- rownames(meta)
libs <- factor(meta[[lib_col]]);    names(libs) <- rownames(meta)

# markers (header with geneID, name and type)
marker_set <- NULL
name_map   <- NULL
if (!is.null(markers_file)) {
    mk <- read.table(markers_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    marker_set <- unique(as.character(mk$geneID))
    name_map   <- setNames(as.character(mk$name), as.character(mk$geneID))
    mk$type <- as.character(mk$type)
    mk_type_unique <- unique(mk[, c("geneID","type")])
    type_agg <- aggregate(type ~ geneID, mk_type_unique,
                        FUN=function(v) paste(unique(v[!is.na(v) & v!=""]), collapse=";"))
    type_map <- setNames(as.character(type_agg$type), as.character(type_agg$geneID))
    cat("Loaded", length(marker_set), "markers from", markers_file, "\n")
} 

## -------- Step2: Read GMT and intersect with matrix genes --------
cat("2) Reading GMT and filtering gene sets by size...\n")
pathways_raw <- gmtPathways(gmt_file)
genes_mat <- rownames(X)

# Keep only genes present in the counts matrix
pathways <- lapply(pathways_raw, function(v) intersect(v, genes_mat))

# Filter by size thresholds (default: 10 <= size <= 600)
lens <- vapply(pathways, length, integer(1))
keep <- lens >= min_gsize & lens <= max_gsize
pathways <- pathways[keep]
cat("Retained gene sets:", length(pathways), "\n\n")
if (length(pathways) == 0) stop("No usable gene sets after filtering.")
genes_in_pathways <- unique(unlist(pathways))

# Build GO-ID -> description map from the GMT (col1 = GO ID, col2 = description)
gmt_df <- read.table(gmt_file, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                     fill=TRUE, quote="", comment.char="", check.names=FALSE)
colnames(gmt_df)[1:2] <- c("pathway", "go_description")
desc_map <- setNames(gmt_df$go_description, gmt_df$pathway)
desc_map <- desc_map[names(pathways)] # only keep retained pathway description


## -------- Step3: Build unified reference panel (TOTAL = mean cluster size) --------
cat("3) Building unified reference panel with TOTAL = mean cluster size...\n")
sz <- as.numeric(table(clu))        # cells per cluster
C  <- length(sz)                    # number of clusters
total_ref_target <- floor(mean(sz)) # mean cluster size

if (total_ref_target < 2)
    stop("Mean cluster size < 2; not enough cells to build a reference.")
if ((total_ref_target %/% C) == 0)
    stop("Base per-cluster quota is 0 under 'drop remainders'. Increase TOTAL or relax the rule.")

per_clu_targets  <- rep(total_ref_target %/% C, C)
names(per_clu_targets) <- levels(clu)

# Sample per cluster into unified ref pool
ref_cells <- unlist(lapply(levels(clu), function(k){
  idx  <- names(clu)[clu == k]
  take <- per_clu_targets[[k]]
  if (take <= 0) return(character(0))
  if (length(idx) >= take) sample(idx, take, replace=FALSE) else sample(idx, take, replace=TRUE)
}), use.names=FALSE)

cat("Reference TOTAL target (approx):", sum(per_clu_targets),
    " | per-cluster quota:", paste(per_clu_targets, collapse=","), "\n\n")

## -------- Helper --------
sumcols <- function(m) if (inherits(m,"Matrix")) Matrix::rowSums(m) else base::rowSums(m)
get_lib_from_colname <- function(nm) sub("^.*_", "", nm)

## -------- Step4: Per-cluster analysis (pseudo-bulk DEG analysis + fgsea) --------
cat("4) Per-cluster DEG analysis + fgsea...\n")
all_clusters <- levels(clu)

for (k in all_clusters){
    cat("- Processing cluster:", k, "\n")
    k_cells <- names(clu)[clu == k]
    if (length(k_cells) < 2) { warning("Cluster ", k, " has <2 cells; skip."); next }

  ## 4.1 Cluster-by-library pseudobulk
    libs_k <- split(k_cells, libs[k_cells])                 
    Yk_list <- lapply(libs_k, function(cs) sumcols(X[, cs, drop=FALSE]))
    if (length(Yk_list) == 0) { warning(" No library found for cluster ", k, "; skip."); next }
    Yk <- do.call(cbind, Yk_list)
    colnames(Yk) <- paste0("cluster_", k, "_", names(libs_k))

  ## 4.2 Unified-ref-by-library pseudobulk  
    libs_ref <- split(ref_cells, libs[ref_cells])          
    Yr_list <- lapply(libs_ref, function(cs) sumcols(X[, cs, drop=FALSE]))
    if (length(Yr_list) == 0) { warning("Unified ref has no libraries; skip."); next }
    Yr <- do.call(cbind, Yr_list)
    colnames(Yr) <- paste0("ref_", k, "_", names(libs_ref))
  
  # 4.3 Combine two datasets and build factors
    Y <- cbind(Yk, Yr)
    group <- factor(c(rep("cluster", ncol(Yk)), rep("ref", ncol(Yr))))
    group <- relevel(group, ref = "ref")  # baseline = ref
  
    library_factor <- factor(c(get_lib_from_colname(colnames(Yk)),
                             get_lib_from_colname(colnames(Yr))))
    tab_gl <- table(group, library_factor)
    cat("Samples per group × library:\n"); print(tab_gl)

  ## 4.4 edgeR pipeline
    dge <- DGEList(counts=Y, group=group)
    keep <- filterByExpr(dge, group=group)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge, method="TMM")

    design <- model.matrix(~ group + library_factor)
  
  # Check if design matrix is full rank; warns if group and library are confounded.
    if (qr(design)$rank < ncol(design)) {
        warning("Design is not full rank (group ~ library confounding). ",
            "Consider dropping/merging levels OR using only common libraries for cluster ", k, ".")
    }

    dge <- estimateDisp(dge, design, robust=TRUE)
    fit <- glmQLFit(dge, design, robust=TRUE)
    qlf  <- glmQLFTest(fit, coef = "groupcluster")

   ## 4.5 Save edgeR full table and marker-only table
    tt  <- topTags(qlf, n=Inf, sort.by="none")$table
    out <- tt[order(tt$PValue), ]
    out$geneID <- rownames(out)
    out <- out[, c("geneID", setdiff(colnames(out), "geneID"))]
    outfile <- file.path(outdir, sprintf("edgeR_allGenes_%s.tsv", make.names(k)))
    fwrite(as.data.frame(out), outfile, sep="\t", quote=FALSE)
    cat("edgeR table saved:", outfile, "\n")

  # export marker-only table if markers were provided
    if (!is.null(marker_set)) {
        keep_idx <- which(out$geneID %in% marker_set)
        if (length(keep_idx) > 0) {
            out_marker <- out[keep_idx, , drop=FALSE]
	    out_marker <- out_marker[order(out_marker$logFC,decreasing = TRUE), ]
            out_marker$name <- name_map[out_marker$geneID]
            out_marker$type <- type_map[out_marker$geneID]
            out_marker <- out_marker[, c("geneID","name","type",setdiff(colnames(out_marker), c("geneID","name","type")))]
            outfile_mk <- file.path(outdir, sprintf("edgeR_markers_%s.tsv", make.names(k)))
            fwrite(out_marker, outfile_mk, sep="\t", quote=FALSE)
            cat("markers table:", outfile_mk, "\n")
        }else {
         warning("No markers matched for cluster ", k)
        }
    }
  

  ## 4.6 fgsea (rank by logFC, restricted to genes present in pathways)
    ranks <- out$logFC; names(ranks) <- rownames(out)
    ranks <- ranks[ intersect(names(ranks), genes_in_pathways) ]
    ranks <- ranks[is.finite(ranks)]
    if (length(ranks) < 10) {
        warning("Too few ranked genes in pathways for fgsea; skip fgsea for cluster ", k); 
        next
    }
    ranks <- sort(ranks, decreasing=TRUE)

    fg <- fgsea(pathways = pathways,
                stats    = ranks,
                nperm    = nperm,
                minSize  = min_gsize,
                maxSize  = max_gsize)

    if ("leadingEdge" %in% names(fg) && is.list(fg$leadingEdge)) {
        fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(x, collapse=","), "")
    }

  # Order by FDR then NES
    fg <- fg[order(fg$padj, -fg$NES), ]
  
  # Attach GO description and place it right after 'pathway'
    fg$go_description <- unname(desc_map[fg$pathway])
    fg <- as.data.frame(fg)
    fg <- fg[, c("pathway","go_description", setdiff(names(fg), c("pathway","go_description")))]
    fg_sig <- fg[fg$padj < 0.05, , drop=FALSE]
    fg_sig_pos <- fg[fg$padj < 0.05 & fg$NES > 0, , drop = FALSE]

  # Save fgsea outputs
    out_full <- file.path(outdir, sprintf("fgsea_%s.tsv", make.names(k)))
    out_sig  <- file.path(outdir, sprintf("fgsea_FDRlt0.05_%s.tsv", make.names(k)))
    out_sig_pos  <- file.path(outdir, sprintf("fgsea_FDRlt0.05_NESgt0_%s.tsv", make.names(k)))
    fwrite(fg, out_full, sep="\t", quote=FALSE)
    fwrite(fg_sig, out_sig,  sep="\t", quote=FALSE)
    fwrite(fg_sig_pos,  out_sig_pos,  sep="\t", quote=FALSE)
    cat("fgsea saved:", out_full, "and", out_sig, "and", out_sig_pos, "\n")
}

cat("DEG and fgsea analysis Done!!\n")
