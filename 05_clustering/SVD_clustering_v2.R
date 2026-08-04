# run Socrates on merged socrates object #
# libraries
library(Socrates)
library(harmony)
library(igraph)
library(Matrix)
library(RcppML)
library(tidyverse)
library(mclust)
library(ggplot2)
library(cluster)

# functions --------------------------------------------------------------------------------------------------------------------------------
# This function uses a modified version of the original reduceDims function from the Socrates package.
# In addition to removing PCs highly correlated with read depth, it also considers correlation with pTSS, FRiP, and pOrg. 
# Any PC that is highly correlated with any one of these metrics will also be removed.
# ------------------------------------------------------------------------------------------------------------------------------------------
reduceDims <- function(obj,
                       method = "SVD",
                       n.pcs = 50,
                       scaleVar = TRUE,
                       num.var = 5000,
                       regNum = 5000,
                       cor.max.depth = 0.75,
                       cor.max.pTSS = 0.75,
                       cor.max.FRiP = 0.75,
                       cor.max.pOrg = 0.75,
                       doL2 = FALSE,
                       doL1 = FALSE,
                       doSTD = TRUE,
                       refit_residuals = FALSE,
                       residuals_slotName = "residuals",
                       svd_slotName = "PCA",
                       verbose = FALSE,
                       ...) {
  
  # sub functions
  l2norm <- function(x) { x / sqrt(sum(x^2)) }
  l1norm <- function(x) { x / sum(x) }
  RowVar <- function(x) {
    spm <- t(x)
    stopifnot(methods::is(spm, "dgCMatrix"))
    ans <- sapply(seq.int(spm@Dim[2]), function(j) {
      if (spm@p[j + 1] == spm@p[j]) return(0)
      mean <- sum(spm@x[(spm@p[j] + 1):spm@p[j + 1]]) / spm@Dim[1]
      sum((spm@x[(spm@p[j] + 1):spm@p[j + 1]] - mean)^2) +
        mean^2 * (spm@Dim[1] - (spm@p[j + 1] - spm@p[j]))
    }) / (spm@Dim[1] - 1)
    names(ans) <- spm@Dimnames[[2]]
    ans
  }

  # check residual slot
  if (is.null(obj[[residuals_slotName]])) {
    message("ERROR: slot ", residuals_slotName, " is missing.")
    stop("Exiting.")
  }

  # variable feature selection
  if (!is.null(num.var)) {
    if (verbose) message(" - identifying variable features...")
    row.var <- RowVar(obj[[residuals_slotName]])
    row.means <- Matrix::rowMeans(obj[[residuals_slotName]])
    adj.row.var <- loess(row.var ~ row.means)$residuals
    names(adj.row.var) <- names(row.var)
    adj.row.var <- adj.row.var[order(adj.row.var, decreasing = TRUE)]

    if (num.var >= 100) {
      topSites <- names(head(adj.row.var, n = num.var))
    } else {
      topSites <- names(adj.row.var[adj.row.var > num.var])
      if (length(topSites) < 50) {
        topSites <- names(adj.row.var[1:100])
      }
    }

    if (verbose) message(" - keeping ", length(topSites), " variable features...")

    if (refit_residuals & obj$norm_method == "tfidf") {
      if (regNum > length(topSites)) regNum <- length(topSites)
      test.dat <- list(counts = obj$counts[topSites, ], meta = obj$meta)
      M <- regModel(test.dat, subpeaks = regNum, verbose = verbose)$residuals
      M <- Matrix(t(apply(M, 1, function(x) x - min(x, na.rm = TRUE))), sparse = TRUE)
      M <- M[Matrix::rowSums(M) > 0, ]
    } else {
      M <- obj[[residuals_slotName]][topSites, ]
      if (obj$norm_method != "tfidf" & method == "NMF") {
        M <- t(apply(M, 1, function(x) x - min(x, na.rm = TRUE)))
      }
      M <- M[Matrix::rowSums(M) > 0, ]
    }
  } else {
    M <- obj[[residuals_slotName]]
  }

  # dimensionality reduction
  if (method == "SVD") {
    if (verbose) message(" - reducing dimensions with SVD...")
    obj$rdMethod <- "SVD"
    pcs <- irlba::irlba(t(M), nv = n.pcs)
    pc <- pcs$u
    if (scaleVar) {
      pc <- pc %*% diag(pcs$d)
    }
    pc[is.na(pc)] <- 0
  } else if (method == "NMF") {
    if (verbose) message(" - reducing dimensions with NMF...")
    obj$rdMethod <- "NMF"
    pcs <- RcppML::nmf(M, n.pcs, verbose = verbose, ...)
    pcs$u <- t(pcs$h)
    pcs$v <- pcs$w
    if (scaleVar) {
      pc <- t(pcs$h) %*% Diagonal(x = 1 / pcs$d)
    } else {
      pc <- t(pcs$h)
    }
    pc[is.na(pc)] <- 0
  }

  # add colnames
  rownames(pc) <- colnames(obj[[residuals_slotName]])
  colnames(pc) <- paste0("PC_", seq_len(ncol(pc)))

  # remove PCs with technical correlation
  if (verbose) message(" - removing PCs correlated with technical covariates...")

  meta <- obj$meta[rownames(pc), ]
  covariates <- list(
    depth = Matrix::colSums(obj$counts[, rownames(pc)]),
    pTSS = meta$pTSS,
    FRiP = meta$FRiP,
    pOrg = meta$pOrg
  )
  cor_mat <- sapply(covariates, function(var) {
    apply(pc, 2, function(u) cor(u, var, method = "spearman"))
  })
  cor_mat <- abs(cor_mat)
  thresh <- c(
    depth = cor.max.depth,
    pTSS = cor.max.pTSS,
    FRiP = cor.max.FRiP,
    pOrg = cor.max.pOrg
  )
  exceed_thresh <- sapply(names(thresh), function(k) cor_mat[, k] > thresh[k])
  remove_pcs <- apply(exceed_thresh, 1, any)

  if (verbose) {
    message("   - PCs removed: ", sum(remove_pcs))
    if (sum(remove_pcs) > 0) {
      print(round(cor_mat[remove_pcs, , drop = FALSE], 3))
    }
  }

  idx.keep <- !remove_pcs
  pc <- pc[, idx.keep, drop = FALSE]

  # normalization
  if (verbose && any(c(doL2, doL1, doSTD))) message(" - normalizing reduced dimensions...")
  if (doL2) {
    pc <- t(apply(pc, 1, l2norm))
  } else if (doL1) {
    pc <- t(apply(pc, 1, l1norm))
  } else if (doSTD) {
    pc <- t(apply(pc, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
  }
  pc[is.na(pc)] <- 0

  # return
  obj[[svd_slotName]] <- pc
  model_slot <- paste0(svd_slotName, "_model")
  obj[[model_slot]] <- list()
  obj[[model_slot]]$d <- pcs$d
  obj[[model_slot]]$v <- pcs$v
  obj[[model_slot]]$u <- pcs$u
  obj[[model_slot]]$keep_pcs <- idx.keep
  obj$hv_sites <- rownames(M)
  return(obj)
}
#--------------------------------------------------------------------------------------------------------


#########################################################################################################
#########################################################################################################
#########################################################################################################

# ================================= step0: get command-line arguments ===========================================
set.seed(1234)
args <- commandArgs(trailingOnly=T)
if(length(args) != 6){stop("Rscript SVD_analysis.R [rds] [feature_rate] [pcs] [knear] [res] [out_dir]")}

rds <- args[1]
feature_rate <- as.numeric(args[2])
pcs <- as.numeric(args[3])
knear <- as.numeric(args[4])
res <- as.numeric(args[5])
out_dir <- args[6]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# =================================== step1: load merged rds file ===============================================
soc.obj <- readRDS(rds)  #args[1]

# =================================== step2: clean sparse counts matrix =========================================
soc.obj <- cleanData(soc.obj,
                     min.c=100,       # minimum number of accessible features per cell
                     min.t=0.0025,    # remove accessible features that are too rare (in less than 0.25% of all nuclei)
                     max.t=0.005,     # remove accessible features that are too common as they provide little discriminatory power across cells (in more than 99.5% of all nuclei)
                     verbose=T)

# =================================== step3: normalize with TFIDF ===============================================
soc.obj <- tfidf(soc.obj, doL2=T)
number.sites <- ceiling(nrow(soc.obj$counts)*feature_rate)  # set highly variable features rate for further dimensionality reduction, feature_rate (#args[2]) ranges from 0 to 1

# =================================== step4: project with SVD using new reduceDims function ===================================================
soc.obj <- reduceDims(soc.obj,
                      method="SVD",
                      n.pcs=pcs,  #args[3]
                      cor.max.depth = 0.5, # remove PCs due to high correlation (>0.5) with nuclear read depth
                      cor.max.pTSS  = 0.5,  # remove PCs due to high correlation (>0.5) with pTSS
                      cor.max.FRiP  = 0.5,  # remove PCs due to high correlation (>0.5) with FRiP
                      cor.max.pOrg  = 0.5,  # remove PCs due to high correlation (>0.5) with pOrg
		      num.var=number.sites,
                      verbose=T,
                      scaleVar=T,
                      doSTD=F,
                      doL1=F,
                      doL2=T,
                      refit_residuals=F)

# =================================== step5: remove batch effects with harmony ==================================
ids <- rownames(soc.obj$PCA)
ref.obj <- RunHarmony(soc.obj$PCA,
		      meta_data=soc.obj$meta,
		      vars_use=c("pool_geno"),  # remove pool and genotype effects
                      theta=c(2),
                      sigma=0.1,
                      lambda=c(0.1),
                      nclust=50,
                      max.iter=30,
                      return_object=T)

# process
soc.obj$h.PCA <- t(ref.obj$Z_corr) 
colnames(soc.obj$h.PCA) <- paste0("h.PC_", 1:(ncol(soc.obj$PCA)))
rownames(soc.obj$h.PCA) <- ids

# ================================== step6: reduce to 2-dimensions with UMAP ===================================
soc.obj$l2.PCA <- t(apply(soc.obj$h.PCA, 1, function(x){(x/sqrt(sum(x^2)))}))  # l2 normalization: scale each cell vector to unit norm for cosine metric

soc.obj <- projectUMAP(soc.obj,
                       metric = "cosine", 
		       k.near = knear, #args[4]
                       svd_slotName = "l2.PCA", 
                       umap_slotName = "h.UMAP")

# ======== step7: identify clusters using Leiden graph-based clustering algorithm at best resolution ==============
soc.obj <- callClusters(soc.obj, 
                        res=res, #args[6]
                        k.near=knear, #args[4]
                        verbose=T,
                        cleanCluster=F,
                        cl.method=4,
                        e.thresh=3,
			threshold=3,
                        min.reads=1e6,
			m.clst=50,
                        svd_slotName="l2.PCA",
                        umap_slotName="h.UMAP",
                        cluster_slotName="h.Clusters")

# ============================ step8: plot cluster membership on UMAP embedding ====================================
pdf(file.path(out_dir, "merged.UMAP.harmony.clusters.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", cex=0.5)
dev.off()

pdf(file.path(out_dir, "merged.UMAP.harmony.library.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="library", cex=0.5, main = "library")
dev.off()

pdf(file.path(out_dir, "merged.UMAP.harmony.pool.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="pool", cex=0.5, main = "pool")
dev.off()

pdf(file.path(out_dir, "merged.UMAP.harmony.genotype.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="genotype", cex=0.5, main = "genotype")
dev.off()

pdf(file.path(out_dir, "merged.UMAP.harmony.pool_geno.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="pool_geno", cex=0.5, main = "pool_geno")
dev.off()


pdf(file.path(out_dir,"merged.UMAP.harmony.QC.pdf"), width=12, height=3)
layout(matrix(c(1:4), nrow=1, ncol = 4))
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="log10nSites", cex=0.2, main = "Total nSites(log10)")
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="pTSS", cex=0.2, main = "Fraction reads in TSS")
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="FRiP", cex=0.2, main = "Fraction reads in ACR")
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="pOrg", cex=0.2, main = "Fraction reads in organelle")
dev.off()


# ============================================== step9: save data ================================================
soc.obj$UMAP_model <- NULL   # remove model to allow save or it will trigger “Error: evaluation nested too deeply”

saveRDS(soc.obj, file = file.path(out_dir, "merged.processed.rds"))  # save full Socrates object

write.table(soc.obj$meta, file = file.path(out_dir, "merged.processed_meta.txt"), quote=F, row.names=T, col.names=T, sep="\t")  # write metadata

# # write reduced dimensions and clustering
write.table(soc.obj$l2.PCA, file = file.path(out_dir, "merged.reduced_dimensions.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(soc.obj$h.Clusters, file = file.path(out_dir,"merged.clustering.txt"), quote=F, row.names=T, col.names=T, sep="\t")

