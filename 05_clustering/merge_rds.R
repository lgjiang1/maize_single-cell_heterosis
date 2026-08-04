## Merge all *filtered.soc.rds files

library(Socrates)
library(Matrix)

# ============== Parse arguments ==================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript merge_rds.R <rds_folder> <out_dir>")
}
rds_dir <- args[1]
out_dir <- args[2]

# ======= Define mergeSocObjects() function =======
mergeSocObjects <- function(obj.list){
  # functions
  .merge.sparse <- function(cnt.list) {
    cnnew <- character()
    rnnew <- character()
    x <- vector()
    i <- numeric()
    j <- numeric()
    
    for (M in cnt.list) {
      cnold <- colnames(M)
      rnold <- rownames(M)
      cnnew <- union(cnnew,cnold)
      rnnew <- union(rnnew,rnold)
      cindnew <- match(cnold,cnnew)
      rindnew <- match(rnold,rnnew)
      ind <- summary(M)
      i <- c(i,rindnew[ind[,1]])
      j <- c(j,cindnew[ind[,2]])
      x <- c(x,ind[,3])
    }
    sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
  }
  # separate counts and meta
  counts <- lapply(obj.list, function(x){
    x$counts
  })
  counts <- .merge.sparse(counts)
  # meta
  metas <- lapply(obj.list, function(x){
    x$meta
  })
  metas <- do.call(rbind, metas)
  rownames(metas) <- metas$cellID
  metas <- metas[colnames(counts),]
  # new object
  new.obj <- list(counts=counts, meta=metas)
  return(new.obj)
}

# ================== Load and merge all filtered.soc.rds files =======================
rds_files <- list.files(rds_dir, pattern = "\\.filtered\\.soc\\.rds$", full.names = TRUE)
obj_list <- lapply(rds_files, readRDS)
cat("Merging", length(obj_list), "RDS files from", rds_dir, "\n")

merged <- mergeSocObjects(obj_list)
merged$counts <- merged$counts[Matrix::rowSums(merged$counts) > 0, ]

stopifnot(identical(colnames(merged$counts), rownames(merged$meta)))  ## consistency check before saving

# ======== Add back quality control metrics (pTSS, FRiP, pOrg) ==================
merged$meta$pTSS <- merged$meta$tss / merged$meta$total
merged$meta$FRiP <- merged$meta$acrs / merged$meta$total
merged$meta$pOrg <- merged$meta$ptmt / merged$meta$total

# =================== Save output ====================
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(merged, file.path(out_dir, "merged.soc.rds"))

meta_out <- file.path(out_dir, "merged_soc.meta.txt")
write.table(merged$meta, file = meta_out, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

