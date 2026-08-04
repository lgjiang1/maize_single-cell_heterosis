####################################################################################
## Per-Family Motif Clustering — JASPAR2026 All-Plant
##
## Clusters 926 JASPAR2026 CORE plant motifs WITHIN each TF family using
## motifStack, preventing cross-family merges (e.g., bHLH E-box vs bZIP E-box).
## Adapted from 13_Arabidopsis_protoplast/5_TF_FP/final/01_motif_signatures/01a_motif_clustering.R
##
## Input (under 5_motif_scanning/):
##   - JASPAR2026_raw/: individual MEME files per motif (1526 total)
##   - signatures/motif_names.txt: 926 selected motifs (latest version per base_id)
##   - signatures/family_assignments.tsv: motif-to-family mapping (926 rows)
##
## Output (under 5_motif_scanning/signatures/):
##   - Plant_Motif_SignatureDB.meme: per-family signature PFMs (MEME format)
##   - Plant_MotifClusters.txt: cluster membership (MOTIFN\t<semicolon-joined members>)
##   - clustering_summary.tsv: per-family clustering summary
##   - Plant_Full_Motif_SignatureDB.rds: R object for later inspection
##
## Usage (Great Lakes):
##   module load GCC/11.2.0 OpenMPI/4.1.1 R/4.3.1
##   Rscript 02a_motif_clustering.R
####################################################################################

suppressMessages(library(tidyverse))
suppressMessages(library(motifStack))
suppressMessages(library(data.table))
library(igraph)
library(Rgraphviz)
library(ade4)
library(universalmotif)
library(TFBSTools)

set.seed(42)

## ---- Config ----

# Paths relative to project root.
# On cluster: /nfs/turbo/lsa-amarand/fabio_home/Projects/15_Heterosis
PROJ_ROOT  <- Sys.getenv("PROJ_ROOT",
                         unset = "/nfs/turbo/lsa-amarand/fabio_home/Projects/15_Heterosis")

JASPAR_DIR <- file.path(PROJ_ROOT, "5_motif_scanning", "JASPAR2026_raw")
SIG_DIR    <- file.path(PROJ_ROOT, "5_motif_scanning", "signatures")

## ---- Helper functions ----

ReadFirstMotif <- function(x) {
  tem <- importMatrix(x, format = "meme", to = "pfm")
  return(tem)
}

## ---- Read input data ----

# Read 926 selected motif names (format: MA0001.3_AGL3)
motif_names_raw <- read.table(file.path(SIG_DIR, "motif_names.txt"), h = FALSE)$V1
# Extract version ID for file matching (MA0001.3)
filesMotif <- sub("_.*", "", motif_names_raw)

# Filter to those available in JASPAR2026_raw
filescontrol <- gsub('.meme', '', list.files(JASPAR_DIR, pattern = "*meme"))
filesMotif <- filesMotif[filesMotif %in% filescontrol]
cat(sprintf("[INFO] %d motifs found in JASPAR2026_raw\n", length(filesMotif)))

# Read all motif PFMs
Plant_Motif <- lapply(paste0(JASPAR_DIR, '/', filesMotif, '.meme'), ReadFirstMotif)
Plant_Motif <- unlist(Plant_Motif)
cat(sprintf("[INFO] %d PFMs loaded\n", length(Plant_Motif)))

# Normalize names: replace dots with underscores for consistency
names(Plant_Motif) <- gsub('[.]', '_', names(Plant_Motif))

# Read family assignments
fam_df <- read.table(file.path(SIG_DIR, "family_assignments.tsv"), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")

# Build lookup: MAxxxx_version -> family
# Use motif_version (e.g. MA0001.3) rather than full motif_name (e.g. MA0001.3_AGL3)
# because TF names in MEME files may differ from the TSV
fam_df$motif_key <- gsub('[.]', '_', fam_df$motif_version)
family_lookup <- setNames(fam_df$tf_family, fam_df$motif_key)

# Map each loaded motif to its family using only the MAxxxx_version prefix
# e.g. names(Plant_Motif) = "MA1253_1_ERF036" -> extract "MA1253_1"
motif_names <- names(Plant_Motif)
motif_version_keys <- sub("^(MA[0-9]+_[0-9]+)_.*", "\\1", motif_names)
# Handle motifs whose full name IS the version (no TF suffix)
motif_version_keys[motif_version_keys == motif_names] <- motif_names[motif_version_keys == motif_names]
motif_families <- family_lookup[motif_version_keys]

# Check for unmapped
n_unmapped <- sum(is.na(motif_families))
if (n_unmapped > 0) {
  cat(sprintf("[WARN] %d motifs have no family assignment:\n", n_unmapped))
  cat(paste("  ", motif_names[is.na(motif_families)], collapse = "\n"), "\n")
  motif_families[is.na(motif_families)] <- "Unknown"
}

families <- unique(motif_families)
cat(sprintf("[INFO] %d unique families\n", length(families)))

## ---- Per-family clustering ----

all_signatures <- list()
cluster_table <- list()  # MOTIFN -> semicolon-joined members
summary_rows <- list()

sig_counter <- 0

for (fam in sort(families)) {
  fam_idx <- which(motif_families == fam)
  fam_motifs <- Plant_Motif[fam_idx]
  n_motifs <- length(fam_motifs)

  if (n_motifs == 0) next

  if (n_motifs == 1) {
    # Singleton: pass through as its own signature
    sig_counter <- sig_counter + 1
    sig_name <- paste0("MOTIF", sig_counter)
    member_name <- names(fam_motifs)[1]

    pfm_obj <- new("pfm", mat = fam_motifs[[1]]@mat, name = member_name)

    all_signatures[[sig_name]] <- pfm_obj
    cluster_table[[sig_name]] <- member_name

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      family = fam, n_input = 1, n_signatures = 1,
      stringsAsFactors = FALSE
    )
    cat(sprintf("[INFO] %s: 1 motif -> 1 signature (singleton)\n", fam))
    next
  }

  if (n_motifs == 2) {
    # Two motifs: cluster but motifSignature needs >= 3 leaves for phylog
    pfms_fam <- mapply(fam_motifs, names(fam_motifs),
                       FUN = function(.pfm, .name) {
                         new("pfm", mat = fam_motifs[[.name]]@mat, name = .name)
                       })

    tryCatch({
      hc <- clusterMotifs(pfms_fam)
      if (max(hc$height) < 0.5) {
        sig_counter <- sig_counter + 1
        sig_name <- paste0("MOTIF", sig_counter)
        m1 <- pfms_fam[[1]]@mat
        m2 <- pfms_fam[[2]]@mat
        maxw <- max(ncol(m1), ncol(m2))
        if (ncol(m1) < maxw) m1 <- cbind(m1, matrix(0.25, 4, maxw - ncol(m1)))
        if (ncol(m2) < maxw) m2 <- cbind(m2, matrix(0.25, 4, maxw - ncol(m2)))
        avg_mat <- (m1 + m2) / 2
        merged_name <- paste0(names(fam_motifs), collapse = ";")
        pfm_merged <- new("pfm", mat = avg_mat, name = merged_name)
        all_signatures[[sig_name]] <- pfm_merged
        cluster_table[[sig_name]] <- merged_name
        cat(sprintf("[INFO] %s: 2 motifs -> 1 signature (merged, height=%.3f)\n",
                    fam, max(hc$height)))
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          family = fam, n_input = 2, n_signatures = 1,
          stringsAsFactors = FALSE
        )
      } else {
        for (j in 1:2) {
          sig_counter <- sig_counter + 1
          sig_name <- paste0("MOTIF", sig_counter)
          all_signatures[[sig_name]] <- pfms_fam[[j]]
          cluster_table[[sig_name]] <- names(fam_motifs)[j]
        }
        cat(sprintf("[INFO] %s: 2 motifs -> 2 signatures (distinct, height=%.3f)\n",
                    fam, max(hc$height)))
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          family = fam, n_input = 2, n_signatures = 2,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      for (j in 1:2) {
        sig_counter <<- sig_counter + 1
        sig_name <- paste0("MOTIF", sig_counter)
        all_signatures[[sig_name]] <<- pfms_fam[[j]]
        cluster_table[[sig_name]] <<- names(fam_motifs)[j]
      }
      cat(sprintf("[WARN] %s: 2 motifs clustering failed (%s), keeping both\n",
                  fam, conditionMessage(e)))
      summary_rows[[length(summary_rows) + 1]] <<- data.frame(
        family = fam, n_input = 2, n_signatures = 2,
        stringsAsFactors = FALSE
      )
    })
    next
  }

  # n_motifs >= 3: full motifStack clustering
  pfms_fam <- mapply(fam_motifs, names(fam_motifs),
                     FUN = function(.pfm, .name) {
                       new("pfm", mat = fam_motifs[[.name]]@mat, name = .name)
                     })

  tryCatch({
    hc <- clusterMotifs(pfms_fam)
    phylog <- hclust2phylog(hc)
    leaves <- names(phylog$leaves)
    pfms_ordered <- pfms_fam[leaves]

    motifSig <- motifSignature(pfms_ordered, phylog,
                               cutoffPval = 0.0001, min.freq = 1)
    sig <- signatures(motifSig)
    n_sig <- length(sig)

    for (j in 1:n_sig) {
      sig_counter <- sig_counter + 1
      sig_name <- paste0("MOTIF", sig_counter)

      all_signatures[[sig_name]] <- sig[[j]]
      cluster_table[[sig_name]] <- sig[[j]]@name
    }

    cat(sprintf("[INFO] %s: %d motifs -> %d signatures\n", fam, n_motifs, n_sig))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      family = fam, n_input = n_motifs, n_signatures = n_sig,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    for (j in seq_along(fam_motifs)) {
      sig_counter <<- sig_counter + 1
      sig_name <- paste0("MOTIF", sig_counter)
      pfm_obj <- new("pfm", mat = fam_motifs[[j]]@mat,
                      name = names(fam_motifs)[j])
      all_signatures[[sig_name]] <<- pfm_obj
      cluster_table[[sig_name]] <<- names(fam_motifs)[j]
    }
    cat(sprintf("[WARN] %s: clustering failed (%s), keeping all %d as singletons\n",
                fam, conditionMessage(e), n_motifs))
    summary_rows[[length(summary_rows) + 1]] <<- data.frame(
      family = fam, n_input = n_motifs, n_signatures = n_motifs,
      stringsAsFactors = FALSE
    )
  })
}

cat(sprintf("\n[INFO] Total: %d signatures from %d motifs across %d families\n",
            length(all_signatures), length(Plant_Motif), length(families)))

## ---- Write output ----

# 1. Write MEME file
um_list <- lapply(names(all_signatures), function(sname) {
  pfm_obj <- all_signatures[[sname]]
  mat <- pfm_obj@mat
  create_motif(mat, type = "PPM", name = pfm_obj@name,
               alphabet = "DNA")
})

write_meme(um_list, file.path(SIG_DIR, "Plant_Motif_SignatureDB.meme"), overwrite = TRUE)
cat(sprintf("[INFO] Wrote Plant_Motif_SignatureDB.meme (%d signatures)\n",
            length(um_list)))

# 2. Write cluster table
cluster_lines <- paste0(names(cluster_table), "\t",
                        unlist(cluster_table))
writeLines(cluster_lines, file.path(SIG_DIR, "Plant_MotifClusters.txt"))
cat(sprintf("[INFO] Wrote Plant_MotifClusters.txt (%d lines)\n",
            length(cluster_lines)))

# 3. Write summary
summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file.path(SIG_DIR, "clustering_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("[INFO] Wrote clustering_summary.tsv\n"))

# 4. Save R objects for later inspection
Object_motifStack <- list(
  all_signatures = all_signatures,
  cluster_table = cluster_table,
  family_lookup = family_lookup,
  summary = summary_df
)
saveRDS(Object_motifStack, file.path(SIG_DIR, "Plant_Full_Motif_SignatureDB.rds"))
cat(sprintf("[INFO] Wrote Plant_Full_Motif_SignatureDB.rds\n"))

cat("\n[DONE] Per-family motif clustering complete.\n")
