library(data.table)

# Usage: Rscript ecdf_filter.R <prefix> <base> <min_prop> <fdr> <min_cell>
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript ecdf_filter.R <prefix> <base> <fdr_cut> <min_prop> <min_cell>")
}
prefix    <- args[1]
base      <- args[2]
fdr_cut   <- as.numeric(args[3])   # e.g., 0.01
min_prop  <- as.numeric(args[4])   # e.g., 0.01
min_cell  <- as.integer(args[5])   # e.g., 10
workdir   <- file.path(base, "work", prefix)


# --------------------------- Background ----------------------------------
# Each control region: fraction of cells with >=1 Tn5 insertion (0..1).
bg <- fread(file.path(workdir, sprintf("%s.bg.p.bed", prefix)), header=FALSE)
setnames(bg, c("chr","start","end","p"))
bg$p <- suppressWarnings(as.numeric(bg$p))
bg   <- bg[is.finite(p)]

# Count unique cells from tn5.bed (robust to empty/missing files).
Ncells <- suppressWarnings(as.numeric(system(
  sprintf("cut -f4 %s/%s.tn5.bed | sort -u | wc -l", file.path(base), prefix),
  intern=TRUE)))
if (!is.finite(Ncells)) Ncells <- 1

# Clamp boundary values into (0,1) for numerical stability.
eps <- max(1e-6, 1/(2*max(1, Ncells)))
bgp <- bg$p
bgp[bgp <= 0] <- eps
bgp[bgp >= 1] <- 1 - eps


# -------------------------------- Peaks -----------------------------------
pk   <- fread(file.path(workdir, sprintf("%s.peaks.p.tsv", prefix)), header=FALSE)
pcol <- ncol(pk)
p_peak_raw <- suppressWarnings(as.numeric(pk[[pcol]]))

p_peak <- p_peak_raw
p_peak[!is.finite(p_peak)] <- NA_real_
p_peak[p_peak <= 0] <- eps
p_peak[p_peak >= 1] <- 1 - eps


# -------------------- ECDF right-tail p-values --------------------
# Null is the empirical CDF of background fractions.
Fhat <- ecdf(bgp)
pval <- 1 - Fhat(p_peak)  # right-tail probability: P(BG >= p_peak)

pval[!is.finite(pval)] <- 1
qval <- p.adjust(pval, method = "fdr")


# ---- dual filter: FDR<fdr_cut & max(min_prop * Ncells, min_cell) cells accessible ----
# Require at least max( ceil(min_prop * Ncells), min_cell ) cells per peak.
min_cells_prop     <- ceiling(min_prop * Ncells)
min_cells_required <- max(min_cells_prop, min_cell)
n_cells_peak <- as.integer(round(p_peak_raw * Ncells))
orig_cols <- 5:(pcol - 2)
keep <- (qval < fdr_cut) & (n_cells_peak >= min_cells_required)


# -------------------- Outputs --------------------
pk_out <- cbind(
  pk[, ..orig_cols],
  p_peak_raw   = p_peak_raw,
  n_cells_peak = n_cells_peak,
  p_value      = pval,
  FDR          = qval
)

fwrite(pk_out,
       file.path(workdir, sprintf("%s.with_p_FDR.tsv", prefix)),
       sep = "\t", col.names = TRUE)

orig_n <- length(orig_cols)

fwrite(pk_out[keep, 1:orig_n, with = FALSE],
       file.path(workdir, sprintf("%s.filtered.bed", prefix)),
       sep = "\t", col.names = FALSE)

# -------------------- Logging & parameters --------------------
total_peaks  <- nrow(pk_out)
passed       <- sum(keep)
filtered_out <- total_peaks - passed
bg_unique    <- uniqueN(bgp)
bg_mean      <- mean(bgp)
pk_mean      <- mean(p_peak_raw, na.rm = TRUE)

message(sprintf(
  "[INFO] %s: method=ecdf; total=%d; passed(FDR<%.3g & cells>=max(%d, ceil(%.3g*%d)=%d))=%d; filtered=%d; Ncells=%d; bg_unique=%d",
  prefix, total_peaks, fdr_cut,
  min_cell, min_prop, Ncells, min_cells_required,
  passed, filtered_out, Ncells, bg_unique))

# Keep a params file; alpha/beta=NA for downstream compatibility.
fwrite(data.table(
         alpha = NA_real_, beta = NA_real_, method = "ecdf",
         Ncells = Ncells,
         total_peaks = total_peaks, passed = passed, filtered_out = filtered_out,
         bg_mean = bg_mean, pk_mean = pk_mean,
         min_prop = min_prop, fdr_cut = fdr_cut,
         min_cell = min_cell,
         min_cells_prop = min_cells_prop,
         min_cells_required = min_cells_required,
         bg_unique = bg_unique),
       file.path(workdir, sprintf("%s.beta_params.tsv", prefix)),
       sep = "\t")
