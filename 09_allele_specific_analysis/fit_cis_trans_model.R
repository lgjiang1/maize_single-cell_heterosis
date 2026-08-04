# ==================================================================================================
# Joint cis/trans model for one cell type
#
# Main classification:
#   non_divergent / cis / trans / cis_trans
#   using joint parent-NB + hybrid-BB model selection by BIC
#
# Modified cis_trans subclassification:
#   x = F1 average log2 allelic ratio
#   y = F0 average log2 fold change
#
#   F0 summary:
#     - use CPM for simple library-size normalization
#     - calculate replicate-wise log2FC
#     - average across replicates
#
#   F1 summary:
#     - calculate replicate-wise log2(hap1/hap2)
#     - average across replicates
#
#   Subtype rules:
#     1) CIS-trans (opposite direction with cis stronger than trans): x*y > 0 and |x| > |y|
#     2) TRANS-cis (opposite direction with trans stronger than cis): x*y < 0
#     3) CIS+trans (same direction with cis stronger than trans): x*y > 0 and |x| < |y| < 2|x|
#     4) TRANS+cis (same direction with trans stronger than cis): x*y > 0 and |y| > 2|x|
# ==================================================================================================
library(data.table)
library(edgeR)

# --------------------------------------------------------------------------------------------------
# Parse arguments
# --------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript step02a_fit_cis_trans_model.R <cell_type>\n",
       "Example: Rscript step02a_fit_cis_trans_model.R 5")
}
ct <- as.integer(args[1])

# --------------------------------------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------------------------------------
hap1_name <- "B73"
hap2_name <- "Ki3"
info_file <- "0_input/All_B-K_celltypes_with_inheritance_plus_extra_ACR_info.tsv"  # This file is used at final step for adding extra information to final cis-trans results

in_dir    <- "1_merged"
out_dir   <- "2_cis_trans"
final_dir <- file.path(out_dir, "final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
in_file <- file.path(in_dir, paste0("merged_", hap1_name, "_", hap2_name, "_counts_without_exonic_ACR.tsv"))

# Filtering thresholds
min_parent_total <- 10            # minimum total reads across both parents per peak
min_hybrid_total <- 20            # minimum total reads across all hybrid samples per peak
min_hybrid_rep_count <- 10        # minimum reads per hybrid replicate to be considered valid
min_hybrid_supported_reps <- 2    # minimum number of valid hybrid replicates required

# Optimization
optim_method <- "BFGS"   # optimization algorithm used for parameter fitting
optim_maxit  <- 1000     # maximum number of optimization iterations

# Numerical protection
eps <- 1e-8

# Summary-effect settings for cis_trans subclassification
summary_pseudocount <- 0.5

# --------------------------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------------------------
clip01 <- function(x, eps = 1e-8) {pmin(pmax(x, eps), 1 - eps)}
safe_exp <- function(x) {exp(pmin(pmax(x, -50), 50))}

odds_from_p <- function(p) {
  p <- clip01(p, eps)
  p / (1 - p)
}

# Negative binomial log-likelihood (NB)
# var = mu + phi * mu^2
nb_loglik_vec <- function(counts, mu, phi) {
  mu <- pmax(mu, eps)
  phi <- pmax(phi, eps)
  size <- 1 / phi
  sum(dnbinom(x = counts, mu = mu, size = size, log = TRUE))
}

# Beta-binomial log-likelihood using alpha/beta parameterization (BB)
bb_loglik_vec_ab <- function(z, n, alpha, beta) {
  alpha <- pmax(alpha, eps)
  beta <- pmax(beta, eps)

  sum(
    lchoose(n, z) +
      lbeta(z + alpha, n - z + beta) -
      lbeta(alpha, beta)
  )
}

calc_bic <- function(loglik, k, n_obs) {
  -2 * loglik + k * log(n_obs)}

# --------------------------------------------------------------------------------------------------
# Model log-likelihoods
# --------------------------------------------------------------------------------------------------
ll_non_divergent <- function(par, p1, p2, sf1, sf2, z, n, phi) {
  eta_p <- par[1]
  log_m <- par[2]
  log_k <- par[3]

  p_shared <- plogis(eta_p)
  m <- safe_exp(log_m)
  k <- safe_exp(log_k)

  od <- odds_from_p(p_shared)

  mu1 <- sf1 * m * od
  mu2 <- sf2 * m * od

  ll_parent <- nb_loglik_vec(p1, mu1, phi) +
    nb_loglik_vec(p2, mu2, phi)

  ll_hybrid <- bb_loglik_vec_ab(z, n, alpha = k, beta = k)

  ll_parent + ll_hybrid
}

ll_cis <- function(par, p1, p2, sf1, sf2, z, n, phi) {
  eta_p_mu <- par[1]
  eta_p_nu <- par[2]
  log_m <- par[3]
  log_c <- par[4]

  p_mu <- plogis(eta_p_mu)
  p_nu <- plogis(eta_p_nu)
  m <- safe_exp(log_m)
  c0 <- safe_exp(log_c)

  od_mu <- odds_from_p(p_mu)
  od_nu <- odds_from_p(p_nu)

  mu1 <- sf1 * m * od_mu
  mu2 <- sf2 * m * od_nu

  alpha <- c0 * od_mu
  beta <- c0 * od_nu

  ll_parent <- nb_loglik_vec(p1, mu1, phi) +
    nb_loglik_vec(p2, mu2, phi)

  ll_hybrid <- bb_loglik_vec_ab(z, n, alpha = alpha, beta = beta)

  ll_parent + ll_hybrid
}

ll_trans <- function(par, p1, p2, sf1, sf2, z, n, phi) {
  eta_p_mu <- par[1]
  eta_p_nu <- par[2]
  log_m <- par[3]
  log_k <- par[4]

  p_mu <- plogis(eta_p_mu)
  p_nu <- plogis(eta_p_nu)
  m <- safe_exp(log_m)
  k <- safe_exp(log_k)

  od_mu <- odds_from_p(p_mu)
  od_nu <- odds_from_p(p_nu)

  mu1 <- sf1 * m * od_mu
  mu2 <- sf2 * m * od_nu

  ll_parent <- nb_loglik_vec(p1, mu1, phi) +
    nb_loglik_vec(p2, mu2, phi)

  ll_hybrid <- bb_loglik_vec_ab(z, n, alpha = k, beta = k)

  ll_parent + ll_hybrid
}

ll_cis_trans <- function(par, p1, p2, sf1, sf2, z, n, phi) {
  eta_p_mu <- par[1]
  eta_p_nu <- par[2]
  log_m <- par[3]
  log_alpha <- par[4]
  log_beta <- par[5]

  p_mu <- plogis(eta_p_mu)
  p_nu <- plogis(eta_p_nu)
  m <- safe_exp(log_m)
  alpha <- safe_exp(log_alpha)
  beta <- safe_exp(log_beta)

  od_mu <- odds_from_p(p_mu)
  od_nu <- odds_from_p(p_nu)

  mu1 <- sf1 * m * od_mu
  mu2 <- sf2 * m * od_nu

  ll_parent <- nb_loglik_vec(p1, mu1, phi) +
    nb_loglik_vec(p2, mu2, phi)

  ll_hybrid <- bb_loglik_vec_ab(z, n, alpha = alpha, beta = beta)

  ll_parent + ll_hybrid
}

# --------------------------------------------------------------------------------------------------
# Fit one model
# --------------------------------------------------------------------------------------------------
fit_one_model <- function(model_name, p1, p2, sf1, sf2, z, n, phi) {
  mean1_norm <- max(mean(p1 / sf1), 0.1)
  mean2_norm <- max(mean(p2 / sf2), 0.1)

  p_mu_init <- clip01(mean1_norm / (mean1_norm + mean2_norm), eps)
  p_nu_init <- clip01(mean2_norm / (mean1_norm + mean2_norm), eps)

  m_init <- max(mean(c(mean1_norm, mean2_norm)), 0.1)

  hyb_frac_init <- if (sum(n) > 0) clip01(sum(z) / sum(n), eps) else 0.5
  alpha_init <- max(hyb_frac_init * 20, 1)
  beta_init <- max((1 - hyb_frac_init) * 20, 1)

  k_init <- 10
  c_init <- 10

  if (model_name == "non_divergent") {
    init <- c(qlogis(0.5), log(m_init), log(k_init))
    fn <- function(par) -ll_non_divergent(par, p1, p2, sf1, sf2, z, n, phi)
    k <- 3
  } else if (model_name == "cis") {
    init <- c(qlogis(p_mu_init), qlogis(p_nu_init), log(m_init), log(c_init))
    fn <- function(par) -ll_cis(par, p1, p2, sf1, sf2, z, n, phi)
    k <- 4
  } else if (model_name == "trans") {
    init <- c(qlogis(p_mu_init), qlogis(p_nu_init), log(m_init), log(k_init))
    fn <- function(par) -ll_trans(par, p1, p2, sf1, sf2, z, n, phi)
    k <- 4
  } else if (model_name == "cis_trans") {
    init <- c(qlogis(p_mu_init), qlogis(p_nu_init), log(m_init), log(alpha_init), log(beta_init))
    fn <- function(par) -ll_cis_trans(par, p1, p2, sf1, sf2, z, n, phi)
    k <- 5
  } else {
    stop("Unknown model: ", model_name)
  }

  fit <- try(optim(par = init, fn = fn, method = optim_method, control = list(maxit = optim_maxit)), silent = TRUE)
  if (inherits(fit, "try-error") || !is.finite(fit$value)) {return(list(ok = FALSE, logLik = NA_real_, BIC = NA_real_, k = k, par = rep(NA_real_, k)))}
  n_obs <- length(p1) + length(p2) + length(z)
  ll <- -fit$value
  bic <- calc_bic(ll, k, n_obs)
  list(ok = TRUE, logLik = ll, BIC = bic, k = k, par = fit$par)
}

# --------------------------------------------------------------------------------------------------
# Decode fitted parameters for output
# --------------------------------------------------------------------------------------------------
decode_params <- function(best_model, fit_obj) {
  if (!isTRUE(fit_obj$ok)) {
    return(list(
      p_mu = NA_real_,
      p_nu = NA_real_,
      parent_odds_mu = NA_real_,
      parent_odds_nu = NA_real_,
      alpha = NA_real_,
      beta = NA_real_,
      hybrid_mu = NA_real_
    ))
  }

  par <- fit_obj$par

  if (best_model == "non_divergent") {
    p_mu <- plogis(par[1])
    p_nu <- p_mu
    alpha <- safe_exp(par[3])
    beta <- alpha
  } else if (best_model == "cis") {
    p_mu <- plogis(par[1])
    p_nu <- plogis(par[2])
    c0 <- safe_exp(par[4])
    alpha <- c0 * odds_from_p(p_mu)
    beta <- c0 * odds_from_p(p_nu)
  } else if (best_model == "trans") {
    p_mu <- plogis(par[1])
    p_nu <- plogis(par[2])
    alpha <- safe_exp(par[4])
    beta <- alpha
  } else if (best_model == "cis_trans") {
    p_mu <- plogis(par[1])
    p_nu <- plogis(par[2])
    alpha <- safe_exp(par[4])
    beta <- safe_exp(par[5])
  } else {
    return(list(p_mu = NA_real_, p_nu = NA_real_, parent_odds_mu = NA_real_, parent_odds_nu = NA_real_, alpha = NA_real_, beta = NA_real_, hybrid_mu = NA_real_))
  }

  od_mu <- odds_from_p(p_mu)
  od_nu <- odds_from_p(p_nu)
  hybrid_mu <- alpha / (alpha + beta)
  list(p_mu = p_mu, p_nu = p_nu, parent_odds_mu = od_mu, parent_odds_nu = od_nu, alpha = alpha, beta = beta, hybrid_mu = hybrid_mu)
}

# --------------------------------------------------------------------------------------------------
# Summary effects for cis_trans subclassification
# x = F1 average log2 allelic ratio
# y = F0 average log2 fold change
# --------------------------------------------------------------------------------------------------
calc_parent_summary_log2fc <- function(p1, p2, lib1, lib2, pseudocount = 0.5) {
  # CPM normalization
  cpm1 <- p1 * 1e6 / lib1
  cpm2 <- p2 * 1e6 / lib2

  # replicate-wise log2FC, then simple average
  mean(log2(cpm1 + pseudocount) - log2(cpm2 + pseudocount))
}

calc_hybrid_summary_log2fc <- function(h1, h2, pseudocount = 0.5) {
  # replicate-wise allele log2 ratio, then simple average
  mean(log2(h1 + pseudocount) - log2(h2 + pseudocount))
}

classify_cis_trans_subtype_simple <- function(x_f1, y_f0, eps = 1e-12) {
  # x = F1 = cis-like
  # y = F0 = overall-like

  if (is.na(x_f1) || is.na(y_f0)) return(NA_character_)
  if (abs(x_f1) < eps || abs(y_f0) < eps) return("ambiguous")

  # 1) CIS-trans
  if ((x_f1 * y_f0 > 0) && (abs(x_f1) > abs(y_f0))) {
    return("CIS-trans")
  }

  # 2) TRANS-cis
  if (x_f1 * y_f0 < 0) {
    return("TRANS-cis")
  }

  # 3) CIS+trans
  if ((x_f1 * y_f0 > 0) && (abs(x_f1) < abs(y_f0)) && (abs(y_f0) < 2 * abs(x_f1))) {
    return("CIS+trans")
  }

  # 4) TRANS+cis
  if ((x_f1 * y_f0 > 0) && (abs(y_f0) > 2 * abs(x_f1))) {
    return("TRANS+cis")
  }

  return("ambiguous")
}

# --------------------------------------------------------------------------------------------------
# Read input
# --------------------------------------------------------------------------------------------------
if (!file.exists(in_file)) {
  stop("Input file not found: ", in_file)
}

dt <- fread(in_file)
message("Loaded input table with ", nrow(dt), " peaks and ", ncol(dt), " columns.")
message("Input file: ", in_file)

coord_cols <- c(
  "peak_id",
  paste0(hap1_name, c("_chr", "_start", "_end")),
  paste0(hap2_name, c("_chr", "_start", "_end"))
)

missing_coord_cols <- setdiff(coord_cols, names(dt))
if (length(missing_coord_cols) > 0) {
  stop("Missing required coordinate columns: ", paste(missing_coord_cols, collapse = ", "))
}

sample_cols <- setdiff(names(dt), coord_cols)

# --------------------------------------------------------------------------------------------------
# Detect reciprocal hybrid names automatically
# --------------------------------------------------------------------------------------------------
hyb_h1_cols_all <- grep(paste0("_on", hap1_name, "$"), sample_cols, value = TRUE)
hyb_h2_cols_all <- grep(paste0("_on", hap2_name, "$"), sample_cols, value = TRUE)

cross_h1 <- unique(sub("_rep.*", "", hyb_h1_cols_all))
cross_h2 <- unique(sub("_rep.*", "", hyb_h2_cols_all))
crosses <- intersect(cross_h1, cross_h2)

if (length(crosses) == 0) {
  stop("No reciprocal hybrid columns detected.")
}
message("Detected reciprocal hybrid crosses: ", paste(crosses, collapse = ", "))

# --------------------------------------------------------------------------------------------------
# Identify columns for this cell type
# --------------------------------------------------------------------------------------------------
p1_cols <- grep(paste0("^", hap1_name, "_rep\\d+-C", ct, "$"), sample_cols, value = TRUE)
p2_cols <- grep(paste0("^", hap2_name, "_rep\\d+-C", ct, "$"), sample_cols, value = TRUE)

h1_cols <- unlist(lapply(crosses, function(x) {
  grep(paste0("^", x, "_rep\\d+-C", ct, "_on", hap1_name, "$"), sample_cols, value = TRUE)
}))

h2_cols <- unlist(lapply(crosses, function(x) {
  grep(paste0("^", x, "_rep\\d+-C", ct, "_on", hap2_name, "$"), sample_cols, value = TRUE)
}))

if (length(p1_cols) == 0 || length(p2_cols) == 0 || length(h1_cols) == 0 || length(h2_cols) == 0) {
  stop("Missing columns for C", ct)
}

# pair hybrid columns explicitly
h1_key <- sub(paste0("_on", hap1_name, "$"), "", h1_cols)
h2_key <- sub(paste0("_on", hap2_name, "$"), "", h2_cols)

if (!setequal(h1_key, h2_key)) {
  stop(
    "Unmatched hybrid columns detected for C", ct,
    ": H1-only = ", paste(setdiff(h1_key, h2_key), collapse = ", "),
    "; H2-only = ", paste(setdiff(h2_key, h1_key), collapse = ", ")
  )
}

common_key <- sort(intersect(h1_key, h2_key))
h1_cols <- h1_cols[match(common_key, h1_key)]
h2_cols <- h2_cols[match(common_key, h2_key)]

stopifnot(
  identical(
    sub(paste0("_on", hap1_name, "$"), "", h1_cols),
    sub(paste0("_on", hap2_name, "$"), "", h2_cols)
  )
)

message("---- Column check for C", ct, " ----")
message("Parent ", hap1_name, ": ", paste(p1_cols, collapse = ", "))
message("Parent ", hap2_name, ": ", paste(p2_cols, collapse = ", "))
message("Hybrid ", hap1_name, " haplotype: ", paste(h1_cols, collapse = ", "))
message("Hybrid ", hap2_name, " haplotype: ", paste(h2_cols, collapse = ", "))

# --------------------------------------------------------------------------------------------------
# Build raw matrices
# --------------------------------------------------------------------------------------------------
P1_raw <- as.matrix(dt[, ..p1_cols]); mode(P1_raw) <- "numeric"
P2_raw <- as.matrix(dt[, ..p2_cols]); mode(P2_raw) <- "numeric"
H1_raw <- as.matrix(dt[, ..h1_cols]); mode(H1_raw) <- "numeric"
H2_raw <- as.matrix(dt[, ..h2_cols]); mode(H2_raw) <- "numeric"

# --------------------------------------------------------------------------------------------------
# Filtering
# --------------------------------------------------------------------------------------------------
parent_total <- rowSums(P1_raw) + rowSums(P2_raw)
hybrid_total <- rowSums(H1_raw) + rowSums(H2_raw)
hyb_rep_total <- H1_raw + H2_raw
hyb_rep_support_n <- rowSums(hyb_rep_total >= min_hybrid_rep_count)

keep <- parent_total >= min_parent_total &
  hybrid_total >= min_hybrid_total &
  hyb_rep_support_n >= min_hybrid_supported_reps

message("Filtering summary for C", ct, ":")
message("  total peaks before filtering = ", nrow(dt))
message("  parent_total >= ", min_parent_total, " : ", sum(parent_total >= min_parent_total))
message("  hybrid_total >= ", min_hybrid_total, " : ", sum(hybrid_total >= min_hybrid_total))
message("  hybrid replicate support pass : ", sum(hyb_rep_support_n >= min_hybrid_supported_reps))
message("  passing all filters = ", sum(keep))

if (sum(keep) == 0) {
  stop("No peaks passed filtering for C", ct)
}

sub_dt <- copy(dt[keep, ..coord_cols])

P1_raw <- P1_raw[keep, , drop = FALSE]
P2_raw <- P2_raw[keep, , drop = FALSE]
H1_raw <- H1_raw[keep, , drop = FALSE]
H2_raw <- H2_raw[keep, , drop = FALSE]

# --------------------------------------------------------------------------------------------------
# Estimate library information and dispersion from parent counts using edgeR
# --------------------------------------------------------------------------------------------------
parent_counts <- cbind(P1_raw, P2_raw)
colnames(parent_counts) <- c(p1_cols, p2_cols)

parent_group <- factor(c(rep(hap1_name, ncol(P1_raw)), rep(hap2_name, ncol(P2_raw))))
dge <- DGEList(counts = parent_counts, group = parent_group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)

# Effective library sizes for the joint model
sf_all <- dge$samples$lib.size * dge$samples$norm.factors
sf_all <- sf_all / exp(mean(log(sf_all)))  # scale geometric mean to ~1

sf1 <- sf_all[seq_len(ncol(P1_raw))]
sf2 <- sf_all[ncol(P1_raw) + seq_len(ncol(P2_raw))]

# Raw library sizes for CPM in F0 summary calculation
lib1 <- dge$samples$lib.size[seq_len(ncol(P1_raw))]
lib2 <- dge$samples$lib.size[ncol(P1_raw) + seq_len(ncol(P2_raw))]

# tagwise dispersion, fallback to common if needed
phi_vec <- dge$tagwise.dispersion
if (is.null(phi_vec) || length(phi_vec) != nrow(P1_raw)) {
  phi_vec <- rep(dge$common.dispersion, nrow(P1_raw))
}

message("Parent size factors estimated by edgeR TMM for joint model.")
message("Parent raw library sizes will be used for CPM-based F0 summary log2FC.")
message("Parent dispersion: using tagwise dispersion with common-dispersion fallback.")

# --------------------------------------------------------------------------------------------------
# Joint fitting
# --------------------------------------------------------------------------------------------------
n_peaks <- nrow(P1_raw)
res_list <- vector("list", n_peaks)
fit_start_time <- Sys.time()

for (i in seq_len(n_peaks)) {
  p1 <- as.numeric(P1_raw[i, ])
  p2 <- as.numeric(P2_raw[i, ])

  phi <- phi_vec[i]
  if (!is.finite(phi) || is.na(phi) || phi <= 0) {
    phi <- dge$common.dispersion
  }
  if (!is.finite(phi) || is.na(phi) || phi <= 0) {
    phi <- 0.1
  }

  z_all <- as.numeric(H1_raw[i, ])
  n_all <- as.numeric(H1_raw[i, ] + H2_raw[i, ])

  keep_h <- n_all >= min_hybrid_rep_count
  z <- z_all[keep_h]
  n <- n_all[keep_h]

  if (length(z) < min_hybrid_supported_reps || sum(n) == 0) {
    res_list[[i]] <- data.table(
      best_model = NA_character_,
      BIC_non_divergent = NA_real_,
      BIC_cis = NA_real_,
      BIC_trans = NA_real_,
      BIC_cis_trans = NA_real_,
      deltaBIC = NA_real_,
      p_mu = NA_real_,
      p_nu = NA_real_,
      parent_odds_mu = NA_real_,
      parent_odds_nu = NA_real_,
      alpha = NA_real_,
      beta = NA_real_,
      hybrid_mu = NA_real_,
      F0_summary_log2FC = NA_real_,
      F1_summary_log2FC = NA_real_,
      trans_summary_log2FC = NA_real_,
      cis_trans_subtype = NA_character_,
      parent_dispersion = phi
    )
    next
  }

  fit_nodiv <- fit_one_model("non_divergent", p1, p2, sf1, sf2, z, n, phi)
  fit_cis <- fit_one_model("cis", p1, p2, sf1, sf2, z, n, phi)
  fit_trans <- fit_one_model("trans", p1, p2, sf1, sf2, z, n, phi)
  fit_ct <- fit_one_model("cis_trans", p1, p2, sf1, sf2, z, n, phi)
  bic_vec <- c(non_divergent = fit_nodiv$BIC, cis = fit_cis$BIC, trans = fit_trans$BIC, cis_trans = fit_ct$BIC)

  if (all(is.na(bic_vec))) {
    res_list[[i]] <- data.table(
      best_model = NA_character_,
      BIC_non_divergent = NA_real_,
      BIC_cis = NA_real_,
      BIC_trans = NA_real_,
      BIC_cis_trans = NA_real_,
      deltaBIC = NA_real_,
      p_mu = NA_real_,
      p_nu = NA_real_,
      parent_odds_mu = NA_real_,
      parent_odds_nu = NA_real_,
      alpha = NA_real_,
      beta = NA_real_,
      hybrid_mu = NA_real_,
      F0_summary_log2FC = NA_real_,
      F1_summary_log2FC = NA_real_,
      trans_summary_log2FC = NA_real_,
      cis_trans_subtype = NA_character_,
      parent_dispersion = phi
    )
    next
  }

  best_model <- names(which.min(bic_vec))
  sorted_bic <- sort(bic_vec, na.last = NA)
  deltaBIC <- if (length(sorted_bic) >= 2) sorted_bic[2] - sorted_bic[1] else NA_real_
  best_fit <- switch(best_model, non_divergent = fit_nodiv, cis = fit_cis, trans = fit_trans, cis_trans = fit_ct)
  dec <- decode_params(best_model, best_fit)

  # summary effects used for cis_trans subclassification
  F0_summary_log2FC <- calc_parent_summary_log2fc(p1 = p1, p2 = p2, lib1 = lib1, lib2 = lib2, pseudocount = summary_pseudocount)
  F1_summary_log2FC <- calc_hybrid_summary_log2fc(h1 = as.numeric(H1_raw[i, ]), h2 = as.numeric(H2_raw[i, ]), pseudocount = summary_pseudocount)
  trans_summary_log2FC <- F0_summary_log2FC - F1_summary_log2FC
  
  subtype <- if (best_model == "cis_trans") {
    classify_cis_trans_subtype_simple(
      x_f1 = F1_summary_log2FC,
      y_f0 = F0_summary_log2FC
    )
  } else {
    NA_character_
  }

  res_list[[i]] <- data.table(
    best_model = best_model,
    BIC_non_divergent = bic_vec["non_divergent"],
    BIC_cis = bic_vec["cis"],
    BIC_trans = bic_vec["trans"],
    BIC_cis_trans = bic_vec["cis_trans"],
    deltaBIC = deltaBIC,
    p_mu = dec$p_mu,
    p_nu = dec$p_nu,
    parent_odds_mu = dec$parent_odds_mu,
    parent_odds_nu = dec$parent_odds_nu,
    alpha = dec$alpha,
    beta = dec$beta,
    hybrid_mu = dec$hybrid_mu,
    F0_summary_log2FC = F0_summary_log2FC,
    F1_summary_log2FC = F1_summary_log2FC,
    trans_summary_log2FC = trans_summary_log2FC,
    cis_trans_subtype = subtype,
    parent_dispersion = phi
  )

  if (i %% 1000 == 0 || i == n_peaks) {
    elapsed_min <- round(as.numeric(difftime(Sys.time(), fit_start_time, units = "mins")), 2)
    message("C", ct, ": processed ", i, " / ", n_peaks, " peaks; elapsed = ", elapsed_min, " min")
  }
}

res_dt <- rbindlist(res_list, fill = TRUE)

# --------------------------------------------------------------------------------------------------
# Build final output
# --------------------------------------------------------------------------------------------------
count_dt <- cbind(as.data.table(P1_raw), as.data.table(P2_raw), as.data.table(H1_raw), as.data.table(H2_raw))
setnames(count_dt, c(p1_cols, p2_cols, h1_cols, h2_cols))

final_dt <- cbind(sub_dt, data.table(cell_type = paste0("C", ct)), res_dt, count_dt)
out_file <- file.path(final_dir, paste0(hap1_name, "-", hap2_name, "_C", ct, "_joint_model.tsv"))
fwrite(final_dt, out_file, sep = "\t", quote = FALSE)

# cis-trans analysis summary
message("Cell type C", ct, " finished.")
message("Output: ", out_file)
message("Model summary:")
print(table(final_dt$best_model, useNA = "ifany"))
message("cis_trans subtype summary:")
print(table(final_dt$cis_trans_subtype, useNA = "ifany"))


# --------------------------------------------------------------------------------------------------
# Append extra information, including genomic context, tau, sequence/accessibility conservation info
# --------------------------------------------------------------------------------------------------
info_dt <- fread(info_file)

append_col_idx <- c(5, (ncol(info_dt) - 10):ncol(info_dt))
append_col_names <- names(info_dt)[append_col_idx]
info_sub <- info_dt[, c("cell_type", "peak_id", append_col_names), with = FALSE]
info_ct <- info_sub[cell_type == paste0("C", ct)]
info_ct <- unique(info_ct, by = "peak_id")
final_dt[, row_id__tmp := .I]
final_with_info_dt <- merge(final_dt, info_ct[, !("cell_type"), with = FALSE], by = "peak_id", all.x = TRUE, sort = FALSE)

setorder(final_with_info_dt, row_id__tmp)
final_with_info_dt[, row_id__tmp := NULL]

original_cols <- setdiff(names(final_dt), "row_id__tmp")
extra_cols <- setdiff(names(final_with_info_dt), c(original_cols, "row_id__tmp"))
final_cols <- c(original_cols, extra_cols)
final_with_info_dt <- final_with_info_dt[, ..final_cols]

out_file_with_info <- file.path(final_dir, paste0(hap1_name, "-", hap2_name, "_C", ct, "_joint_model_with_extra_info.tsv"))
fwrite(final_with_info_dt, out_file_with_info, sep = "\t", quote = FALSE, na = "NA")
message("Output with extra info: ", out_file_with_info)


