# =============================================================================
# heckman.R
# Correction (individual-level): Heckman two-step selection correction
#
# Step 1: probit model of selection, fitted on everyone in the sampling frame.
# Step 2: outcome regression among participants, including the inverse Mills
#         ratio from step 1, which absorbs the effect of selection on
#         unmeasured factors correlated with the outcome.
#
# Requirements: covariates for non-participants (within-cohort nonresponse);
#               an exclusion variable (E) that predicts selection but not Y.
# Predictors of selection that also affect Y (Z) must be in the outcome model,
# so estimates are DIRECT SNP effects conditional on Z (compare to truth_direct).
#
# If SNPs affect selection (e.g. through Y, module a2), the selection equation
# must include the SNP (per_snp = TRUE). Fitting it once without genotype
# (per_snp = FALSE) is valid only if SNP effects on selection are negligible.
#
# Module: (a). References: Heckman (1979); Puhani (2000); Mignogna et al. (2023).
# Software: sampleSelection (two-step estimates with corrected standard errors).
# =============================================================================

source("1_functions/gwas_helpers.R")

#' Heckman-corrected GWAS
#'
#' @param S 0/1 selection indicator for everyone in the sampling frame.
#' @param Y Outcome (values for non-participants are ignored).
#' @param G Genotype matrix for everyone in the sampling frame.
#' @param sel_covar Data frame of selection predictors (must include the
#'   exclusion variable).
#' @param out_covar Data frame of outcome-model covariates (selection
#'   predictors that also affect Y); may be NULL.
#' @param per_snp Include the SNP in the selection equation.
#' @return Data frame: snp, beta, se, z, p, n, plus attribute "imr_coef"
#'   (mean inverse Mills ratio coefficient across SNPs).
heckman_gwas <- function(S, Y, G, sel_covar, out_covar = NULL, per_snp = TRUE) {
  if (!requireNamespace("sampleSelection", quietly = TRUE))
    stop("Install the sampleSelection package (see 0_install/install.R).")

  Y_fill <- ifelse(S == 1, Y, 0)       # outcome is only used where S == 1
  base   <- data.frame(S = S, Y = Y_fill, sel_covar)
  if (!is.null(out_covar)) base <- cbind(base, out_covar)

  sel_terms <- names(sel_covar)
  out_terms <- if (is.null(out_covar)) character(0) else names(out_covar)

  sel_f <- stats::reformulate(c(if (per_snp) "g", sel_terms), response = "S")
  out_f <- stats::reformulate(c("g", out_terms), response = "Y")

  est <- t(sapply(seq_len(ncol(G)), function(j) {
    dat   <- base
    dat$g <- G[, j]
    if (stats::var(dat$g[S == 1]) == 0) return(c(NA, NA, NA))
    fit <- tryCatch(
      sampleSelection::selection(sel_f, out_f, data = dat, method = "2step"),
      error = function(e) NULL)
    if (is.null(fit)) return(c(NA, NA, NA))
    b   <- stats::coef(fit, part = "outcome")
    V   <- stats::vcov(fit, part = "outcome")
    imr <- if ("invMillsRatio" %in% names(b)) b[["invMillsRatio"]] else NA
    c(b[["g"]], sqrt(V["g", "g"]), imr)
  }))

  out <- data.frame(snp = colnames(G), beta = est[, 1], se = est[, 2],
                    row.names = NULL)
  out$z <- out$beta / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  out$n <- sum(S == 1)
  attr(out, "imr_coef") <- mean(est[, 3], na.rm = TRUE)
  out
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0) {
  d <- readRDS("data/module_a.rds")
  p <- d$pheno

  results <- list()
  for (sc in c("S_a1", "S_a2")) {
    sel <- p[[sc]] == 1
    cat("\n====", sc, "====\n")

    # Comparison: naive regression among participants, adjusted for Z
    naive_adjZ <- run_gwas(p$Y[sel], d$G[sel, ], covar = p[sel, "Z", drop = FALSE])

    # Heckman with and without the SNP in the selection equation
    hk_snp  <- heckman_gwas(p[[sc]], p$Y, d$G,
                            sel_covar = p[, c("Z", "E")],
                            out_covar = p[, "Z", drop = FALSE], per_snp = TRUE)
    hk_once <- heckman_gwas(p[[sc]], p$Y, d$G,
                            sel_covar = p[, c("Z", "E")],
                            out_covar = p[, "Z", drop = FALSE], per_snp = FALSE)

    cat(sprintf("Mean inverse Mills ratio coefficient: %.3f (0 = no selection on unobservables)\n",
                attr(hk_snp, "imr_coef")))

    results[[sc]] <- rbind(
      compare_by_group(naive_adjZ, d$truth_direct, d$groups, paste(sc, "naive (adj. Z)")),
      compare_by_group(hk_once,    d$truth_direct, d$groups, paste(sc, "Heckman, SNP not in selection")),
      compare_by_group(hk_snp,     d$truth_direct, d$groups, paste(sc, "Heckman, SNP in selection"))
    )
  }

  out <- do.call(rbind, results)
  print(out[, c("analysis", "snp_group", "bias_causal", "rel_bias_causal",
                "mean_abs_bias", "coverage", "type1_error")], row.names = FALSE)
}
