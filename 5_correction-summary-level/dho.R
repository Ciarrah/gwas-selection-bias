# =============================================================================
# dho.R
# Correction (summary-level): Dudbridge-Hedges-Olkin (DHO) index event correction
#
# Uses the authors' package indexevent (DudbridgeLab/indexevent).
# Among cases, the SNP effect on progression is biased by -b * beta_X. b is the
# slope of beta_Y on beta_X across independent SNPs (weights 1/se_Y^2),
# corrected for regression dilution; adjusted effect = beta_Y - b * beta_X.
#
# method = "Hedges-Olkin" is DHO (Dudbridge et al. 2019). The package default,
# "CWLS" (Cai et al. 2022), is an alternative dilution correction.
# Assumes SNP effects on incidence and progression are uncorrelated across the
# SNPs used (violated in scenario c2).
# Note: The indexevent package now defaults to corrected weighted least squares (Cai et al. 2022);
# Hedges–Olkin is selected with method = “Hedges-Olkin”
# =============================================================================

source("1_functions/gwas_helpers.R")

#' DHO correction via indexevent
#'
#' @param snp SNP identifiers.
#' @param beta_x,se_x Incidence GWAS estimates and SEs (aligned to beta_y).
#' @param beta_y,se_y Progression GWAS estimates and SEs (cases).
#' @param prune Indices of approximately independent SNPs for the slope
#'   (NULL = all SNPs).
#' @param method "Hedges-Olkin" (DHO), "CWLS" or "Simex".
#' @param n Sample size to record.
#' @return Data frame: snp, beta, se, z, p, n (adjusted effects), with the
#'   indexevent object as attribute "fit".
dho_correct <- function(snp, beta_x, se_x, beta_y, se_y, prune = NULL,
                        method = "Hedges-Olkin", n = NA) {
  if (!requireNamespace("indexevent", quietly = TRUE))
    stop("Install indexevent: remotes::install_github('DudbridgeLab/indexevent')")

  fit <- indexevent::indexevent(beta_x, se_x, beta_y, se_y,
                                weighted = TRUE, prune = prune, method = method)

  out <- data.frame(snp = snp, beta = fit$ybeta.adj, se = fit$yse.adj)
  out$z <- out$beta / out$se
  out$p <- fit$yp.adj
  out$n <- n
  attr(out, "fit") <- fit
  out
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0) {
  for (sc in c("c1", "c2", "c3")) {
    f <- sprintf("data/module_%s.rds", sc)
    if (!file.exists(f)) next
    d <- readRDS(f)
    s <- d$sumstats

    adj <- dho_correct(s$snp, s$beta_x, s$se_x, s$beta_y, s$se_y,
                       n = d$n_cases)
    fit <- attr(adj, "fit")

    # Oracle: Hedges-Olkin slope fitted on the true X-only SNPs
    oracle <- indexevent::indexevent(s$beta_x, s$se_x, s$beta_y, s$se_y,
                                     prune = which(s$group == "X_only"),
                                     method = "Hedges-Olkin")$b

    cat(sprintf(paste0("\nModule %s (%s): DHO slope %.3f (SE %.3f; raw %.3f); ",
                       "oracle (X-only) %.3f\n"),
                sc, d$params$label, fit$b, fit$b.se, fit$b.raw, oracle))
    print(rbind(
      compare_by_group(d$gwas_y, d$truth_y, d$groups, "naive (cases)"),
      compare_by_group(adj,      d$truth_y, d$groups, "DHO")
    )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
          "type1_error")], row.names = FALSE)
  }
}
