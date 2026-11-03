# =============================================================================
# proxy_outcomes.R
# Diagnostic: proxy (parental) outcomes for survival to recruitment
#
# Participants' genotypes are regressed on their parents' outcomes (GWAX).
# Parental outcomes are recorded whether or not the parent survived, so the
# proxy GWAS is not distorted by the participant's survival. Each participant
# carries half of a parent's alleles, so proxy effects are scaled by 2.
#
# Diagnostic: a SNP whose own-outcome estimate differs from its proxy estimate
# indicates survival (or other) selection on the participant's outcome.
# Proxy estimates can also replace the biased estimates, at the cost of wider
# standard errors.
#
# Assumes own and proxy estimates are independent given genotype (no shared
# environment between parent and offspring).
# Module: (e). References: Liu et al. (2017) Nat Genet (GWAX);
#         Mitchell et al. (2023) for selection on survival.
# =============================================================================

source("1_functions/gwas_helpers.R")

#' GWAS of a proxy (parental) outcome, rescaled to the participant scale
#' @param y_proxy Parental outcome (one parent, or the mean of both parents).
#' @param G Participants' genotypes.
#' @param scale Rescaling factor: 2 for one parent or the parental mean.
#' @param covar Optional covariates.
#' @return run_gwas() output with beta and se multiplied by `scale`.
proxy_gwas <- function(y_proxy, G, scale = 2, covar = NULL) {
  res <- run_gwas(y_proxy, G, covar = covar)
  res$beta <- res$beta * scale
  res$se   <- res$se * scale
  res
}

#' Per-SNP test of own-outcome vs proxy estimates
#' @param own,proxy run_gwas()-style data frames (same SNP order).
#' @return Data frame: snp, diff, se, z, p.
compare_own_proxy <- function(own, proxy) {
  stopifnot(identical(own$snp, proxy$snp))
  diff <- own$beta - proxy$beta
  se   <- sqrt(own$se^2 + proxy$se^2)
  data.frame(snp = own$snp, diff = diff, se = se, z = diff / se,
             p = 2 * stats::pnorm(-abs(diff / se)))
}

#' Summarise own-vs-proxy differences
#' @param cmp compare_own_proxy() output.
#' @param idx Optional subset of SNP rows.
#' @return One-row data frame: n_snps, lambda, prop_p_lt_0.05, omnibus_p
#'   (sum of z^2 against chi-square with n_snps df).
summarise_own_proxy <- function(cmp, idx = NULL) {
  if (!is.null(idx)) cmp <- cmp[idx, , drop = FALSE]
  z <- cmp$z[!is.na(cmp$z)]
  data.frame(n_snps = length(z),
             lambda = genomic_lambda(cmp$p[!is.na(cmp$p)]),
             prop_p_lt_0.05 = mean(cmp$p < 0.05, na.rm = TRUE),
             omnibus_p = stats::pchisq(sum(z^2), df = length(z),
                                       lower.tail = FALSE))
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0 && file.exists("data/module_e.rds")) {
  d <- readRDS("data/module_e.rds")
  p <- d$pheno
  alive <- p$R == 1

  full  <- run_gwas(p$Y, d$G)
  own   <- run_gwas(p$Y[alive], d$G[alive, ])
  proxy <- proxy_gwas(p$Y_P[alive], d$G[alive, ])

  cat("Estimates of SNP effects on Y\n")
  print(rbind(
    compare_by_group(full,  d$truth, d$groups, "full cohort (reference)"),
    compare_by_group(own,   d$truth, d$groups, "survivors, own outcome"),
    compare_by_group(proxy, d$truth, d$groups, "survivors, parental proxy")
  )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
        "power", "type1_error")], row.names = FALSE)

  cmp <- compare_own_proxy(own, proxy)
  cat("\nDiagnostic: own-outcome vs proxy estimates\n")
  print(rbind(
    cbind(snp_group = "all", summarise_own_proxy(cmp)),
    do.call(rbind, lapply(names(d$groups), function(g)
      cbind(snp_group = g, summarise_own_proxy(cmp, d$groups[[g]]))))
  ), row.names = FALSE)
}
