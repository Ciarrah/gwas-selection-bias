# =============================================================================
# sex_stratified.R
# Correction (summary-level): sex-stratified GWAS + fixed-effects meta-analysis
#
# GWAS are run separately in men and women and combined by inverse-variance
# weighted (IVW) fixed-effects meta-analysis, as in Lee et al. (2022).
# Removes bias from pooling sexes when participation depends on sex and is
# rare (logistic selection approximately multiplicative); does not remove
# within-sex collider bias when participation is common.
#
# Module: (d). Reference: Lee et al. (2022) Genome Biol;
#         code: github.com/hanbin973/autosomalsexgwas
# =============================================================================

source("1_functions/gwas_helpers.R")

#' Fixed-effects IVW meta-analysis of two sets of GWAS results
#' @param res1,res2 run_gwas()-style data frames (same SNP order).
#' @return Data frame: snp, beta, se, z, p, n, het_p (Cochran's Q, 1 df).
meta_ivw <- function(res1, res2) {
  stopifnot(identical(res1$snp, res2$snp))
  w1 <- 1 / res1$se^2
  w2 <- 1 / res2$se^2
  beta <- (w1 * res1$beta + w2 * res2$beta) / (w1 + w2)
  se   <- sqrt(1 / (w1 + w2))
  Q    <- w1 * (res1$beta - beta)^2 + w2 * (res2$beta - beta)^2
  out  <- data.frame(snp = res1$snp, beta = beta, se = se)
  out$z <- out$beta / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  out$n <- res1$n + res2$n
  out$het_p <- stats::pchisq(Q, df = 1, lower.tail = FALSE)
  out
}

#' Sex-stratified GWAS combined by fixed-effects meta-analysis
#' @param y Outcome.
#' @param G Genotype matrix.
#' @param sex 0/1 sex.
#' @param covar Optional covariates (applied within each sex).
#' @param family "gaussian" or "binomial".
#' @return meta_ivw() output, with the sex-specific results as attribute
#'   "strata" (list(female, male)).
sex_stratified_gwas <- function(y, G, sex, covar = NULL,
                                family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  fit <- function(k) {
    i <- sex == k
    cv <- if (is.null(covar)) NULL else as.matrix(covar)[i, , drop = FALSE]
    run_gwas(y[i], G[i, , drop = FALSE], covar = cv, family = family)
  }
  female <- fit(0)
  male   <- fit(1)
  out <- meta_ivw(female, male)
  attr(out, "strata") <- list(female = female, male = male)
  out
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0) {
  for (sc in c("d_low", "d_high")) {
    f <- sprintf("data/module_%s.rds", sc)
    if (!file.exists(f)) next
    d <- readRDS(f)
    p <- d$pheno

    pooled <- run_gwas(p$Y, d$G)
    adjsex <- run_gwas(p$Y, d$G, covar = p[, "Sex", drop = FALSE])
    strat  <- sex_stratified_gwas(p$Y, d$G, p$Sex)

    cat(sprintf("\n==== Module %s (%s) ====\n", sc, d$params$label))
    print(rbind(
      compare_by_group(pooled, d$truth, d$groups, "pooled"),
      compare_by_group(adjsex, d$truth, d$groups, "pooled, sex-adjusted"),
      compare_by_group(strat,  d$truth, d$groups, "sex-stratified meta")
    )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
          "type1_error")], row.names = FALSE)
  }
}
