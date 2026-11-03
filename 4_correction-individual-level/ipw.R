# =============================================================================
# ipw.R
# Correction (individual-level): inverse probability weighting (IPW)
#
# Part 1 (module a): non-representative participation
#   1. A selection model predicts each person's probability of being in the
#      sample from measured covariates (fitted where non-participants'
#      covariates are known: a sampling frame, or a stacked reference sample).
#   2. Participants are weighted by the (stabilised) inverse of that probability.
#   3. The GWAS is run with these weights and robust standard errors.
#   Corrects selection on measured covariates (a1), not on the outcome (a2).
#
# Part 2 (module b): case-control secondary phenotypes
#   Sampling probabilities are known from the design (cases over-sampled), so
#   weights are population share / sample share. With one secondary phenotype
#   this is a weighted GWAS; with several, the scaled marginal model (SMM) fits
#   a common SNP effect on standardised phenotypes by weighted GEE.
#
# References: Cole and Hernan (2008); Schoeler et al. (2023);
#             van Alten et al. (2024); Schifano et al. (2013).
# =============================================================================

source("1_functions/gwas_helpers.R")

# ---- Part 1: participation weights -------------------------------------------

#' Estimate inverse probability weights from a selection model
#'
#' @param S 0/1 selection indicator for everyone in the sampling frame.
#' @param covar Data frame of selection predictors (same rows as S).
#' @param link "probit" or "logit".
#' @param stabilise Multiply by the marginal selection probability.
#' @param truncate Optional quantiles to truncate weights at, e.g. c(0.01, 0.99).
#' @return Weights for selected individuals (S == 1), with attributes
#'   "p_hat" (all individuals), "ess" (effective sample size), "model".
ipw_weights <- function(S, covar, link = c("probit", "logit"),
                        stabilise = TRUE, truncate = NULL) {
  link  <- match.arg(link)
  dat   <- data.frame(S = S, covar)
  fit   <- stats::glm(S ~ ., data = dat, family = stats::binomial(link = link))
  p_hat <- stats::fitted(fit)

  w <- 1 / p_hat[S == 1]
  if (stabilise) w <- w * mean(S)
  if (!is.null(truncate)) {
    q <- stats::quantile(w, truncate)
    w <- pmin(pmax(w, q[1]), q[2])
  }

  attr(w, "p_hat") <- p_hat
  attr(w, "ess")   <- sum(w)^2 / sum(w^2)
  attr(w, "model") <- fit
  w
}

#' IPW-weighted GWAS
#'
#' @param y Outcome for the analysed sample.
#' @param G Genotypes for the analysed sample.
#' @param w Weights (from ipw_weights() or case-control sampling weights).
#' @param covar Optional outcome-model covariates.
#' @return run_gwas() output with robust standard errors.
ipw_gwas <- function(y, G, w, covar = NULL) {
  run_gwas(y, G, covar = covar, weights = as.numeric(w), robust = TRUE)
}

# ---- Part 2: case-control secondary phenotypes ------------------------------

#' Case-control sampling weights
#'
#' @param D 0/1 disease status in the case-control sample.
#' @param K Population prevalence of D.
#' @return Weight per individual: population share / sample share.
case_control_weights <- function(D, K) {
  share_case <- mean(D)
  ifelse(D == 1, K / share_case, (1 - K) / (1 - share_case))
}

#' Scaled marginal model (IPW-weighted GEE) for several secondary phenotypes
#'
#' Each phenotype is standardised using weighted (population-representative)
#' means and SDs; a common SNP effect alpha is estimated across phenotypes,
#' with phenotype-specific intercepts, weighted GEE and robust SEs.
#'
#' @param Ymat Matrix or data frame of secondary phenotypes (n x M).
#' @param G Genotype matrix (n x m).
#' @param w Sampling weights (length n).
#' @param corstr GEE working correlation ("exchangeable", "independence").
#' @return Data frame: snp, beta (common effect on the standardised scale),
#'   se, z, p, n.
smm_gwas <- function(Ymat, G, w, corstr = c("exchangeable", "independence")) {
  corstr <- match.arg(corstr)
  if (!requireNamespace("geepack", quietly = TRUE))
    stop("Install the geepack package (see 0_install/install.R).")

  Ymat <- as.matrix(Ymat)
  n <- nrow(Ymat)
  M <- ncol(Ymat)

  # Weighted standardisation (population-representative scale)
  Ys <- apply(Ymat, 2, function(y) {
    mu <- sum(w * y) / sum(w)
    sdv <- sqrt(sum(w * (y - mu)^2) / sum(w))
    (y - mu) / sdv
  })

  long <- data.frame(
    id    = rep(seq_len(n), times = M),
    pheno = factor(rep(seq_len(M), each = n)),
    y     = as.numeric(Ys),
    w     = rep(w, times = M)
  )
  long <- long[order(long$id, long$pheno), ]
  ord  <- long$id

  est <- t(apply(G, 2, function(g) {
    long$g <- g[ord]
    fit <- tryCatch(
      geepack::geeglm(y ~ 0 + pheno + g, id = id, weights = w, data = long,
                      family = stats::gaussian(), corstr = corstr),
      error = function(e) NULL)
    if (is.null(fit)) return(c(NA, NA))
    cf <- summary(fit)$coefficients
    c(cf["g", "Estimate"], cf["g", "Std.err"])
  }))

  out <- data.frame(snp = colnames(G), beta = est[, 1], se = est[, 2],
                    row.names = NULL)
  out$z <- out$beta / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  out$n <- n
  out
}

# ---- Examples (run only when this file is run directly) --------------------
if (sys.nframe() == 0) {

  # Module (a): participation
  if (file.exists("data/module_a.rds")) {
    d <- readRDS("data/module_a.rds")
    p <- d$pheno
    res_a <- list()
    for (sc in c("S_a1", "S_a2")) {
      sel <- p[[sc]] == 1
      w   <- ipw_weights(p[[sc]], covar = p[, c("Z", "E")])
      cat(sprintf("\nModule a, %s: participants %d, effective sample size %.0f\n",
                  sc, sum(sel), attr(w, "ess")))
      res_a[[sc]] <- rbind(
        compare_by_group(run_gwas(p$Y[sel], d$G[sel, ]), d$truth, d$groups,
                         paste(sc, "naive")),
        compare_by_group(ipw_gwas(p$Y[sel], d$G[sel, ], w), d$truth, d$groups,
                         paste(sc, "IPW"))
      )
    }
    print(do.call(rbind, res_a)[, c("analysis", "snp_group", "rel_bias_causal",
                                     "coverage", "power", "type1_error")],
          row.names = FALSE)
  }

  # Module (b): case-control secondary phenotypes
  if (file.exists("data/module_b.rds")) {
    d <- readRDS("data/module_b.rds")
    p <- d$pheno
    w <- case_control_weights(p$D, d$K)

    # Effects on the standardised scale (true SD of Y is ~1 in the simulation)
    naive <- run_gwas(p$Y1, d$G)
    ipw1  <- ipw_gwas(p$Y1, d$G, w)
    smm   <- smm_gwas(p[, c("Y1", "Y2")], d$G, w)

    cat("\nModule b: secondary phenotype analyses\n")
    print(rbind(
      compare_by_group(naive, d$truth, d$groups, "naive (Y1)"),
      compare_by_group(ipw1,  d$truth, d$groups, "IPW (Y1)"),
      compare_by_group(smm,   d$truth, d$groups, "SMM (Y1 + Y2)")
    )[, c("analysis", "snp_group", "rel_bias_causal", "mean_abs_bias",
          "coverage", "type1_error")], row.names = FALSE)
  }
}
