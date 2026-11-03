# =============================================================================
# gwas_helpers.R
# Shared GWAS and evaluation functions used by every module and method.
# Contains functions only; no code runs when the file is sourced.
# =============================================================================

#' Per-SNP GWAS (linear or logistic), optionally weighted
#'
#' Fits y ~ SNP + covariates separately for each SNP. With weights (e.g. IPW),
#' robust (sandwich, HC0) standard errors are used by default.
#'
#' @param y Outcome vector (continuous, or 0/1 for family = "binomial").
#' @param G Genotype matrix (n x m) with column names.
#' @param covar Optional covariate matrix or data frame (n rows).
#' @param weights Optional observation weights (length n).
#' @param family "gaussian" (linear) or "binomial" (logistic).
#' @param robust Use sandwich standard errors (default TRUE when weighted).
#' @return Data frame: snp, beta, se, z, p, n.
run_gwas <- function(y, G, covar = NULL, weights = NULL,
                     family = c("gaussian", "binomial"),
                     robust = !is.null(weights)) {
  family <- match.arg(family)
  n <- length(y)
  stopifnot(nrow(G) == n)
  
  C <- if (is.null(covar)) NULL else as.matrix(covar)
  if (!is.null(C)) stopifnot(nrow(C) == n)
  w <- if (is.null(weights)) rep(1, n) else weights
  stopifnot(length(w) == n, all(w > 0, na.rm = TRUE))
  
  keep <- !is.na(y) & !is.na(w)
  if (!is.null(C)) keep <- keep & stats::complete.cases(C)
  yk <- y[keep]
  wk <- w[keep]
  Ck <- if (is.null(C)) NULL else C[keep, , drop = FALSE]
  
  fit_one <- function(g) {
    if (stats::var(g) == 0) return(c(NA_real_, NA_real_))
    X <- cbind(1, g, Ck)
    
    if (family == "gaussian") {
      fit <- stats::lm.wfit(X, yk, wk)
      if (fit$rank < ncol(X)) return(c(NA_real_, NA_real_))
      e     <- fit$residuals
      bread <- chol2inv(qr.R(fit$qr))          # (X'WX)^-1
      if (robust) {
        meat <- crossprod(X * (wk * e))         # sum_i w_i^2 e_i^2 x_i x_i'
        V    <- bread %*% meat %*% bread
      } else {
        sigma2 <- sum(wk * e^2) / (length(yk) - ncol(X))
        V      <- bread * sigma2
      }
      b <- fit$coefficients[2]
      
    } else {
      fit <- suppressWarnings(stats::glm.fit(X, yk, weights = wk,
                                             family = stats::binomial()))
      if (!fit$converged || fit$rank < ncol(X)) return(c(NA_real_, NA_real_))
      bread <- chol2inv(qr.R(fit$qr))          # (X' W_working X)^-1
      if (robust) {
        meat <- crossprod(X * (wk * (yk - fit$fitted.values)))
        V    <- bread %*% meat %*% bread
      } else {
        V <- bread
      }
      b <- fit$coefficients[2]
    }
    c(b, sqrt(V[2, 2]))
  }
  
  est <- t(apply(G[keep, , drop = FALSE], 2, fit_one))
  out <- data.frame(
    snp  = colnames(G),
    beta = est[, 1],
    se   = est[, 2],
    row.names = NULL
  )
  out$z <- out$beta / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))     # normal approximation (large n)
  out$n <- sum(keep)
  out
}

#' Genomic inflation factor (lambda)
#'
#' @param p Vector of association p-values.
#' @return Median chi-square statistic divided by its expected value under
#'   the null (0.4549 for 1 df).
genomic_lambda <- function(p) {
  p <- p[!is.na(p)]
  stats::median(stats::qchisq(p, df = 1, lower.tail = FALSE)) /
    stats::qchisq(0.5, df = 1)
}

#' Compare GWAS estimates with the known true SNP effects
#'
#' @param gwas Output of run_gwas() (or any data frame with beta, se, p).
#' @param truth Vector of true effects, same order as gwas rows.
#' @param label Optional name for the analysis (e.g. "naive", "IPW").
#' @param alpha Significance level for power and type I error.
#' @return One-row data frame of performance metrics.
compare_to_truth <- function(gwas, truth, label = NA_character_, alpha = 0.05) {
  stopifnot(nrow(gwas) == length(truth))
  causal <- truth != 0
  lo <- gwas$beta - stats::qnorm(0.975) * gwas$se
  hi <- gwas$beta + stats::qnorm(0.975) * gwas$se
  
  data.frame(
    analysis      = label,
    bias_causal   = mean((gwas$beta - truth)[causal], na.rm = TRUE),
    mean_abs_bias = mean(abs(gwas$beta - truth), na.rm = TRUE),
    rmse          = sqrt(mean((gwas$beta - truth)^2, na.rm = TRUE)),
    coverage      = mean(lo <= truth & truth <= hi, na.rm = TRUE),
    power         = mean(gwas$p[causal]  < alpha, na.rm = TRUE),
    type1_error   = mean(gwas$p[!causal] < alpha, na.rm = TRUE),
    lambda_null   = genomic_lambda(gwas$p[!causal])
  )
}
