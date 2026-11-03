# =============================================================================
# slope_hunter.R
# Correction (summary-level): Slope-Hunter index event correction
#
# Uses the authors' package SlopeHunter (Osmahmoud/SlopeHunter).
# Among SNPs associated with incidence (p_x <= xp_thresh), a two-component
# Gaussian mixture separates SNPs affecting incidence only ("Hunted") from
# pleiotropic SNPs. The slope of beta_Y on beta_X in the Hunted cluster
# estimates the collider bias b; adjusted effect = beta_Y - b * beta_X.
#
# Assumes the incidence-only SNPs contribute at least as much to the
# heritability of incidence as the pleiotropic SNPs. Unlike DHO, it allows
# correlated pleiotropy (scenario c2). No regression-dilution correction.
#
# Reference: Mahmoud et al. (2022) Nat Commun.
# =============================================================================

source("1_functions/gwas_helpers.R")

#' Slope-Hunter correction via SlopeHunter::hunt()
#'
#' @param snp SNP identifiers.
#' @param beta_x,se_x Incidence GWAS estimates and SEs (aligned to beta_y).
#' @param beta_y,se_y Progression GWAS estimates and SEs (cases).
#' @param xp_thresh Incidence p-value threshold for SNPs entering the mixture.
#' @param M Bootstrap samples for the slope SE (0 = no bootstrap; SE then 0).
#' @param seed Seed for the bootstrap.
#' @param n Sample size to record.
#' @return Data frame: snp, beta, se, z, p, n (adjusted effects), with the
#'   hunt() object as attribute "fit".
slope_hunter_correct <- function(snp, beta_x, se_x, beta_y, se_y,
                                 xp_thresh = 0.001, M = 100, seed = 777,
                                 n = NA) {
  if (!requireNamespace("SlopeHunter", quietly = TRUE))
    stop("Install SlopeHunter: remotes::install_github('Osmahmoud/SlopeHunter')")

  dat <- data.frame(SNP = snp,
                    BETA.incidence = beta_x, SE.incidence = se_x,
                    Pval.incidence = 2 * stats::pnorm(-abs(beta_x / se_x)),
                    BETA.prognosis = beta_y, SE.prognosis = se_y,
                    Pval.prognosis = 2 * stats::pnorm(-abs(beta_y / se_y)))

  # hunt() prints every bootstrap iteration; silence it
  utils::capture.output(
    fit <- SlopeHunter::hunt(dat, xp_thresh = xp_thresh,
                             Bootstrapping = M > 0, M = max(M, 1), seed = seed,
                             Plot = FALSE, show_adjustments = TRUE),
    file = nullfile())

  est <- fit$est[match(snp, fit$est$SNP), ]
  out <- data.frame(snp = snp, beta = est$ybeta_adj, se = est$yse_adj)
  out$z <- out$beta / out$se
  out$p <- est$yp_adj
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

    adj <- slope_hunter_correct(s$snp, s$beta_x, s$se_x, s$beta_y, s$se_y,
                                n = d$n_cases)
    fit <- attr(adj, "fit")

    # Agreement of the Hunted cluster with the true X-only SNPs
    hunted <- as.character(fit$Fit$SNP[fit$Fit$clusters == "Hunted"])
    truth  <- s$group[match(hunted, s$snp)]

    cat(sprintf(paste0("\nModule %s (%s): Slope-Hunter slope %.3f (SE %.3f); ",
                       "pi %.2f; %d SNPs hunted (%.0f%% truly X-only)\n"),
                sc, d$params$label, fit$b, fit$bse, fit$pi,
                length(hunted), 100 * mean(truth == "X_only")))
    print(rbind(
      compare_by_group(d$gwas_y, d$truth_y, d$groups, "naive (cases)"),
      compare_by_group(adj,      d$truth_y, d$groups, "Slope-Hunter")
    )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
          "type1_error")], row.names = FALSE)
  }
}
