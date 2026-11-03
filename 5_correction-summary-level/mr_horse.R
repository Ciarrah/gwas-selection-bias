# =============================================================================
# mr_horse.R
# Correction (summary-level): MR-Horse adapted for index event bias
#
# Bayesian model (mr_horse_model.jags) fitted with rjags. The collider slope
# theta is estimated allowing sparse direct effects on progression (horseshoe)
# that may be correlated with incidence effects (rho ~ 2*Beta(3,1) - 1).
# Adjusted effect = beta_Y - theta * beta_X, as for DHO and Slope-Hunter.
#
# Model code follows Grant and Burgess (aj-grant/mrhorse); prior on rho as
# adapted by Donovan et al. Requires JAGS (https://mcmc-jags.sourceforge.io).
#
# References: Grant and Burgess (2024) Am J Hum Genet; Donovan et al.
# =============================================================================

source("1_functions/gwas_helpers.R")

#' Estimate the collider slope with MR-Horse
#'
#' @param beta_x,se_x,beta_y,se_y Summary statistics.
#' @param use Logical vector of SNPs to fit the model on.
#' @param n_iter,n_burnin Iterations kept and discarded per chain.
#' @param n_chains Number of chains.
#' @param seed Base RNG seed (chain k uses seed + k).
#' @param model_file Path to the JAGS model.
#' @return List: b (posterior mean), se (posterior SD), ci, rhat, n_snps, samples.
mr_horse_slope <- function(beta_x, se_x, beta_y, se_y, use,
                           n_iter = 10000, n_burnin = 10000, n_chains = 3,
                           seed = 2025,
                           model_file = "3_methods/summary_level/mr_horse_model.jags") {
  if (!requireNamespace("rjags", quietly = TRUE))
    stop("Install JAGS and the rjags package (see 0_install/install.R).")

  dat   <- list(N = sum(use), bx = beta_x[use], sx = se_x[use],
                by = beta_y[use], sy = se_y[use])
  inits <- lapply(seq_len(n_chains), function(k)
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + k))

  jm <- rjags::jags.model(model_file, data = dat, inits = inits,
                          n.chains = n_chains, quiet = TRUE)
  stats::update(jm, n_burnin, progress.bar = "none")
  smp <- rjags::coda.samples(jm, "theta", n.iter = n_iter,
                             progress.bar = "none")
  th  <- as.matrix(smp)[, "theta"]

  list(b = mean(th), se = stats::sd(th),
       ci = stats::quantile(th, c(0.025, 0.975)),
       rhat = coda::gelman.diag(smp)$psrf[1, 1],
       n_snps = dat$N, samples = smp)
}

#' Apply the MR-Horse slope to all SNPs
#'
#' SE formula as in indexevent and SlopeHunter.
#' @return Data frame: snp, beta, se, z, p, n.
mr_horse_correct <- function(snp, beta_x, se_x, beta_y, se_y, slope, n = NA) {
  b <- slope$b; bse <- slope$se
  out <- data.frame(snp = snp, beta = beta_y - b * beta_x,
                    se = sqrt(se_y^2 + b^2 * se_x^2 + bse^2 * beta_x^2 +
                              bse^2 * se_x^2))
  out$z <- out$beta / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  out$n <- n
  out
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0) {
  for (sc in c("c1", "c2", "c3")) {
    f <- sprintf("data/module_%s.rds", sc)
    if (!file.exists(f)) next
    d <- readRDS(f)
    s <- d$sumstats

    slope <- mr_horse_slope(s$beta_x, s$se_x, s$beta_y, s$se_y,
                            use = s$p_x < 1e-3)
    adj   <- mr_horse_correct(s$snp, s$beta_x, s$se_x, s$beta_y, s$se_y,
                              slope, n = d$n_cases)

    cat(sprintf(paste0("\nModule %s (%s): MR-Horse slope %.3f (SD %.3f; ",
                       "95%% CrI %.3f to %.3f; Rhat %.3f; %d SNPs)\n"),
                sc, d$params$label, slope$b, slope$se, slope$ci[1],
                slope$ci[2], slope$rhat, slope$n_snps))
    print(rbind(
      compare_by_group(d$gwas_y, d$truth_y, d$groups, "naive (cases)"),
      compare_by_group(adj,      d$truth_y, d$groups, "MR-Horse")
    )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
          "type1_error")], row.names = FALSE)
  }
}
