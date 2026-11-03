# =============================================================================
# c_index_event.R
# Simulation module (c): index event (collider) bias
#
# DAG (Figure 2c): G -> X (incidence), G -> Y (progression), U -> X, U -> Y.
# Y is analysed only in cases (X = 1). Conditioning on X induces an
# association between X-affecting SNPs and Y, approximately -b * beta_X.
#
# Scenarios:
#   c1  independent: SNP effects on X and Y uncorrelated (DHO assumption holds)
#   c2  correlated: shared SNPs have correlated effects (rho = 0.6), plus an
#       incidence-only cluster (Slope-Hunter)
#   c3  sparse: few SNPs affect Y (MR-Horse horseshoe prior)
#
# Output (summary statistics only; methods for (c) are summary-level):
#   data/module_c1.rds, data/module_c2.rds, data/module_c3.rds
# Memory: ~1-2 GB while running (50,000 x 1,000 genotypes).
# =============================================================================

source("1_functions/simulate_genotypes.R")
source("1_functions/gwas_helpers.R")

base_c <- list(n = 50000, m = 1000, h2_x = 0.30, h2_y = 0.20,
               b_ux = 0.60, b_uy = 0.50, K_x = 0.30, n_neg_control = 20)

params_c <- list(
  c1 = modifyList(base_c, list(label = "independent", n_x_only = 50,
                               n_y_only = 50, n_both = 50, rho_both = 0)),
  c2 = modifyList(base_c, list(label = "correlated", n_x_only = 50,
                               n_y_only = 50, n_both = 50, rho_both = 0.6)),
  c3 = modifyList(base_c, list(label = "sparse", n_x_only = 100,
                               n_y_only = 10, n_both = 10, rho_both = 0))
)

#' Draw SNP effects on X (liability) and Y with optional correlation
draw_effects_c <- function(m, maf, n_x_only, n_y_only, n_both, rho_both,
                           h2_x, h2_y) {
  idx <- sample(m, n_x_only + n_y_only + n_both)
  ix  <- idx[seq_len(n_x_only)]
  iy  <- idx[n_x_only + seq_len(n_y_only)]
  ib  <- idx[n_x_only + n_y_only + seq_len(n_both)]
  
  zb <- matrix(stats::rnorm(2 * n_both), ncol = 2) %*%
    chol(matrix(c(1, rho_both, rho_both, 1), 2))
  zx <- zy <- numeric(m)
  zx[ix] <- stats::rnorm(n_x_only); zx[ib] <- zb[, 1]
  zy[iy] <- stats::rnorm(n_y_only); zy[ib] <- zb[, 2]
  
  v <- 2 * maf * (1 - maf)
  list(b_x = zx * sqrt(h2_x / (n_x_only + n_both) / v),
       b_y = zy * sqrt(h2_y / (n_y_only + n_both) / v),
       groups = list(X_only = sort(ix), Y_only = sort(iy), both = sort(ib),
                     null = sort(setdiff(seq_len(m), idx))))
}

#' Simulate one index event scenario and return summary statistics
simulate_c <- function(params, seed) {
  set.seed(seed)
  with(params, {
    G   <- simulate_genotypes(n, m)
    maf <- attr(G, "maf")
    eff <- draw_effects_c(m, maf, n_x_only, n_y_only, n_both, rho_both,
                          h2_x, h2_y)
    
    U    <- stats::rnorm(n)
    liab <- genetic_score(G, eff$b_x) + b_ux * U + stats::rnorm(n)
    X    <- as.integer(liab > stats::quantile(liab, 1 - K_x))
    Y    <- genetic_score(G, eff$b_y) + b_uy * U + stats::rnorm(n)
    case <- X == 1
    
    gwas_x     <- run_gwas(X, G, family = "binomial")   # incidence, full cohort
    gwas_y     <- run_gwas(Y[case], G[case, ])          # progression, cases only
    gwas_y_pop <- run_gwas(Y, G)                        # reference, unselected
    
    neg_control <- head(eff$groups$X_only, n_neg_control)
    group <- rep("null", m)
    for (g in names(eff$groups)) group[eff$groups[[g]]] <- g
    
    sumstats <- data.frame(
      snp = colnames(G), maf = maf, group = group,
      neg_control = seq_len(m) %in% neg_control,
      beta_x = gwas_x$beta, se_x = gwas_x$se, p_x = gwas_x$p,
      beta_y = gwas_y$beta, se_y = gwas_y$se, p_y = gwas_y$p,
      beta_y_pop = gwas_y_pop$beta, se_y_pop = gwas_y_pop$se,
      truth_x_liab = eff$b_x, truth_y = eff$b_y
    )
    
    # Oracle collider slope: beta_Y on beta_X across X-only SNPs (true b_y = 0)
    xo <- sumstats$group == "X_only"
    oracle_slope <- unname(stats::coef(
      stats::lm(beta_y ~ 0 + beta_x, data = sumstats[xo, ])))
    
    list(sumstats = sumstats, gwas_x = gwas_x, gwas_y = gwas_y,
         gwas_y_pop = gwas_y_pop, truth_y = eff$b_y, groups = eff$groups,
         neg_control = neg_control, oracle_slope = oracle_slope,
         n_cases = sum(case), params = params, seed = seed)
  })
}

# ---- Run all scenarios -------------------------------------------------------
dir.create("data", showWarnings = FALSE)
seeds <- c(c1 = 3001, c2 = 3002, c3 = 3003)

for (sc in names(params_c)) {
  d <- simulate_c(params_c[[sc]], seeds[[sc]])
  saveRDS(d, sprintf("data/module_%s.rds", sc))
  cat(sprintf("\nModule %s (%s) saved: %d cases of %d; oracle collider slope %.3f\n",
              sc, d$params$label, d$n_cases, d$params$n, d$oracle_slope))
  print(rbind(
    compare_by_group(d$gwas_y_pop, d$truth_y, d$groups, "population (reference)"),
    compare_by_group(d$gwas_y,     d$truth_y, d$groups, "naive (cases only)")
  )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
        "type1_error")], row.names = FALSE)
}
