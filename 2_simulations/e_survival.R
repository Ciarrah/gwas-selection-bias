# =============================================================================
# e_survival.R
# Simulation module (e): survival to recruitment (Figure 2e)
#
# DAG: G_P -> G (transmission), G_P -> Y_P; G -> Y, G -> R, U -> Y, U -> R.
# Only people alive at recruitment (R = 1) are studied. Conditioning on R
# induces an association between survival SNPs and U, and hence with Y.
# Parental outcomes Y_P (recorded regardless of parental survival) depend on
# parental genotype and the parent's own U_P, so they are not affected by the
# participant's survival and can serve as a proxy outcome.
#
# Output: data/module_e.rds (full cohort, with survival indicator R)
# =============================================================================

source("1_functions/simulate_genotypes.R")
source("1_functions/gwas_helpers.R")

params_e <- list(
  n        = 50000, m = 200,
  n_y_only = 20,    # SNPs affecting Y only
  n_r_only = 20,    # SNPs affecting survival only
  n_both   = 20,    # SNPs affecting Y and survival
  h2_y     = 0.10,  # variance in Y explained by SNPs
  h2_r     = 0.20,  # survival liability explained by SNPs
  b_uy     = 0.50,  # effect of U on Y
  g_u      = 1.00,  # effect of U on survival (logit scale)
  surv     = 0.50   # proportion alive at recruitment
)

simulate_e <- function(params = params_e, seed = 5001) {
  set.seed(seed)
  with(params, {
    po  <- simulate_parent_offspring(n, m)
    G   <- po$G
    G_P <- po$G_P
    maf <- po$maf

    idx <- sample.int(m)
    iy  <- sort(idx[seq_len(n_y_only)])
    ir  <- sort(idx[n_y_only + seq_len(n_r_only)])
    ib  <- sort(idx[n_y_only + n_r_only + seq_len(n_both)])
    i0  <- setdiff(seq_len(m), c(iy, ir, ib))

    b_y <- make_snp_effects(m, n_y_only + n_both, h2 = h2_y, maf = maf,
                            causal_idx = sort(c(iy, ib)))
    b_r <- make_snp_effects(m, n_r_only + n_both, h2 = h2_r, maf = maf,
                            causal_idx = sort(c(ir, ib)))

    U   <- stats::rnorm(n)
    U_P <- stats::rnorm(n)
    Y   <- genetic_score(G,   b_y) + b_uy * U   + stats::rnorm(n, sd = sqrt(1 - h2_y))
    Y_P <- genetic_score(G_P, b_y) + b_uy * U_P + stats::rnorm(n, sd = sqrt(1 - h2_y))

    # Survival to recruitment (logistic; intercept calibrated to target rate)
    lin <- genetic_score(G, b_r) + g_u * U
    a0  <- stats::uniroot(function(a) mean(stats::plogis(a + lin)) - surv,
                          interval = c(-30, 30))$root
    R   <- stats::rbinom(n, 1, stats::plogis(a0 + lin))

    list(pheno  = data.frame(id = seq_len(n), Y = Y, Y_P = Y_P, R = R,
                             U = U),            # U unobserved in practice
         G      = G,
         maf    = maf,
         truth  = as.numeric(b_y),              # effects on Y (participant)
         truth_r = as.numeric(b_r),
         groups = list(Y_only = iy, R_only = ir, both = ib, null = i0),
         params = params, seed = seed)
  })
}

# ---- Generate and save -------------------------------------------------------
module_e <- simulate_e()
dir.create("data", showWarnings = FALSE)
saveRDS(module_e, "data/module_e.rds")

d <- module_e
p <- d$pheno
alive <- p$R == 1
cat(sprintf("Module (e) saved: %d of %d alive at recruitment (%.0f%%)\n",
            sum(alive), nrow(p), 100 * mean(alive)))

full <- run_gwas(p$Y, d$G)
surv <- run_gwas(p$Y[alive], d$G[alive, ])
print(rbind(
  compare_by_group(full, d$truth, d$groups, "full cohort (reference)"),
  compare_by_group(surv, d$truth, d$groups, "survivors (naive)")
)[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
      "type1_error")], row.names = FALSE)
