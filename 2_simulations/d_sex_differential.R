# =============================================================================
# d_sex_differential.R
# Simulation module (d): sex-differential participation (Figure 2d)
#
# DAG: G -> Y, Sex -> Y; participation S depends on G (sex-specific effects),
# Sex, and Y (opposite direction in men and women), following Lee et al.
# (2022, Genome Biol; hanbin973/autosomalsexgwas).
#
# Pooling sexes induces associations between participation SNPs and Sex, and
# hence with Y. Sex-stratified GWAS + fixed-effects meta-analysis removes this
# when participation is rare (logistic selection ~ multiplicative), but not
# when participation is common.
#
# Scenarios: d_low (5% participation), d_high (50% participation),
# both with ~15,000 participants.
# Output: data/module_d_low.rds, data/module_d_high.rds (participants only)
# =============================================================================

source("1_functions/simulate_genotypes.R")
source("1_functions/gwas_helpers.R")

base_d <- list(n_part = 15000, m = 200, n_causal_y = 20, n_causal_p = 50,
               h2_y = 0.10,    # SNP heritability of Y
               b_sex_y = 0.50, # effect of Sex on Y
               h2_p = 0.30,    # participation liability explained by SNPs
               rg_p = 0,       # between-sex correlation of participation SNP effects
               g_y_sex = 0.50, # Y effect on participation: +g (men), -g (women)
               g_sex = 0.30)   # main effect of Sex on participation

params_d <- list(
  d_low  = modifyList(base_d, list(label = "low participation",  rate = 0.05)),
  d_high = modifyList(base_d, list(label = "high participation", rate = 0.50))
)

#' Logistic intercept giving a target participation rate
calibrate_logit <- function(lin, rate) {
  stats::uniroot(function(a) mean(stats::plogis(a + lin)) - rate,
                 interval = c(-30, 30))$root
}

#' Simulate one sex-differential participation scenario
simulate_d <- function(params, seed) {
  set.seed(seed)
  with(params, {
    n   <- round(n_part / rate)
    G   <- simulate_genotypes(n, m)
    maf <- attr(G, "maf")

    idx   <- sample.int(m)
    idx_y <- sort(idx[seq_len(n_causal_y)])
    idx_p <- sort(idx[n_causal_y + seq_len(n_causal_p)])
    idx_0 <- setdiff(seq_len(m), c(idx_y, idx_p))

    b_y <- make_snp_effects(m, n_causal_y, h2 = h2_y, maf = maf,
                            causal_idx = idx_y)

    # Sex-specific participation effects (correlation rg_p between sexes)
    z  <- matrix(stats::rnorm(2 * n_causal_p), ncol = 2) %*%
          chol(matrix(c(1, rg_p, rg_p, 1), 2))
    sc <- sqrt(h2_p / n_causal_p / (2 * maf[idx_p] * (1 - maf[idx_p])))
    b_p_f <- b_p_m <- numeric(m)
    b_p_f[idx_p] <- z[, 1] * sc
    b_p_m[idx_p] <- z[, 2] * sc

    Sex <- stats::rbinom(n, 1, 0.5)                       # 1 = male
    Y   <- genetic_score(G, b_y) + b_sex_y * Sex +
           stats::rnorm(n, sd = sqrt(1 - h2_y))

    gp  <- ifelse(Sex == 1, genetic_score(G, b_p_m), genetic_score(G, b_p_f))
    lin <- gp + g_sex * Sex + g_y_sex * (2 * Sex - 1) * Y
    a0  <- calibrate_logit(lin, rate)
    S   <- stats::rbinom(n, 1, stats::plogis(a0 + lin)) == 1

    list(pheno  = data.frame(id = which(S), Sex = Sex[S], Y = Y[S]),
         G      = G[S, ],
         maf    = maf,
         truth  = as.numeric(b_y),
         groups = list(Y_snps = idx_y, P_snps = idx_p, null = idx_0),
         n_pop  = n, rate_obs = mean(S), params = params, seed = seed)
  })
}

# ---- Run both scenarios ------------------------------------------------------
dir.create("data", showWarnings = FALSE)
seeds <- c(d_low = 4001, d_high = 4002)

for (sc in names(params_d)) {
  d <- simulate_d(params_d[[sc]], seeds[[sc]])
  saveRDS(d, sprintf("data/module_%s.rds", sc))
  p <- d$pheno

  cat(sprintf("\nModule %s (%s) saved: %d participants of %d (%.1f%%); %.0f%% male\n",
              sc, d$params$label, nrow(p), d$n_pop, 100 * d$rate_obs,
              100 * mean(p$Sex)))

  pooled <- run_gwas(p$Y, d$G)
  adjsex <- run_gwas(p$Y, d$G, covar = p[, "Sex", drop = FALSE])
  print(rbind(
    compare_by_group(pooled, d$truth, d$groups, "pooled"),
    compare_by_group(adjsex, d$truth, d$groups, "pooled, sex-adjusted")
  )[, c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
        "type1_error")], row.names = FALSE)
}
