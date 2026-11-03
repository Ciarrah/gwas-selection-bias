# =============================================================================
# a_participation.R
# Module (a): non-representative participation (Figure 2a)
#
#   Outcome:    G -> Y, G -> Z, Z -> Y, U -> Y
#   Selection:  a1: Z -> S, E -> S            (measured; missing at random)
#               a2: Z -> S, E -> S, Y -> S    (+ selection on the outcome; MNAR)
#
# Z: measured heritable predictor of selection (e.g. education)
# E: predicts selection but not Y (Heckman exclusion variable)
# U: unmeasured factor affecting Y
# N: negative control phenotype (no genetic effects; shares causes with Y and Z)
#
# Selection is within-cohort nonresponse: Z, E and N are observed for everyone,
# Y only for responders (as in Mignogna et al., 2023).
#
# Truth: total SNP effects (naive GWAS, IPW) and direct effects given Z (Heckman).
# simulate_a() can be called with different seeds for replicate analyses.
# Output: data/module_a.rds (one replicate)
# =============================================================================

source("1_functions/simulate_genotypes.R")
source("1_functions/gwas_helpers.R")
source("1_functions/selection_helpers.R")

params_a <- list(
  n          = 50000,  # cohort size
  m          = 200,    # SNPs (20 Y-acting, 20 Z-acting, 160 null)
  n_causal_y = 20,     # SNPs with direct effects on Y
  n_causal_z = 20,     # SNPs affecting Z (disjoint from Y SNPs)
  h2_y       = 0.10,   # variance in Y explained by direct SNP effects
  h2_z       = 0.30,   # SNP heritability of Z
  b_zy       = 0.50,   # effect of Z on Y
  b_uy       = 0.50,   # effect of U on Y
  g_z        = 0.80,   # effect of Z on selection (probit scale)
  g_e        = 0.60,   # effect of E on selection
  g_y        = 0.60,   # effect of Y on selection (a2 only)
  rate       = 0.30    # target response rate
)

simulate_a <- function(params = params_a, seed = 97) {
  set.seed(seed)

  # ---- Genotypes and SNP effects --------------------------------------------
  G   <- simulate_genotypes(params$n, params$m)
  maf <- attr(G, "maf")

  idx   <- sample.int(params$m)
  idx_y <- sort(idx[seq_len(params$n_causal_y)])
  idx_z <- sort(idx[params$n_causal_y + seq_len(params$n_causal_z)])
  idx_0 <- setdiff(seq_len(params$m), c(idx_y, idx_z))

  b_y <- make_snp_effects(params$m, params$n_causal_y, h2 = params$h2_y,
                          maf = maf, causal_idx = idx_y)
  b_z <- make_snp_effects(params$m, params$n_causal_z, h2 = params$h2_z,
                          maf = maf, causal_idx = idx_z)

  # ---- Phenotypes -----------------------------------------------------------
  eps_z <- rnorm(params$n, sd = sqrt(1 - params$h2_z))
  Z     <- genetic_score(G, b_z) + eps_z
  U     <- rnorm(params$n)
  E     <- rnorm(params$n)

  Y <- genetic_score(G, b_y) + params$b_zy * Z + params$b_uy * U +
       rnorm(params$n, sd = sqrt(1 - params$h2_y))

  # Negative control: no genetic effects; shares non-genetic causes with Z and Y
  N <- 0.5 * eps_z + 0.5 * U + rnorm(params$n, sd = sqrt(0.5))

  # ---- Selection ------------------------------------------------------------
  lin_a1 <- params$g_z * Z + params$g_e * E
  lin_a2 <- lin_a1 + params$g_y * Y

  S_a1 <- draw_selection(lin_a1, params$rate)
  S_a2 <- draw_selection(lin_a2, params$rate)

  # ---- Output ---------------------------------------------------------------
  list(
    pheno = data.frame(
      id = seq_len(params$n), Y = Y, Z = Z, E = E, N = N,
      U = U,                                  # unobserved in practice
      S_a1 = as.integer(S_a1), S_a2 = as.integer(S_a2)
    ),
    G            = G,
    maf          = maf,
    truth        = as.numeric(b_y + params$b_zy * b_z),  # total (naive, IPW)
    truth_direct = as.numeric(b_y),                      # direct given Z (Heckman)
    truth_N      = rep(0, params$m),
    groups       = list(Y_snps = idx_y, Z_snps = idx_z, null = idx_0),
    params       = params,
    seed         = seed
  )
}

# ---- Generate and save one replicate ----------------------------------------
module_a <- simulate_a()
dir.create("data", showWarnings = FALSE)
saveRDS(module_a, "data/module_a.rds")

message(sprintf("Module (a) saved. Response rate: a1 = %.2f, a2 = %.2f",
                mean(module_a$pheno$S_a1), mean(module_a$pheno$S_a2)))
