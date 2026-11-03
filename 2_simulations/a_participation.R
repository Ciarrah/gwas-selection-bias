# =============================================================================
# a_participation.R
# Module (a): non-representative participation (Figure 2a)
#
#   Outcome:    G -> Y, G -> Z, Z -> Y, U -> Y
#   Selection:  a1: Z -> S, E -> S            (measured; missing at random)
#               a2: Z -> S, E -> S, U -> S    (+ unmeasured; missing not at random)
#
# Z: measured heritable predictor of selection (e.g. education)
# E: predicts selection but not Y (Heckman exclusion variable)
# U: unmeasured factor affecting Y and (in a2) selection
# N: negative control phenotype (no genetic effects; shares causes with selection)
#
# Selection is within-cohort nonresponse: Z, E and N are observed for everyone,
# Y only for responders (as in Mignogna et al., 2023).
#
# Truth: total SNP effects (naive GWAS, IPW) and direct effects given Z (Heckman).
# Output: data/module_a.rds
# =============================================================================

source("1_functions/simulate_genotypes.R")
source("1_functions/gwas_helpers.R")
source("1_functions/selection_helpers.R")

set.seed(2025)

# ---- Parameters -------------------------------------------------------------
params <- list(
  n          = 50000,  # cohort size
  m          = 1000,   # SNPs
  n_causal_y = 20,     # SNPs with direct effects on Y
  n_causal_z = 20,     # SNPs affecting Z (disjoint from Y SNPs)
  h2_y       = 0.10,   # variance in Y explained by direct SNP effects
  h2_z       = 0.20,   # SNP heritability of Z
  b_zy       = 0.30,   # effect of Z on Y
  b_uy       = 0.30,   # effect of U on Y
  g_z        = 0.60,   # effect of Z on selection (probit scale)
  g_e        = 0.60,   # effect of E on selection
  g_u        = 0.60,   # effect of U on selection (a2 only)
  rate       = 0.30    # target response rate
)

# ---- Genotypes and SNP effects ----------------------------------------------
G   <- simulate_genotypes(params$n, params$m)
maf <- attr(G, "maf")

idx   <- sample.int(params$m)
idx_y <- sort(idx[seq_len(params$n_causal_y)])
idx_z <- sort(idx[params$n_causal_y + seq_len(params$n_causal_z)])

b_y <- make_snp_effects(params$m, params$n_causal_y, h2 = params$h2_y,
                        maf = maf, causal_idx = idx_y)
b_z <- make_snp_effects(params$m, params$n_causal_z, h2 = params$h2_z,
                        maf = maf, causal_idx = idx_z)

# ---- Phenotypes -------------------------------------------------------------
eps_z <- rnorm(params$n, sd = sqrt(1 - params$h2_z))
Z     <- genetic_score(G, b_z) + eps_z                  # Var(Z) ~ 1
U     <- rnorm(params$n)
E     <- rnorm(params$n)

Y <- genetic_score(G, b_y) + params$b_zy * Z + params$b_uy * U +
     rnorm(params$n, sd = sqrt(1 - params$h2_y))

# Negative control: no genetic effects; shares the non-genetic part of Z and U
N <- 0.5 * eps_z + 0.5 * U + rnorm(params$n, sd = sqrt(0.5))

# True population SNP effects on Y
truth        <- as.numeric(b_y + params$b_zy * b_z)  # total effect (naive GWAS, IPW)
truth_direct <- as.numeric(b_y)                      # direct effect given Z (Heckman)

# ---- Selection (probit; intercept calibrated to the target response rate) --
lin_a1 <- params$g_z * Z + params$g_e * E
lin_a2 <- lin_a1 + params$g_u * U

S_a1 <- draw_selection(lin_a1, params$rate)
S_a2 <- draw_selection(lin_a2, params$rate)

# ---- Assemble and save ------------------------------------------------------
pheno <- data.frame(
  id   = seq_len(params$n),
  Y    = Y,
  Z    = Z,
  E    = E,
  N    = N,
  U    = U,                 # unobserved in practice; kept only for checking
  S_a1 = as.integer(S_a1),
  S_a2 = as.integer(S_a2)
)

module_a <- list(
  pheno        = pheno,
  G            = G,
  maf          = maf,
  truth        = truth,               # total effects on Y (naive, IPW)
  truth_direct = truth_direct,        # direct effects on Y given Z (Heckman)
  truth_N      = rep(0, params$m),    # true effects on the negative control
  params       = params,
  causal_y     = idx_y,
  causal_z     = idx_z
)

dir.create("data", showWarnings = FALSE)
saveRDS(module_a, "data/module_a.rds")

message(sprintf("Module (a) saved. Response rate: a1 = %.2f, a2 = %.2f",
                mean(S_a1), mean(S_a2)))
