# =============================================================================
# b_case_control.R
# Module (b): case-control ascertainment with secondary phenotypes (Figure 2b)
#
#   G -> Y (secondary phenotype), G -> D, Y -> D; sampling on disease D.
#
# Disease D follows a liability-threshold model with population prevalence K.
# Cases are over-sampled relative to the population, so analysing Y in the
# case-control sample conditions on D and distorts SNP-Y associations:
# SNPs affecting D only acquire spurious associations with Y, and SNPs
# affecting Y are biased.
#
# Two secondary phenotypes (Y1, Y2) share the same SNP effects on the
# standardised scale (the assumption of the IPW scaled marginal model);
# only Y1 affects disease.
#
# Truth: population SNP effects on Y1 and Y2 (zero for disease-only SNPs).
# simulate_b() can be called with different seeds for replicate analyses.
# Output: data/module_b.rds (case-control sample only)
# =============================================================================

source("1_functions/simulate_genotypes.R")
source("1_functions/gwas_helpers.R")

params_b <- list(
  n_pop      = 200000,  # source population
  m          = 200,     # SNPs (20 Y-acting, 20 D-acting, 160 null)
  n_causal_y = 20,      # SNPs affecting Y1 and Y2 (same standardised effects)
  n_causal_d = 20,      # SNPs affecting disease liability only
  h2_y       = 0.10,    # SNP heritability of Y1 and Y2
  h2_d       = 0.20,    # liability variance explained by D-acting SNPs
  g_yd       = 0.50,    # effect of Y1 on disease liability
  r_y        = 0.50,    # residual correlation between Y1 and Y2
  K          = 0.05,    # population prevalence of D
  n_cases    = 5000,
  n_controls = 5000
)

simulate_b <- function(params = params_b, seed = 2025) {
  set.seed(seed)
  n <- params$n_pop

  # ---- Genotypes and SNP effects --------------------------------------------
  G   <- simulate_genotypes(n, params$m)
  maf <- attr(G, "maf")

  idx   <- sample.int(params$m)
  idx_y <- sort(idx[seq_len(params$n_causal_y)])
  idx_d <- sort(idx[params$n_causal_y + seq_len(params$n_causal_d)])
  idx_0 <- setdiff(seq_len(params$m), c(idx_y, idx_d))

  b_y <- make_snp_effects(params$m, params$n_causal_y, h2 = params$h2_y,
                          maf = maf, causal_idx = idx_y)
  b_d <- make_snp_effects(params$m, params$n_causal_d, h2 = params$h2_d,
                          maf = maf, causal_idx = idx_d)

  # ---- Secondary phenotypes (variance ~1, correlated residuals) ------------
  e1 <- rnorm(n)
  e2 <- params$r_y * e1 + sqrt(1 - params$r_y^2) * rnorm(n)
  gy <- genetic_score(G, b_y)
  Y1 <- gy + sqrt(1 - params$h2_y) * e1
  Y2 <- gy + sqrt(1 - params$h2_y) * e2

  # ---- Disease: liability threshold -----------------------------------------
  L <- genetic_score(G, b_d) + params$g_yd * Y1 + rnorm(n)
  D <- as.integer(L > stats::quantile(L, 1 - params$K))

  # ---- Case-control sampling --------------------------------------------------
  cases    <- sample(which(D == 1), params$n_cases)
  controls <- sample(which(D == 0), params$n_controls)
  s        <- c(cases, controls)

  # Inverse sampling-probability weights (population / sample share)
  share_case <- params$n_cases / length(s)
  w <- ifelse(D[s] == 1, params$K / share_case,
                         (1 - params$K) / (1 - share_case))

  list(
    pheno = data.frame(
      id = s, D = D[s], Y1 = Y1[s], Y2 = Y2[s], w_sampling = w
    ),
    G       = G[s, ],
    maf     = maf,
    truth   = as.numeric(b_y),   # population effects on Y1 and Y2
    groups  = list(Y_snps = idx_y, D_snps = idx_d, null = idx_0),
    K       = params$K,
    params  = params,
    seed    = seed
  )
}

# ---- Generate and save one replicate ----------------------------------------
module_b <- simulate_b()
dir.create("data", showWarnings = FALSE)
saveRDS(module_b, "data/module_b.rds")

message(sprintf("Module (b) saved: %d cases, %d controls (population prevalence %.2f).",
                sum(module_b$pheno$D == 1), sum(module_b$pheno$D == 0), module_b$K))
