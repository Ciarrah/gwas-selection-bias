# =============================================================================
# participation_gwas.R
# Diagnostic: GWAS of participation (selected vs not selected)
#
# Identifies SNPs associated with study inclusion. An excess of associated
# SNPs relative to the null expectation indicates genetic selection, which
# will distort SNP-trait associations among participants.
#
# Three versions:
#   1. Participation indicator available (e.g. register linkage): logistic
#      GWAS of S on each SNP in the full sampling frame.
#   2. No indicator: compare allele frequencies in the sample with a
#      population reference (here the full simulated cohort).
#   3. No indicator, sex-differential participation: GWAS of sex among
#      participants (autosomal SNPs cannot cause sex).
#
# Modules: (a), (d). References: Pirastu et al. (2021); Tyrrell et al. (2021);
#          Benonisdottir and Kong (2023); Lee et al. (2022).
# =============================================================================

source("1_functions/gwas_helpers.R")

#' Participation GWAS using a selection indicator
#'
#' @param S 0/1 participation indicator for everyone in the sampling frame.
#' @param G Genotype matrix for the same individuals.
#' @param covar Optional covariates (e.g. sex, age, principal components).
#' @return run_gwas() output (log-odds per allele).
participation_gwas <- function(S, G, covar = NULL) {
  run_gwas(S, G, covar = covar, family = "binomial")
}

#' Allele-frequency comparison: sample vs population reference
#'
#' Two-sample test of allele frequency per SNP. Use when no participation
#' indicator exists but reference allele frequencies are available.
#'
#' @param G_sample Genotypes of participants.
#' @param G_ref Genotypes (or a reference panel) for the target population.
#' @param nested TRUE if the participants are part of the reference sample
#'   (e.g. the full sampling frame); FALSE for an independent reference panel.
#' @return Data frame: snp, af_sample, af_ref, diff, z, p.
allele_freq_test <- function(G_sample, G_ref, nested = FALSE) {
  af_s <- colMeans(G_sample) / 2
  af_r <- colMeans(G_ref) / 2
  n_s  <- 2 * nrow(G_sample)
  n_r  <- 2 * nrow(G_ref)

  if (nested) {
    stopifnot(n_s < n_r)
    p0 <- af_r                                         # reference contains the sample
    se <- sqrt(p0 * (1 - p0) * (1 / n_s - 1 / n_r))
  } else {
    p0 <- (af_s * n_s + af_r * n_r) / (n_s + n_r)      # pooled frequency
    se <- sqrt(p0 * (1 - p0) * (1 / n_s + 1 / n_r))
  }
  z <- (af_s - af_r) / se
  data.frame(snp = colnames(G_sample), af_sample = af_s, af_ref = af_r,
             diff = af_s - af_r, z = z, p = 2 * stats::pnorm(-abs(z)),
             row.names = NULL)
}

#' Participation GWAS without a participation indicator: GWAS of sex
#'
#' Autosomal SNPs cannot cause biological sex, so SNP-sex associations among
#' participants indicate sex-differential participation (Pirastu et al. 2021;
#' Lee et al. 2022). Summarise with summarise_participation().
#'
#' @param sex 0/1 sex (1 = male) among participants.
#' @param G Genotype matrix for participants.
#' @param covar Optional covariates (e.g. principal components).
#' @return run_gwas() output (log-odds of male per allele).
sex_gwas <- function(sex, G, covar = NULL) {
  run_gwas(sex, G, covar = covar, family = "binomial")
}

#' Summarise a participation scan
#'
#' Compares the number of associated SNPs with the number expected by chance.
#'
#' @param res Output of participation_gwas(), allele_freq_test() or sex_gwas().
#' @param thresholds P-value thresholds to report.
#' @return Data frame: threshold, observed, expected, ratio; attribute "lambda".
summarise_participation <- function(res, thresholds = c(0.05, 1e-3, 5e-8)) {
  p   <- res$p[!is.na(res$p)]
  out <- data.frame(
    threshold = thresholds,
    observed  = sapply(thresholds, function(t) sum(p < t)),
    expected  = thresholds * length(p)
  )
  out$ratio <- out$observed / out$expected
  attr(out, "lambda") <- genomic_lambda(p)
  out
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0) {

  # Versions 1 and 2 (module a)
  if (file.exists("data/module_a.rds")) {
    d <- readRDS("data/module_a.rds")
    p <- d$pheno

    for (sc in c("S_a1", "S_a2")) {
      cat("\n====", sc, "====\n")

      # Version 1: participation indicator
      pg <- participation_gwas(p[[sc]], d$G)
      s1 <- summarise_participation(pg)
      cat("Participation GWAS (lambda =", round(attr(s1, "lambda"), 2), ")\n")
      print(s1)

      # Which SNP groups carry the signal?
      by_group <- sapply(d$groups, function(i) mean(pg$p[i] < 0.05))
      cat("Proportion of SNPs with p < 0.05 by group:\n")
      print(round(by_group, 2))

      # Version 2: allele-frequency comparison with the full cohort as reference
      af <- allele_freq_test(d$G[p[[sc]] == 1, ], d$G, nested = TRUE)
      s2 <- summarise_participation(af)
      cat("Allele-frequency comparison (lambda =", round(attr(s2, "lambda"), 2), ")\n")
      print(s2)
    }
  }

  # Version 3: GWAS of sex (module d)
  for (sc in c("d_low", "d_high")) {
    f <- sprintf("data/module_%s.rds", sc)
    if (!file.exists(f)) next
    d  <- readRDS(f)
    sg <- sex_gwas(d$pheno$Sex, d$G)
    s3 <- summarise_participation(sg)
    cat(sprintf("\n==== %s: GWAS of sex (lambda = %.2f) ====\n",
                sc, attr(s3, "lambda")))
    print(s3)
    cat("Proportion of SNPs with p < 0.05 by group:\n")
    print(round(sapply(d$groups, function(i) mean(sg$p[i] < 0.05)), 2))
  }
}
