# =============================================================================
# simulate_genotypes.R
# Shared genotype and effect-size generators used by every simulation module.
# Contains functions only; no code runs when the file is sourced.
# =============================================================================

#' Simulate unrelated individuals' genotypes
#'
#' SNPs are independent (linkage equilibrium) and in Hardy-Weinberg
#' equilibrium; genotypes are coded 0/1/2 copies of the effect allele.
#'
#' @param n Number of individuals.
#' @param m Number of SNPs.
#' @param maf Optional vector of minor allele frequencies (length m).
#' @param maf_range Range for uniform MAF draws if `maf` is not supplied.
#' @return Integer matrix (n x m) with attribute "maf".
simulate_genotypes <- function(n, m, maf = NULL, maf_range = c(0.05, 0.5)) {
  if (is.null(maf)) maf <- runif(m, maf_range[1], maf_range[2])
  stopifnot(length(maf) == m)
  
  # Matrix fills column-wise, so rep(maf, each = n) gives SNP j its own MAF
  G <- matrix(rbinom(n * m, 2, rep(maf, each = n)), nrow = n, ncol = m)
  colnames(G) <- sprintf("SNP%04d", seq_len(m))
  attr(G, "maf") <- maf
  G
}

#' Simulate parent-offspring pairs (one parent genotyped)
#'
#' Used for proxy outcomes (module e). Each offspring inherits one allele per
#' SNP from the genotyped parent (chosen at random) and one from an untyped
#' parent drawn from the population.
#'
#' @inheritParams simulate_genotypes
#' @return List with G (offspring), G_P (parent) and maf.
simulate_parent_offspring <- function(n, m, maf = NULL, maf_range = c(0.05, 0.5)) {
  if (is.null(maf)) maf <- runif(m, maf_range[1], maf_range[2])
  stopifnot(length(maf) == m)
  p <- rep(maf, each = n)
  
  # Genotyped parent's two haplotypes
  h1 <- matrix(rbinom(n * m, 1, p), nrow = n, ncol = m)
  h2 <- matrix(rbinom(n * m, 1, p), nrow = n, ncol = m)
  G_P <- h1 + h2
  
  # Transmitted allele from the genotyped parent (Mendelian segregation)
  transmit    <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  from_parent <- ifelse(transmit == 1, h1, h2)
  
  # Allele from the untyped parent
  from_other <- matrix(rbinom(n * m, 1, p), nrow = n, ncol = m)
  
  G <- from_parent + from_other
  
  colnames(G) <- colnames(G_P) <- sprintf("SNP%04d", seq_len(m))
  list(G = G, G_P = G_P, maf = maf)
}

#' Generate SNP effect sizes
#'
#' Either a fixed per-allele effect for every causal SNP, or effects scaled so
#' the causal SNPs jointly explain a target heritability (h2), each
#' contributing equally. Signs are random. Non-causal SNPs have effect 0.
#'
#' @param m Number of SNPs.
#' @param n_causal Number of causal SNPs.
#' @param beta Fixed per-allele effect (used if `h2` is NULL).
#' @param h2 Target variance explained by all causal SNPs (phenotype variance 1).
#' @param maf MAF vector, required if `h2` is supplied.
#' @param causal_idx Optional indices of causal SNPs; random if NULL.
#' @return Numeric vector of length m with attribute "causal_idx".
make_snp_effects <- function(m, n_causal, beta = NULL, h2 = NULL,
                             maf = NULL, causal_idx = NULL) {
  if (is.null(causal_idx)) causal_idx <- sort(sample.int(m, n_causal))
  stopifnot(length(causal_idx) == n_causal)
  
  b <- numeric(m)
  if (!is.null(h2)) {
    stopifnot(!is.null(maf))
    p <- maf[causal_idx]
    # Var(beta_j * G_j) = 2 p (1 - p) beta_j^2 = h2 / n_causal
    b[causal_idx] <- sqrt(h2 / n_causal / (2 * p * (1 - p)))
  } else {
    stopifnot(!is.null(beta))
    b[causal_idx] <- beta
  }
  b[causal_idx] <- b[causal_idx] * sample(c(-1, 1), n_causal, replace = TRUE)
  
  attr(b, "causal_idx") <- causal_idx
  b
}

#' Genetic score (linear predictor) for a set of SNP effects
#'
#' @param G Genotype matrix (n x m).
#' @param b Effect vector of length m.
#' @return Numeric vector of length n.
genetic_score <- function(G, b) {
  stopifnot(ncol(G) == length(b))
  as.numeric(G %*% b)
}
