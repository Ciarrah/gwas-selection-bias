# =============================================================================
# negative_controls.R
# Diagnostic: negative controls
#
# 1. Negative control outcomes: phenotypes with no genetic effects. Among
#    participants, SNPs linked to selection become associated with them
#    (module a). Compare with the same analysis in an unselected reference.
# 2. Negative control variants: SNPs with no effect on the outcome. In
#    case-only GWAS they acquire spurious associations through index-event
#    bias (module c; example added with that module).
#
# A null result does not rule out selection: power depends on how strongly the
# negative control shares causes with selection.
#
# References: Mitchell et al. (2023); Dahabreh and Kent (2011);
#             Devlin and Roeder (1999) for the genomic inflation factor.
# =============================================================================

source("1_functions/gwas_helpers.R")

#' GWAS of a negative control outcome in the sample and, optionally, a reference
#'
#' @param N Negative control phenotype (length n).
#' @param G Genotype matrix (n x m).
#' @param selected Logical or 0/1 vector: in the analysed sample.
#' @param reference Optional logical/0-1 vector defining an unselected reference
#'   sample (e.g. full cohort or an independent population sample).
#' @param covar Optional covariates.
#' @return List with `sample` and (if requested) `reference` run_gwas() results.
negative_control_outcome <- function(N, G, selected, reference = NULL,
                                     covar = NULL) {
  sel <- as.logical(selected)
  cv  <- function(idx) if (is.null(covar)) NULL else as.matrix(covar)[idx, , drop = FALSE]

  out <- list(sample = run_gwas(N[sel], G[sel, , drop = FALSE], covar = cv(sel)))
  if (!is.null(reference)) {
    ref <- as.logical(reference)
    out$reference <- run_gwas(N[ref], G[ref, , drop = FALSE], covar = cv(ref))
  }
  out
}

#' Summarise negative control associations
#'
#' @param res run_gwas() output for the negative control.
#' @param idx Optional subset of SNP rows (e.g. negative control variants).
#' @param alpha Nominal significance level.
#' @return One-row data frame: n_snps, lambda, prop_p_lt_alpha, n_gw_sig.
summarise_negative_control <- function(res, idx = NULL, alpha = 0.05) {
  if (!is.null(idx)) res <- res[idx, , drop = FALSE]
  p <- res$p[!is.na(res$p)]
  data.frame(
    n_snps          = length(p),
    lambda          = genomic_lambda(p),
    prop_p_lt_alpha = mean(p < alpha),   # expected: alpha
    n_gw_sig        = sum(p < 5e-8)      # expected: ~0
  )
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0) {
  d <- readRDS("data/module_a.rds")
  p <- d$pheno

  for (sc in c("S_a1", "S_a2")) {
    cat("\n====", sc, ": negative control outcome N ====\n")

    nc <- negative_control_outcome(p$N, d$G, selected = p[[sc]] == 1,
                                   reference = rep(TRUE, nrow(p)))

    cat("All SNPs:\n")
    print(rbind(
      cbind(sample = "participants", summarise_negative_control(nc$sample)),
      cbind(sample = "reference",    summarise_negative_control(nc$reference))
    ))

    cat("By SNP group (participants):\n")
    print(do.call(rbind, lapply(names(d$groups), function(g)
      cbind(snp_group = g, summarise_negative_control(nc$sample, d$groups[[g]])))))
  }
}
