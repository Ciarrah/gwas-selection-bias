# =============================================================================
# spreg_prepare.R
# Correction (individual-level): SPREG for case-control secondary phenotypes
#
# SPREG (Lin and Zeng 2009) is external C software, not reimplemented here:
#   https://dlin.web.unc.edu/software/spreg-2/
# It jointly models disease status and the secondary phenotype, using the
# population disease prevalence to correct for case-control ascertainment.
# The executable is built for 64-bit x86 Linux; on macOS or Windows, run it
# on a Linux server or in a Linux container (e.g. Docker).
#
# This script:
#   1. writes module (b) data in SPREG's input format;
#   2. reads SPREG's output back and compares it with the truth.
#
# SPREG input (one text matrix, subjects as COLUMNS):
#   row 1          case-control status (1 = case, 0 = control)
#   row 2          secondary phenotype
#   rows 3..d+2    environmental covariates (optional)
#   remaining rows genotypes, one row per SNP (0/1/2)
#   missing values coded NA
# Command:  spreg infile outfile rate nenv
#   rate = disease prevalence (any value < 0.01 if the disease is rare)
#   nenv = number of covariate rows (0 here)
# Output columns: Gene_number, Estimate, Std_Error, Z_stat, p_value
#
# Module: (b). Reference: Lin and Zeng (2009) Genet Epidemiol 33:256-265.
# =============================================================================

source("1_functions/gwas_helpers.R")

#' Write module (b) data in SPREG input format
#'
#' @param d Module (b) object (readRDS("data/module_b.rds")).
#' @param outdir Output folder.
#' @param pheno Secondary phenotype column to export.
#' @param covar Optional data frame of covariates (one column per covariate).
#' @return Invisibly, a list with the input path, SNP order and the command
#'   to run.
export_spreg <- function(d, outdir = "data/spreg", pheno = "Y1",
                         covar = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  p <- d$pheno

  mat <- rbind(
    p$D,
    p[[pheno]],
    if (!is.null(covar)) t(as.matrix(covar)),
    t(d$G)                                   # SNPs as rows, subjects as columns
  )

  infile <- file.path(outdir, "spreg_input.dat")
  utils::write.table(mat, infile, quote = FALSE, row.names = FALSE,
                     col.names = FALSE, na = "NA")

  # SNP order, so Gene_number in the output can be mapped back to SNP IDs
  utils::write.table(data.frame(Gene_number = seq_len(ncol(d$G)),
                                snp = colnames(d$G)),
                     file.path(outdir, "snp_order.txt"),
                     quote = FALSE, row.names = FALSE, sep = "\t")

  nenv <- if (is.null(covar)) 0 else ncol(covar)
  cmd  <- sprintf("spreg spreg_input.dat spreg_output.out %s %d", d$K, nenv)
  message("SPREG input written to ", infile, "\nRun in ", outdir, ":  ", cmd)
  invisible(list(infile = infile, snp = colnames(d$G), command = cmd))
}

#' Read SPREG output into the repository's GWAS format
#'
#' @param file SPREG output file.
#' @param snp SNP identifiers in input order (Gene_number 1, 2, ...).
#' @param n Sample size to record.
#' @return Data frame: snp, beta, se, z, p, n.
read_spreg <- function(file, snp, n = NA) {
  r <- utils::read.table(file, header = TRUE, stringsAsFactors = FALSE)
  names(r)[1:5] <- c("Gene_number", "Estimate", "Std_Error", "Z_stat", "p_value")
  r <- r[order(r$Gene_number), ]
  stopifnot(nrow(r) == length(snp))

  data.frame(snp  = snp[r$Gene_number],
             beta = as.numeric(r$Estimate),
             se   = as.numeric(r$Std_Error),
             z    = as.numeric(r$Z_stat),
             p    = as.numeric(r$p_value),
             n    = n)
}

# ---- Example (runs only when this file is run directly) ---------------------
if (sys.nframe() == 0 && file.exists("data/module_b.rds")) {
  d  <- readRDS("data/module_b.rds")
  ex <- export_spreg(d)

  # After running SPREG in data/spreg/ (see the command printed above):
  res_file <- "data/spreg/spreg_output.out"
  if (file.exists(res_file)) {
    spreg <- read_spreg(res_file, snp = ex$snp, n = nrow(d$pheno))
    print(compare_by_group(spreg, d$truth, d$groups, "SPREG (Y1)")[
      , c("analysis", "snp_group", "mean_abs_bias", "rmse", "coverage",
          "type1_error")], row.names = FALSE)
  } else {
    message("SPREG output not found. Run SPREG, then re-run this script.")
  }
}
