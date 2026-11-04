# =============================================================================
# 
# Installs all R packages used in this repository, check that JAGS
# is available (required for the MR-Horse adaptation).
#
# Run once from the repository root:  source("0_install/Rpackages.R")
# =============================================================================

# ---- CRAN packages ----------------------------------------------------------
cran_pkgs <- c(
  # Simulation, data handling, plotting
  "data.table", "dplyr", "ggplot2",
  # Diagnostics and corrections
  "ipw",              # inverse probability weighting
  "survey",           # weighted regression with robust standard errors
  "geepack",          # weighted GEE for the IPW scaled marginal model
  "sampleSelection",  # Heckman correction
  "rjags", "coda",    # MR-Horse adaptation (requires JAGS, see below)
  # Vignettes, reproducibility, GitHub installs
  "rmarkdown", "knitr", "renv", "remotes",
  # Figures and Shiny app
  "DiagrammeR", "DiagrammeRsvg", "rsvg", "shiny", "shinythemes", "readxl"
)

missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

# ---- GitHub packages --------------------------------------------------------
github_pkgs <- c(
  SlopeHunter = "Osmahmoud/SlopeHunter",     # Slope-Hunter (Mahmoud et al., 2022)
  indexevent  = "DudbridgeLab/indexevent"    # DHO reference implementation (validation only), DHO correction
)

for (pkg in names(github_pkgs)) {
  if (!pkg %in% rownames(installed.packages())) {
    message("Installing ", pkg, " from GitHub (", github_pkgs[[pkg]], ")")
    remotes::install_github(github_pkgs[[pkg]], upgrade = "never")
  }
}

# ---- JAGS check -------------------------------------------------------------
# rjags only loads if the JAGS system library is installed.
#   macOS:   brew install jags   (or installer from https://sourceforge.net/projects/mcmc-jags/)
#   Windows: installer from https://sourceforge.net/projects/mcmc-jags/
#   Linux:   sudo apt install jags
jags_ok <- tryCatch({
  suppressMessages(library(rjags))
  TRUE
}, error = function(e) FALSE)

if (jags_ok) {
  message("JAGS found: MR-Horse adaptation scripts will run.")
} else {
  warning("JAGS not found. Install JAGS, then rerun install.R. ",
          "All other methods will still run.")
}

# ---- SPREG note -------------------------------------------------------------
message("Note: SPREG is external software and is not installed by this script. ",
        "Download it from https://dlin.web.unc.edu/software/spreg-2/ ; ",
        "4_correction-individual-level/spreg_prepare.R. exports input files for it.")

# ---- Record package versions ------------------------------------------------
# After the first successful run, uncomment to create renv.lock:
# renv::init(bare = TRUE); renv::snapshot()

message("Setup complete.")
