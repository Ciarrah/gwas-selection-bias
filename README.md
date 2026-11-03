# gwas-selection-bias
Companion repository to “Navigating Selection Bias in Genome-Wide Association Studies”. R implementations of diagnostics (negative controls, proxy outcomes, participation GWAS) and corrections (IPW, Heckman, SPREG, sex-stratified meta-analysis, DHO, Slope-Hunter, MR-Horse) applied to simulated data.
Each simulation module corresponds to a panel of Figure 2 in the paper.

## Methods

| Method | Purpose | Data | Module | Script | Software |
|---|---|---|---|---|---|
| Participation GWAS (incl. GWAS of sex) | Diagnostic | Individual | a, d | `3_diagnostics/participation_gwas.R` | — |
| Negative controls (outcomes, variants) | Diagnostic | Individual / summary | a, c | `3_diagnostics/negative_controls.R` | — |
| Proxy outcomes | Diagnostic | Individual | e | `3_diagnostics/proxy_outcomes.R` | — |
| Inverse probability weighting (incl. scaled marginal model) | Correction | Individual | a, b | `4_correction-individual-level/ipw.R` | geepack |
| Heckman correction | Correction | Individual | a | `4_correction-individual-level/heckman.R` | sampleSelection |
| SPREG | Correction | Individual | b | `4_correction-individual-level/spreg_prepare.R` | [SPREG](https://dlin.web.unc.edu/software/spreg-2/) (external) |
| Sex-stratified GWAS meta-analysis | Correction | Summary | d | `5_correction-summary-level/sex_stratified.R` | — |
| DHO | Correction | Summary | c | `5_correction-summary-level/dho.R` | indexevent |
| Slope-Hunter | Correction | Summary | c | `5_correction-summary-level/slope_hunter.R` | SlopeHunter |
| MR-Horse adaptation | Correction | Summary | c | `5_correction-summary-level/mr_horse.R` | rjags (JAGS 4.x) |

## Simulation modules

| Module | Selection structure | Script |
|---|---|---|
| a | Non-representative participation (a1: on measured predictors; a2: also on the outcome) | `2_simulations/a_participation.R` |
| b | Case-control ascertainment, secondary phenotypes | `2_simulations/b_case_control.R` |
| c | Index event: progression among cases (c1 independent, c2 correlated pleiotropy, c3 sparse) | `2_simulations/c_index_event.R` |
| d | Sex-differential participation (low and high participation) | `2_simulations/d_sex_differential.R` |
| e | Survival to recruitment, with parental outcomes | `2_simulations/e_survival.R` |

## Getting started

Run everything from the repository root.

```r
# 1. Install packages (once)
source("0_install/install.R")

# 2. Simulate a module (writes data/module_*.rds)
source("2_simulations/a_participation.R")

# 3. Run a method's worked example
#    (examples run only when a script is run directly, not when sourced)
system("Rscript 3_diagnostics/participation_gwas.R")
```

Method scripts define functions that can be sourced and applied to your own
data; their example blocks run only when the script is run directly
(`Rscript <script>`).

Simulations use fixed seeds. Module (e) simulates 200,000 individuals and needs
about 1 GB of memory.

## Software notes

**JAGS (MR-Horse).** Install JAGS 4.x from
https://sourceforge.net/projects/mcmc-jags/ before `rjags`. The current CRAN
`rjags` does not work with JAGS 5.x, which Homebrew installs on macOS.

**SPREG (external).** SPREG is C software for 64-bit Linux and is not run
inside this repository.

1. Download it from https://dlin.web.unc.edu/software/spreg-2/
2. Run `4_correction-individual-level/spreg_prepare.R`, which writes
   `data/spreg/spreg_input.dat`.
3. In `data/spreg/`, run `spreg spreg_input.dat spreg_output.out 0.05 0`.
4. Re-run `spreg_prepare.R` to compare SPREG's estimates with the truth.

The IPW and scaled marginal model corrections for the same scenario run
entirely in R (`ipw.R`).

## Repository structure

```
0_install/                      package installation
1_functions/                    shared simulation and GWAS helpers
2_simulations/                  one script per module (a-e)
3_diagnostics/                  diagnostic methods
4_correction-individual-level/  corrections needing individual-level data
5_correction-summary-level/     corrections using summary statistics
6_vignettes/                    worked examples by module
7_figures/                      code for Figures 1 and 2
8_shiny/                        method-selection app
data/                           simulated data (not tracked; regenerate)
results/                        outputs
```

## Citation

If you use this code, please cite the paper (details to follow) and the
original method publications listed in each script header.
