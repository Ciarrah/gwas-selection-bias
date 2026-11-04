Module (b): Case-control ascertainment and secondary phenotypes
================

- [Scenario](#scenario)
- [1. The problem: naive GWAS in the case-control
  sample](#1-the-problem-naive-gwas-in-the-case-control-sample)
- [2. Diagnostics](#2-diagnostics)
- [3. Corrections](#3-corrections)
  - [Inverse probability weighting](#inverse-probability-weighting)
  - [Scaled marginal model](#scaled-marginal-model)
  - [SPREG](#spreg)
- [4. Summary](#4-summary)

``` r
library(ggplot2)
source("1_functions/gwas_helpers.R")
source("4_correction-individual-level/ipw.R")
source("4_correction-individual-level/spreg_prepare.R")

cols <- c("analysis", "snp_group", "rel_bias_causal", "mean_abs_bias",
          "coverage", "type1_error")
```

## Scenario

Figure 2(b). Genotype G affects a secondary phenotype Y and disease D,
and Y also affects D. Cases are over-sampled: the study has equal
numbers of cases and controls, although only 5% of the population has
the disease. Analysing Y in this sample conditions on D, which distorts
the SNP-Y associations:

- SNPs affecting only disease (`D_snps`) become associated with Y;
- SNPs affecting Y (`Y_snps`) are biased.

Two secondary phenotypes (Y1, Y2) have the same SNP effects on the
standardised scale. Only Y1 affects disease.

``` r
if (!file.exists("data/module_b.rds")) source("2_simulations/b_case_control.R")
d <- readRDS("data/module_b.rds")
p <- d$pheno

knitr::kable(data.frame(parameter = names(d$params),
                        value = format(unlist(d$params), scientific = FALSE,
                                       drop0trailing = TRUE)),
             row.names = FALSE)
```

| parameter  | value  |
|:-----------|:-------|
| n_pop      | 200000 |
| m          | 200    |
| n_causal_y | 20     |
| n_causal_d | 20     |
| h2_y       | 0.1    |
| h2_d       | 0.2    |
| g_yd       | 0.5    |
| r_y        | 0.5    |
| K          | 0.05   |
| n_cases    | 5000   |
| n_controls | 5000   |

## 1. The problem: naive GWAS in the case-control sample

``` r
naive <- run_gwas(p$Y1, d$G)

knitr::kable(compare_by_group(naive, d$truth, d$groups, "naive (Y1)")[, cols],
             digits = 3, row.names = FALSE)
```

| analysis   | snp_group | rel_bias_causal | mean_abs_bias | coverage | type1_error |
|:-----------|:----------|----------------:|--------------:|---------:|------------:|
| naive (Y1) | Y_snps    |           0.115 |         0.021 |    0.800 |         NaN |
| naive (Y1) | D_snps    |             NaN |         0.032 |    0.550 |       0.450 |
| naive (Y1) | null      |             NaN |         0.013 |    0.938 |       0.062 |

Disease-only SNPs show spurious associations with Y1 well above the
nominal 5%, and the Y-acting SNPs are biased.

## 2. Diagnostics

There is no diagnostic specific to this design: the sampling scheme is
known, so bias should be assumed whenever a secondary phenotype is
analysed in a case-control sample. Variants known to affect only the
disease can serve as negative control variants (see module c).

## 3. Corrections

All three corrections use the known population prevalence of the disease
(K = 0.05).

### Inverse probability weighting

Each case and control is weighted by its population share divided by its
sample share.

``` r
w    <- case_control_weights(p$D, d$K)
ipw1 <- ipw_gwas(p$Y1, d$G, w)

knitr::kable(data.frame(group = c("cases", "controls"),
                        weight = round(c(w[p$D == 1][1], w[p$D == 0][1]), 3)),
             row.names = FALSE)
```

| group    | weight |
|:---------|-------:|
| cases    |    0.1 |
| controls |    1.9 |

### Scaled marginal model

When several secondary phenotypes share SNP effects, the scaled marginal
model estimates one common effect across them by weighted GEE, gaining
precision over analysing Y1 alone.

``` r
smm <- smm_gwas(p[, c("Y1", "Y2")], d$G, w)
```

``` r
knitr::kable(rbind(
  compare_by_group(naive, d$truth, d$groups, "naive (Y1)"),
  compare_by_group(ipw1,  d$truth, d$groups, "IPW (Y1)"),
  compare_by_group(smm,   d$truth, d$groups, "SMM (Y1 + Y2)")
)[, cols], digits = 3, row.names = FALSE)
```

| analysis      | snp_group | rel_bias_causal | mean_abs_bias | coverage | type1_error |
|:--------------|:----------|----------------:|--------------:|---------:|------------:|
| naive (Y1)    | Y_snps    |           0.115 |         0.021 |    0.800 |         NaN |
| naive (Y1)    | D_snps    |             NaN |         0.032 |    0.550 |       0.450 |
| naive (Y1)    | null      |             NaN |         0.013 |    0.938 |       0.062 |
| IPW (Y1)      | Y_snps    |           0.064 |         0.022 |    0.850 |         NaN |
| IPW (Y1)      | D_snps    |             NaN |         0.020 |    0.900 |       0.100 |
| IPW (Y1)      | null      |             NaN |         0.017 |    0.963 |       0.038 |
| SMM (Y1 + Y2) | Y_snps    |           0.040 |         0.016 |    0.850 |         NaN |
| SMM (Y1 + Y2) | D_snps    |             NaN |         0.015 |    0.950 |       0.050 |
| SMM (Y1 + Y2) | null      |             NaN |         0.016 |    0.956 |       0.044 |

### SPREG

SPREG (Lin and Zeng 2009) is external software for 64-bit Linux. Run
`4_correction-individual-level/spreg_prepare.R` to write the input file,
then run SPREG as described in the README. Its results are shown here if
the output file exists.

``` r
res_file <- "data/spreg/spreg_output.out"
if (file.exists(res_file)) {
  spreg <- read_spreg(res_file, snp = colnames(d$G), n = nrow(p))
  knitr::kable(compare_by_group(spreg, d$truth, d$groups, "SPREG (Y1)")[, cols],
               digits = 3, row.names = FALSE)
} else {
  cat("SPREG output not found; see the README to run SPREG.\n")
}
```

    ## SPREG output not found; see the README to run SPREG.

## 4. Summary

``` r
grp <- rep(names(d$groups), lengths(d$groups))[order(unlist(d$groups))]

est <- function(res, label) {
  data.frame(analysis = label, snp_group = grp, truth = d$truth,
             beta = res$beta)
}

plot_df <- rbind(est(naive, "naive (Y1)"), est(ipw1, "IPW (Y1)"),
                 est(smm, "SMM (Y1 + Y2)"))
plot_df$analysis <- factor(plot_df$analysis, levels = unique(plot_df$analysis))

ggplot(plot_df, aes(truth, beta, colour = snp_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.7) +
  facet_wrap(~ analysis, nrow = 1) +
  labs(x = "True SNP effect", y = "Estimated SNP effect",
       colour = "SNP group") +
  theme_bw()
```

![](../vignette_sources/b_case_control_files/figure-gfm/summary-figure-1.png)<!-- -->

Points on the dashed line are unbiased. In the naive analysis,
disease-only SNPs (true effect 0) spread away from zero. Weighting by
the known sampling fractions brings them back, and the scaled marginal
model does so with the smallest errors by pooling the two secondary
phenotypes.

``` r
sessionInfo()
```

    ## R version 4.6.0 (2026-04-24)
    ## Platform: aarch64-apple-darwin23
    ## Running under: macOS Tahoe 26.7.1
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/Copenhagen
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_4.0.3
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6           xfun_0.59              htmlwidgets_1.6.4      lattice_0.22-9         vctrs_0.7.3           
    ##  [6] tools_4.6.0            generics_0.1.4         stats4_4.6.0           sandwich_3.1-2         tibble_3.3.1          
    ## [11] pkgconfig_2.0.3        indexevent_0.2.0       Matrix_1.7-5           data.table_1.18.2.1    RColorBrewer_1.1-3    
    ## [16] S7_0.2.2               lifecycle_1.0.5        compiler_4.6.0         farver_2.1.2           maxLik_1.6-10         
    ## [21] codetools_0.2-20       systemfit_1.1-30       carData_3.0-6          htmltools_0.5.9        yaml_2.3.12           
    ## [26] lazyeval_0.2.3         Formula_1.2-6          plotly_4.12.0          pillar_1.11.1          car_3.1-5             
    ## [31] tidyr_1.3.2            MASS_7.3-65            abind_1.4-8            mclust_6.1.3           tidyselect_1.2.1      
    ## [36] digest_0.6.39          mvtnorm_1.4-2          dplyr_1.2.1            purrr_1.2.2            labeling_0.4.3        
    ## [41] geepack_1.3.13         splines_4.6.0          miscTools_0.6-30       fastmap_1.2.0          grid_4.6.0            
    ## [46] cli_3.6.6              magrittr_2.0.5         broom_1.0.12           SlopeHunter_1.1.0      withr_3.0.3           
    ## [51] scales_1.4.0           backports_1.5.1        rmarkdown_2.31         httr_1.4.8             otel_0.2.0            
    ## [56] sampleSelection_1.2-16 rjags_4-18             zoo_1.8-15             VGAM_1.1-14            coda_0.19-4.1         
    ## [61] evaluate_1.0.5         knitr_1.51             lmtest_0.9-40          viridisLite_0.4.3      rlang_1.3.0           
    ## [66] glue_1.8.1             rstudioapi_0.18.0      jsonlite_2.0.0         ieugwasr_1.1.0.9000    R6_2.6.1
