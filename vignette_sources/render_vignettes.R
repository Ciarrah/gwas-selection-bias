# =============================================================================
# render_vignettes.R
# Renders the vignette sources and places the walkthroughs in 6_vignettes/.
# Figures stay in vignette_sources/<name>_files/ and are linked from there.
# Run from the repository root:  source("vignette_sources/render_vignettes.R")
# =============================================================================

src <- "vignette_sources"
out <- "6_vignettes"
vignettes <- c("a_participation", "b_case_control", "c_index_event",
               "d_sex_differential", "e_survival")

dir.create(out, showWarnings = FALSE)

for (f in vignettes) {
  rmarkdown::render(file.path(src, paste0(f, ".Rmd")),
                    output_format = rmarkdown::github_document(
                      toc = TRUE, html_preview = FALSE))
  
  md_src <- file.path(src, paste0(f, ".md"))
  md     <- readLines(md_src)
  # Point image links at the figures in vignette_sources/
  md <- gsub(paste0("(\\]\\(|src=\")", f, "_files/"),
             paste0("\\1../", src, "/", f, "_files/"), md)
  writeLines(md, file.path(out, paste0(f, ".md")))
  file.remove(md_src)
}