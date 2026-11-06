# =============================================================================
# Figure 1: Decision framework for addressing selection bias in GWAS
# =============================================================================

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ---- Diagram ----------------------------------------------------------------
# Data requirements for each method are given in Table 1, not in the figure.
# Both exits share group = "exits" so they sit in one column on the right.
fig1  grViz('
digraph fig1 {
  graph [rankdir = TB, fontname = "Arial", nodesep = 0.25, ranksep = 0.4,
         newrank = true, compound = true]
  node  [shape = box, style = "rounded,filled", fontname = "Arial", fontsize = 11,
         color = "#4D4D4D", penwidth = 1.2, margin = "0.12,0.06"]
  edge  [color = "#4D4D4D", arrowsize = 0.7, fontname = "Arial", fontsize = 9]

  start [label = "GWAS in a potentially\\nselected sample", fillcolor = "#FFFFFF"]

  // Step 1: diagnostics
  subgraph cluster_diag {
    label = "Step 1: Diagnose (use in combination)"
    labeljust = "l"; fontname = "Arial"; fontsize = 12; style = "rounded,filled";
    fillcolor = "#EAF2FB"; color = "#9DBBD9"
    nc [label = "Negative controls",  fillcolor = "#FFFFFF"]
    po [label = "Proxy outcomes",     fillcolor = "#FFFFFF"]
    pg [label = "Participation GWAS", fillcolor = "#FFFFFF"]
  }

  // Decision fed by Step 1
  ev    [label = "Evidence of selection,\\nor plausible by design?", shape = diamond,
         fillcolor = "#FFFFFF", width = 3.4, height = 0.9, fixedsize = true]
  stop1 [label = "No evidence detected:\\nreport diagnostics and\\nremaining uncertainty",
         fillcolor = "#F2F2F2", group = "exits"]

  // Step 2: corrections
  subgraph cluster_corr {
    label = "Step 2: Correct (data requirements in Table 1)"
    labeljust = "l"; fontname = "Arial"; fontsize = 12; style = "rounded,filled";
    fillcolor = "#FDF1E6"; color = "#E0B48A"
    data [label = "What data\\nare available?", shape = diamond,
          fillcolor = "#FFFFFF", width = 2.2, height = 0.8, fixedsize = true]
    ind  [label = "Individual-level data", fillcolor = "#FFF8F0"]
    sum  [label = "Summary statistics\\n(or generated from\\nindividual-level data)",
          fillcolor = "#FFF8F0"]
    nrp  [label = "Non-representative\\nparticipation",        fillcolor = "#FFFFFF"]
    ccs  [label = "Case-control\\nsecondary phenotype",        fillcolor = "#FFFFFF"]
    sdp  [label = "Sex-differential\\nparticipation",          fillcolor = "#FFFFFF"]
    ieb  [label = "Index event\\n(progression\\namong cases)", fillcolor = "#FFFFFF"]
    ipw [label = "IPW",                                      fillcolor = "#F3F3F3"]
    hk  [label = "Heckman\\ncorrection",                     fillcolor = "#F3F3F3"]
    sp  [label = "SPREG",                                    fillcolor = "#F3F3F3"]
    ss  [label = "Sex-stratified\\nmeta-analysis",           fillcolor = "#F3F3F3"]
    idx [label = "DHO\\nSlope-Hunter\\nMR-Horse adaptation", fillcolor = "#F3F3F3"]
  }

  // Exit when no correction is defensible
  nodef [label = "No defensible correction:\\nsensitivity analysis and\\ntransparent reporting",
         fillcolor = "#F2F2F2", group = "exits"]

  // Flow
  {rank = same; nc; po; pg}
  nc = po = pg [style = invis]
  start = po [lhead = cluster_diag]
  po = ev    [ltail = cluster_diag]
  ev = stop1 [label = " No"]
  {rank = same; ev; stop1}
  ev = data  [label = " Yes", lhead = cluster_corr, weight = 20]

  data = ind
  data = sum
  ind = nrp
  ind = ccs
  sum = sdp
  sum = ieb
  nrp = ipw
  nrp = hk
  ccs = sp
  sdp = ss
  ieb = idx
  data = sp [style = invis, weight = 10]

  // Both exits in one column on the right; the second leaves the Step 2 box
  stop1 = nodef [style = invis, weight = 50]
  ss = nodef [ltail = cluster_corr, style = dashed,
               label = " Data or\\n assumptions\\n not met"]

  {rank = same; nrp; ccs; sdp; ieb}
  {rank = same; ipw; hk; sp; ss; idx; nodef}
}
')

fig1

# ---- Export -----------------------------------------------------------------
svg_txt = export_svg(fig1)
writeLines(svg_txt, "Figure1.svg")
rsvg_png(charToRaw(svg_txt), "Figure1.png", width = 3000)
