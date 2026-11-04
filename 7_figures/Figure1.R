# =============================================================================
# Figure 1: Decision framework for addressing selection bias in GWAS
# =============================================================================

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ---- Diagram ----------------------------------------------------------------
fig1 = grViz('
digraph fig1 {
  graph [rankdir = TB, fontname = "Arial", nodesep = 0.22, ranksep = 0.4,
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

    nc [fillcolor = "#FFFFFF", label = <Negative controls<BR/>any selection<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: a negative</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>control phenotype</I></FONT>>]
    po [fillcolor = "#FFFFFF", label = <Proxy outcomes<BR/>survival to recruitment<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: familial</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>outcome data</I></FONT>>]
    pg [fillcolor = "#FFFFFF", label = <Participation GWAS<BR/>non-representative<BR/>participation<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: participation indicator</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>or population reference</I></FONT>>]
  }

  // Decision fed by Step 1
  ev    [label = "Evidence of selection,\\nor plausible by design?", shape = diamond,
         fillcolor = "#FFFFFF", margin = "0.03,0.03"]
  stop1 [label = "No: proceed and\\nreport diagnostics", fillcolor = "#F2F2F2"]

  // Step 2: corrections
  subgraph cluster_corr {
    label = "Step 2: Correct"
    labeljust = "l"; fontname = "Arial"; fontsize = 12; style = "rounded,filled";
    fillcolor = "#FDF1E6"; color = "#E0B48A"

    data [label = "What data\\nare available?", shape = diamond,
          fillcolor = "#FFFFFF", margin = "0.03,0.03"]
    ind  [label = "Individual-level data", fillcolor = "#FFF8F0"]
    sum  [label = "Summary statistics",    fillcolor = "#FFF8F0"]

    nrp  [label = "Non-representative\\nparticipation",          fillcolor = "#FFFFFF"]
    ccs  [label = "Case-control\\nsecondary phenotype",          fillcolor = "#FFFFFF"]
    sdp  [label = "Sex-differential\\nparticipation",            fillcolor = "#FFFFFF"]
    ieb  [label = "Index event\\n(progression\\namong cases)",   fillcolor = "#FFFFFF"]

    ipw [fillcolor = "#F3F3F3", label = <IPW<BR/>selection on<BR/>measured factors<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: population</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>reference data to</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>calibrate weights</I></FONT>>]
    hk  [fillcolor = "#F3F3F3", label = <Heckman correction<BR/>selection on<BR/>unmeasured factors<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: covariates on</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>non-responders; an</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>exclusion variable</I></FONT>>]
    sp  [fillcolor = "#F3F3F3", label = <SPREG<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: population</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>disease prevalence</I></FONT>>]
    ss  [fillcolor = "#F3F3F3", label = <Sex-stratified<BR/>GWAS<BR/>meta-analysis<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: sex-</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>stratified GWAS</I></FONT>>]
    idx [fillcolor = "#F3F3F3", label = <DHO<BR/>Slope-Hunter<BR/>MR-Horse<BR/>adaptation<BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>Needs: incidence</I></FONT><BR/><FONT POINT-SIZE="9" COLOR="#B35900"><I>GWAS</I></FONT>>]
  }

  // Flow
  {rank = same; nc; po; pg}
  nc = po = pg [style = invis]

  start = po [lhead = cluster_diag]
  po = ev    [ltail = cluster_diag]
  ev = stop1 [label = " No"]
  {rank = same; ev; stop1}
  ev = data  [label = " Yes", lhead = cluster_corr]

  data = ind
  data = sum
  ind  = sum [style = dashed, label = " Can generate", constraint = false]

  ind = nrp
  ind = ccs
  sum = sdp
  sum = ieb

  nrp = ipw
  nrp = hk
  ccs = sp
  sdp = ss
  ieb = idx

  // Centres the top layers over Step 2 (invisible anchor to the middle method)
  data = hk [style = invis, weight = 20]

  {rank = same; nrp; ccs; sdp; ieb}
  {rank = same; ipw; hk; sp; ss; idx}
}
')

fig1

# ---- Export -----------------------------------------------------------------
svg_txt = export_svg(fig1)
writeLines(svg_txt, "Figure1.svg")
rsvg_png(charToRaw(svg_txt), "Figure1.png", width = 3000)
