# Figure 2: DAGs of the selection structures addressed in this review
# 

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# nodes touch -- Graphviz switch to straight edges
fig2 = grViz('
digraph fig2 {
  graph [layout  neato, splines  true, overlap  true, fontname  "Arial",
         outputorder  edgesfirst]
  node  [fontname  "Arial", fontsize  12, shape  circle, style  filled,
         fillcolor  "#FFFFFF", color  "#4D4D4D", penwidth  1.2,
         width  0.45, fixedsize  true]
  edge  [color  "#4D4D4D", arrowsize  0.7, penwidth  1.1, fontname  "Arial",
         fontsize  9]

  // ---------- (a) Non-representative participation (module a) ----------
  ta [shape  plaintext, fixedsize  false, style  "", label  "(a) Non-representative participation", pos  "1.5,3.0!"]
  aG [label  "G", pos  "0,1.4!"]
  aY [label  "Y", pos  "1.2,1.4!"]
  aS [label  "S", shape  box, pos  "2.6,1.4!", fillcolor  "#F2F2F2", penwidth  1.8]
  aZ [label  "Z", pos  "1.3,2.35!"]
  aE [label  "E", pos  "2.9,2.35!"]
  aU [label  "U", pos  "1.9,0.3!", style  "filled,dashed", fillcolor  "#F7F7F7"]
  aG -> aY; aG -> aZ; aZ -> aY; aZ -> aS; aE -> aS; aU -> aY
  aY -> aS [label  "a2"]
  aG -> aU [dir  none, style  dashed, color  "#B35900"]

  // ---------- (b) Case-control secondary phenotype (module b) ----------
  tb [shape  plaintext, fixedsize  false, style  "", label  "(b) Case-control secondary phenotype", pos  "5.5,3.0!"]
  bG [label  "G", pos  "4.4,2.2!"]
  bY [label  "Y", pos  "6.2,2.2!"]
  bD [label  "D", shape  box, pos  "5.3,1.0!", fillcolor  "#F2F2F2", penwidth  1.8]
  bG -> bY; bG -> bD; bY -> bD
  bG -> bY [dir  none, style  dashed, color  "#B35900"]

  // ---------- (c) Index event: progression among cases (module c) ----------
  tc [shape  plaintext, fixedsize  false, style  "", label  "(c) Index event: progression among cases", pos  "9.2,3.0!"]
  cG [label  "G", pos  "8.0,1.9!"]
  cX [label  "X", shape  box, pos  "9.2,1.9!", fillcolor  "#F2F2F2", penwidth  1.8]
  cU [label  "U", pos  "10.4,1.9!", style  "filled,dashed", fillcolor  "#F7F7F7"]
  cY [label  "Y", pos  "9.2,0.4!"]
  cG -> cX; cU -> cX; cG -> cY; cU -> cY
  cG -> cU [dir  none, style  dashed, color  "#B35900"]

  // ---------- (d) Sex-differential participation (module d) ----------
  td [shape  plaintext, fixedsize  false, style  "", label  "(d) Sex-differential participation", pos  "3.2,-0.9!"]
  dG   [label  "G",   pos  "2.0,-1.8!"]
  dY   [label  "Y",   pos  "3.2,-1.8!"]
  dSex [label  "Sex", pos  "4.4,-1.8!", width  0.55]
  dS   [label  "S", shape  box, pos  "3.2,-3.1!", fillcolor  "#F2F2F2", penwidth  1.8]
  dG -> dY; dSex -> dY; dSex -> dS
  dG -> dS [label  "by sex"]
  dY -> dS [label  "by sex"]
  dG -> dSex [dir  none, style  dashed, color  "#B35900"]

  // ---------- (e) Survival to recruitment (module e) ----------
  te [shape  plaintext, fixedsize  false, style  "", label  "(e) Survival to recruitment", pos  "7.6,-0.9!"]
  eGP [label  <G<SUB>P</SUB>>, pos  "6.2,-1.6!"]
  eP  [label  <Y<SUB>P</SUB>>, pos  "6.2,-2.9!"]
  eG  [label  "G", pos  "7.4,-1.6!"]
  eY  [label  "Y", pos  "8.6,-1.6!"]
  eR  [label  "R", shape  box, pos  "8.0,-2.9!", fillcolor  "#F2F2F2", penwidth  1.8]
  eU  [label  "U", pos  "9.3,-2.9!", style  "filled,dashed", fillcolor  "#F7F7F7"]
  eGP -> eG; eGP -> eP; eG -> eY; eG -> eR; eU -> eR; eU -> eY
  eG -> eU [dir  none, style  dashed, color  "#B35900"]
}
')

fig2   # view in the RStudio Viewer

# ---- Export -----------------------------------------------------------------
svg_txt <- export_svg(fig2)
writeLines(svg_txt, "7_figures/Figure2.svg")
rsvg_pdf(charToRaw(svg_txt), "7_figures/Figure2.pdf")
rsvg_png(charToRaw(svg_txt), "7_figures/Figure2.png", width  3000)

# ---- Export -----------------------------------------------------------------
svg_txt  export_svg(fig2)
writeLines(svg_txt, "Figure2.svg")
rsvg_png(charToRaw(svg_txt), "Figure2.png", width  3000)
