# Icon-led methods variant: schematic artwork only, no sample-data analysis.
group_root <- if (file.exists("Sarya_George.Rproj")) "." else "Sarya_George"
tryCatch(source(file.path(group_root, "R/setup.R")), error = function(e) {
  message("Using base R for vector artwork only.")
})
step_root <- file.path(group_root, "steps/10_integration_poster")
out <- file.path(step_root, "figures/methods_svg")
needed <- c("ink", "secondary", "blue", "teal", "mint", "pale", "pale_blue",
            "esc", "txt", "rect", "circle", "path", "group", "arrow", "defs",
            "database_icon", "contrast_icon", "enrich_icon", "network")
expressions <- parse(file.path(step_root, "scripts/methods_illustrations/01_render_design_analysis_svgs.R"))
for (x in as.list(expressions)) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% needed) eval(x)
}
stopifnot(all(vapply(needed, exists, logical(1), envir = globalenv(), inherits = FALSE)))
p <- "analysis-visual"
marker <- 'marker-end="url(#analysis-visual-arrow)"'

a <- c(
  txt(20, 28, "Project analysis", 25, weight = 700),
  txt(950, 27, "TSB vs SN", 16, secondary, anchor = "end"),
  rect(205, 67, 555, 135, pale, r = 22),
  rect(205, 211, 555, 121, pale_blue, r = 22),
  # Input mapping and the shared split into two data layers.
  rect(20, 94, 132, 210, pale_blue, r = 22),
  txt(86, 119, "ID mapping", 17, weight = 700, anchor = "middle"),
  database_icon(54, 135, .75),
  txt(86, 224, "JE2 IDs", 14, anchor = "middle"),
  path("M 86 230 L 86 243", teal, sw = 1.5, extra = marker),
  rect(49, 248, 74, 24, "#FFFFFF", blue, r = 8, sw = 1.3),
  txt(86, 265, "UniProt", 14, blue, 700, "middle"),
  txt(86, 291, "anc / stk / stp", 12, secondary, anchor = "middle"),
  path("M 155 178 L 178 178 L 178 129 L 201 129", teal, sw = 2, extra = marker),
  path("M 155 190 L 178 190 L 178 269 L 201 269", teal, sw = 2, extra = marker),
  # Protein branch: contrasts feed enrichment and independently select STRING nodes.
  txt(280, 92, "Total proteins", 18, weight = 700, anchor = "middle"),
  contrast_icon(216, 112, .56),
  txt(274, 126, "stk vs anc", 13),
  txt(274, 142, "stp vs anc", 13),
  txt(274, 158, "stp vs stk", 13),
  txt(280, 191, "Interactions", 14, secondary, anchor = "middle"),
  arrow(353, 137, 418, 137, p, sw = 2.5),
  txt(485, 92, "GO / KEGG", 18, weight = 700, anchor = "middle"),
  enrich_icon(445, 108, .72),
  txt(485, 191, "ORA / GSEA", 14, secondary, anchor = "middle"),
  path("M 280 65 L 280 57 Q 280 52 286 52 L 679 52 Q 685 52 685 58 L 685 65", teal, sw = 1.8, extra = marker),
  txt(485, 44, "Contrast-selected", 12, teal, 700, "middle"),
  txt(685, 92, "STRING", 18, weight = 700, anchor = "middle"),
  network(634, 101, .5),
  txt(685, 191, "Stk1 / Stp1", 14, secondary, anchor = "middle"),
  # Phosphosite branch uses within-medium comparisons, not interaction-selected sets.
  txt(280, 235, "Phosphosites", 18, weight = 700, anchor = "middle"),
  path("M 237 267 C 250 249 258 290 275 271 S 306 249 323 274", blue, sw = 7),
  circle(278, 265, 12, mint, teal, 2),
  txt(278, 270, "P", 13, ink, 700, "middle"),
  txt(280, 308, "Mutant vs anc", 14, secondary, anchor = "middle"),
  txt(280, 325, "Within TSB / SN", 11, secondary, anchor = "middle"),
  arrow(353, 272, 418, 272, p, sw = 2.5),
  txt(485, 235, "Overlap", 18, weight = 700, anchor = "middle"),
  circle(466, 272, 25, "#DDEEF5", blue, 2),
  circle(501, 272, 25, "#BAE2DA", teal, 2),
  txt(485, 316, "Venn + residues", 14, secondary, anchor = "middle"),
  arrow(548, 272, 621, 272, p, sw = 2.5),
  txt(685, 235, "Sequences", 18, weight = 700, anchor = "middle"),
  # User-supplied A/T/S/P methods symbol; heights are schematic, not frequencies.
  txt(628, 296, "A", 35, blue, 700),
  txt(653, 296, "T", 47, teal, 700),
  txt(680, 296, "S", 64, ink, 700),
  txt(723, 296, "P", 36, blue, 700),
  txt(685, 316, "JE2 logos", 14, secondary, anchor = "middle"),
  # Convergence is comparison/prioritisation, not a validated direct-target claim.
  path("M 763 137 Q 784 137 797 158", teal, sw = 2.5, extra = marker),
  path("M 763 272 Q 784 272 797 224", teal, sw = 2.5, extra = marker),
  circle(862, 190, 61, pale_blue, "none"),
  txt(862, 110, "Protein + site", 18, weight = 700, anchor = "middle"),
  circle(852, 180, 29, "#FFFFFF", teal, 4),
  path("M 875 202 L 903 230", teal, sw = 8),
  path("M 833 181 C 842 170 847 193 856 183 S 866 175 872 183", blue, sw = 4.5),
  circle(854, 178, 7, mint, teal, 1.4),
  txt(854, 181, "P", 8, ink, 700, "middle"),
  txt(862, 278, "Candidates", 19, teal, 700, "middle"),
  txt(862, 299, "Compare signals", 12, secondary, anchor = "middle"),
  txt(950, 346, "Network edges: functional associations", 11, secondary, anchor = "end")
)
svg <- paste0('<svg xmlns="http://www.w3.org/2000/svg" width="970" height="354" viewBox="0 0 970 354" role="img" aria-labelledby="',p,'-title ',p,'-desc">',
  '<title id="',p,'-title">Visual project analysis methods</title><desc id="',p,'-desc">Icon-led version of the existing workflow. Protein interactions compare stk versus anc, stp versus anc and stp versus stk. They feed enrichment; network nodes are independently selected from contrasts. Phosphosite mutant-versus-ancestor comparisons within media feed overlap, residue and JE2 sequence summaries, alongside separate interaction tests in the analysis. Protein/site signals are compared for prioritisation, without claiming completed occupancy correction or direct substrates. All icons, plots and networks are schematic method symbols. The A/T/S/P letters reproduce the requested methods icon: their heights are illustrative, not measured frequencies or a data-derived sequence motif. Detailed settings remain in the original methods panel and documentation.</desc>',defs(p),paste(a,collapse=""),'</svg>')
output_path <- file.path(out, "03_analysis_workflow_visual.svg")
writeLines(enc2utf8(svg), output_path, useBytes = TRUE)
message(normalizePath(output_path, winslash = "/"))
