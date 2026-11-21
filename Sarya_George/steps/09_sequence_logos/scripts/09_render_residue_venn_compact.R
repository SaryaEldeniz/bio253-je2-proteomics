# Run from the repository root or Sarya_George. No DEA is repeated.
group_root <- if (file.exists("steps/09_sequence_logos/reports/step09_phospho_residue_venn_counts.tsv")) "." else "Sarya_George"
step_root <- file.path(group_root, "steps/09_sequence_logos")
tryCatch(source(file.path(group_root, "R/setup.R")), error = function(e) {
  message("Project packages unavailable; using base R for this graphics-only export.")
})

counts <- read.delim(file.path(step_root, "reports/step09_phospho_residue_venn_counts.tsv"))
panels <- subset(counts, stk_only_observed + stp_only_observed > 0L)
stopifnot(identical(paste(panels$medium, panels$residue),
                    c("TSB S", "TSB T", "TSB Y", "SN S", "SN T")))

draw_compact_venn <- function(panels, output_path) {
  grDevices::svg(output_path, width = 6.5, height = 3.45,
                 pointsize = 11, bg = "white")
  on.exit(invisible(dev.off()), add = TRUE)
  par(mar = rep(0, 4), family = "sans", xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 468), ylim = c(0, 248.4), asp = 1)
  ink <- "#152D40"
  blue <- "#0B5278"
  text(234, 235, "Exclusive observed phosphosites", font = 2, cex = 1.3, col = ink)

  centres <- data.frame(x = c(78, 234, 390, 156, 312), y = c(168, 168, 168, 69, 69))
  for (i in seq_len(nrow(panels))) {
    p <- panels[i, ]
    x <- centres$x[i]
    y <- centres$y[i]
    text(x, y + 44, paste0(p$medium, " / p", p$residue),
         font = 2, cex = 1.1, col = ink)
    symbols(c(x - 18, x + 18), c(y, y), circles = c(33, 33),
            inches = FALSE, add = TRUE,
            bg = adjustcolor("#CBD8E2", alpha.f = 0.35),
            fg = "#718394", lwd = 0.8)
    text(c(x - 30, x + 30), y + 4,
         c(p$stk_only_observed, p$stp_only_observed),
         col = blue, font = 2, cex = 1.5)
    text(c(x - 30, x + 30), y - 14, c("stk", "stp"), col = ink, cex = 0.8)
  }
  text(234, 19,
       "Mutants vs ancestor; source observed only. Shared sites omitted; circle areas schematic.",
       cex = 0.73, col = "#495866")
  text(234, 7, "FDR < 0.10; |log2 fold change| > 1.", cex = 0.73, col = "#495866")
}

output_path <- file.path(step_root, "figures/figure09_phospho_residue_venn_compact.svg")
draw_compact_venn(panels, output_path)
message(normalizePath(output_path, winslash = "/"))
