# Run from the repository root or the Sarya_George project directory.
group_root <- if (file.exists("scripts/05_figures/09_render_central_residue_logo.R")) "." else "Sarya_George"
counts <- read.delim(file.path(group_root, "reports", "step09_phospho_residue_venn_counts.tsv"))
output_path <- file.path(group_root, "figures", "step09", "figure09_central_residue_logo.png")

draw_central_logo <- function(counts, output_path) {
  ink <- "#152D40"
  colours <- c(S = "#087F71", T = "#185F96", Y = "#974680")
  png(output_path, width = 1840, height = 800, res = 160, bg = "white")
  on.exit(invisible(dev.off()), add = TRUE)

  grid::grid.text("Phosphorylated residue at position 0", x = 0.07, y = 0.91, just = "left",
                  gp = grid::gpar(fontsize = 26, fontface = "bold", col = ink))
  grid::grid.text("Exclusive significant STYH hits labelled Source observed",
                  x = 0.07, y = 0.82, just = "left",
                  gp = grid::gpar(fontsize = 13, col = "#495866"))
  grid::grid.text("Frequency (%)", x = 0.025, y = 0.44, rot = 90,
                  gp = grid::gpar(fontsize = 11, col = ink))

  groups <- expand.grid(genotype = c("stk", "stp"), medium = c("TSB", "SN"),
                        stringsAsFactors = FALSE)
  panel_centres <- seq(0.19, 0.85, length.out = nrow(groups))
  for (i in seq_len(nrow(groups))) {
    medium <- groups$medium[i]
    genotype <- groups$genotype[i]
    rows <- counts[counts$medium == medium, ]
    values <- setNames(rows[[paste0(genotype, "_only_observed")]], rows$residue)
    total <- sum(values)
    grid::grid.text(paste(medium, "/", "Only", genotype), x = panel_centres[i], y = 0.72,
                    gp = grid::gpar(fontsize = 13, fontface = "bold", col = ink))
    grid::grid.text(paste0("n = ", total, " hits"), x = panel_centres[i], y = 0.665,
                    gp = grid::gpar(fontsize = 10, col = "#495866"))

    grid::pushViewport(grid::viewport(x = panel_centres[i], y = 0.40, width = 0.18,
                                      height = 0.42, yscale = c(0, 1)))
    grid::grid.rect(gp = grid::gpar(fill = "#F3F7FA", col = NA))
    if (i == 1L) {
      grid::grid.yaxis(at = c(0, 0.5, 1), label = c("0", "50", "100"),
                       gp = grid::gpar(fontsize = 9, col = ink))
    }

    fractions <- sort(values[values > 0] / total)
    top <- cumsum(fractions)
    for (j in seq_along(fractions)) {
      letter <- names(fractions)[j]
      fraction <- fractions[j]
      glyph <- grid::textGrob(letter, gp = grid::gpar(fontsize = 100, fontface = "bold"))
      glyph_height <- grid::convertHeight(grid::grobHeight(glyph), "inches", valueOnly = TRUE)
      band_height <- grid::convertHeight(grid::unit(0.94 * fraction, "native"),
                                         "inches", valueOnly = TRUE)
      grid::grid.text(letter, y = grid::unit(top[j] - fraction / 2, "native"),
                      gp = grid::gpar(fontsize = 100 * band_height / glyph_height,
                                      fontface = "bold", col = colours[letter]))
    }
    grid::popViewport()
    grid::grid.text("0", x = panel_centres[i], y = 0.155, gp = grid::gpar(fontsize = 10, col = ink))
  }

  grid::grid.text("Letter height represents the residue frequency within each group.",
                  x = 0.07, y = 0.075, just = "left",
                  gp = grid::gpar(fontsize = 10.5, col = "#495866"))
}

draw_central_logo(counts, output_path)
