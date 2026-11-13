# Run from the repository root or from the Sarya_George project directory.
group_root <- if (file.exists("scripts/05_figures/09_render_residue_summary.R")) "." else "Sarya_George"
data_dir <- Sys.getenv("JE2_DATA_DIR")
data_config <- file.path(group_root, "config", "data_dir.txt")
if (!nzchar(data_dir) && file.exists(data_config)) {
  data_dir <- trimws(readLines(data_config, n = 1L, warn = FALSE)[1L])
}
if (is.na(data_dir) || !nzchar(data_dir)) data_dir <- file.path(group_root, "..", "data")

contrasts <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB",
               "stk_vs_anc_at_SN", "stp_vs_anc_at_SN")
residues <- c("S", "T", "Y", "H")
estimate_types <- c("observed", "lod_imputed")
fdr_cutoff <- 0.10
log2fc_cutoff <- 1

phospho <- read.delim(file.path(data_dir, "derived", "dea_phospho_STYH.tsv.gz"))
tested <- subset(phospho, contrast %in% contrasts & !is_contaminant &
                   is.finite(diff) & is.finite(FDR))
stopifnot(setequal(tested$contrast, contrasts), !anyNA(tested$site_id),
          !anyDuplicated(tested[c("contrast", "site_id")]),
          all(tested$residue %in% residues), all(tested$estimate_type %in% estimate_types))
tested$contrast <- factor(tested$contrast, levels = contrasts)
tested$residue <- factor(tested$residue, levels = residues)
tested$estimate_type <- factor(tested$estimate_type, levels = estimate_types)
hits <- subset(tested, FDR < fdr_cutoff & abs(diff) > log2fc_cutoff)

counts <- as.data.frame(xtabs(~ residue + contrast, data = tested), stringsAsFactors = FALSE)
names(counts)[3] <- "tested"
hit_counts <- xtabs(~ residue + contrast + estimate_type, data = hits)
counts$observed <- as.integer(hit_counts[, , "observed"])
counts$lod_imputed <- as.integer(hit_counts[, , "lod_imputed"])
counts$significant <- counts$observed + counts$lod_imputed
counts <- counts[c("contrast", "residue", "tested", "significant", "observed", "lod_imputed")]

draw_residue_counts <- function(counts, output_path) {
  titles <- c("stk vs anc in TSB", "stp vs anc in TSB", "stk vs anc in SN", "stp vs anc in SN")
  colours <- c("#0B5278", "#E59A47")
  ink <- "#152D40"
  png(output_path, width = 1640, height = 1230, res = 165, bg = "white")
  on.exit(invisible(dev.off()), add = TRUE)
  layout(rbind(1:2, 3:4, c(5, 5)), heights = c(1, 1, 0.26))
  par(oma = c(0, 2.5, 4.5, 0.7), family = "sans")
  axis_max <- max(10, 10 * ceiling((max(counts$significant) + 3) / 10))

  for (i in seq_along(contrasts)) {
    panel <- counts[counts$contrast == contrasts[i], ]
    par(mar = c(2.6, 2.7, 2.6, 0.8))
    centers <- barplot(t(as.matrix(panel[c("observed", "lod_imputed")])),
                       names.arg = paste0("p", panel$residue), col = colours,
                       border = NA, ylim = c(0, axis_max), axes = FALSE, cex.names = 0.95)
    axis(2, at = seq(0, axis_max, 10), las = 1, tck = -0.02, cex.axis = 0.84, col = ink)
    abline(h = 0, col = ink, xpd = FALSE)
    text(centers, panel$significant, labels = panel$significant, pos = 3, offset = 0.28,
         col = ink, cex = 0.91)
    title(main = titles[i], adj = 0, line = 0.55, col.main = ink, cex.main = 1.03)
  }

  par(mar = rep(0, 4))
  plot.new()
  legend(0.5, 0.68, legend = c("Source observed", "LOD-imputed"), fill = colours,
         horiz = TRUE, bty = "n", border = NA, xjust = 0.5, yjust = 0.5)
  text(0.5, 0.21, sprintf("FDR < %.2f and |log2 fold change| > %g.", fdr_cutoff, log2fc_cutoff),
       cex = 0.89, col = "#495866")
  mtext("Significant STYH phosphosites by residue", side = 3, outer = TRUE,
        line = 2.25, font = 2, cex = 1.32, col = ink)
  mtext("Mutant-versus-ancestor contrasts within TSB and SN", side = 3, outer = TRUE,
        line = 0.75, cex = 0.94, col = "#495866")
  mtext("Number of significant sites", side = 2, outer = TRUE, line = 0.45, col = ink)
}

draw_synthetic_logo <- function(output_path) {
  # A fixed toy alignment made for this illustration.
  sequences <- c("AVTASGLVK", "GLTASGIVR", "VLTGSGLLK", "ALVASALVR",
                 "GMVGSGIVK", "VMLGSGLLR", "AVVASALVK", "GLTGSGIVR",
                 "VMTGSAILK", "ALVGSGVVR", "GLLASGLLK", "VMVGSALVR")
  alignment <- do.call(rbind, strsplit(sequences, "", fixed = TRUE))
  positions <- seq_len(ncol(alignment)) - 5L
  palette <- c(A = "#286694", G = "#457b8f", I = "#427a51", K = "#8556a0",
               L = "#264d7a", M = "#287d62", R = "#9d5587", S = "#b65c14",
               T = "#008179", V = "#526c35")

  png(output_path, width = 1820, height = 670, res = 165, bg = "white")
  on.exit(invisible(dev.off()), add = TRUE)
  ink <- "#152D40"
  grid::grid.text("How to read a sequence logo", x = 0.08, y = 0.91, just = "left",
                  gp = grid::gpar(fontsize = 26, fontface = "bold", col = ink))
  grid::grid.text("Synthetic alignment for methods illustration", x = 0.08, y = 0.79, just = "left",
                  gp = grid::gpar(fontsize = 12, fontface = "bold", col = "#8E4C1E"))

  grid::pushViewport(grid::viewport(x = 0.54, y = 0.48, width = 0.82, height = 0.26,
                                    xscale = c(-4.5, 4.5), yscale = c(0, 1)))
  grid::grid.rect(x = grid::unit(0, "native"), width = grid::unit(0.92, "native"),
                  gp = grid::gpar(fill = "#FFF0DF", col = NA))
  grid::grid.yaxis(at = c(0, 0.25, 0.5, 0.75, 1), label = c("0", "25", "50", "75", "100"),
                   gp = grid::gpar(fontsize = 10, col = ink))
  grid::grid.xaxis(at = positions, gp = grid::gpar(fontsize = 10, col = ink))
  grid::grid.text("Frequency (%)", x = grid::unit(-0.10, "npc"), rot = 90,
                  gp = grid::gpar(fontsize = 11, col = ink))

  for (i in seq_along(positions)) {
    proportions <- sort(prop.table(table(alignment[, i])))
    lower_edge <- 0
    for (j in seq_along(proportions)) {
      letter <- names(proportions)[j]
      fraction <- as.numeric(proportions[j])
      # Size each glyph to fit its frequency band and its column.
      measure <- grid::textGrob(letter, gp = grid::gpar(fontsize = 100, fontface = "bold"))
      glyph_height <- grid::convertHeight(grid::grobHeight(measure), "inches", valueOnly = TRUE)
      glyph_width <- grid::convertWidth(grid::grobWidth(measure), "inches", valueOnly = TRUE)
      band_height <- grid::convertHeight(grid::unit(0.96 * fraction, "native"), "inches", valueOnly = TRUE)
      column_width <- grid::convertWidth(grid::unit(0.96, "native"), "inches", valueOnly = TRUE)
      fontsize <- 100 * band_height / glyph_height
      stopifnot(glyph_width * fontsize / 100 <= column_width)
      grid::grid.text(letter, x = grid::unit(positions[i], "native"),
                      y = grid::unit(lower_edge + fraction / 2, "native"),
                      gp = grid::gpar(fontsize = fontsize, fontface = "bold", col = palette[letter]))
      lower_edge <- lower_edge + fraction
    }
  }
  grid::popViewport()
  grid::grid.text("Position in the alignment", x = 0.54, y = 0.19,
                  gp = grid::gpar(fontsize = 12, col = ink))
  grid::grid.text("Letter height represents frequency among 12 invented sequences. Position 0 is the aligned centre.",
                  x = 0.08, y = 0.063, just = "left",
                  gp = grid::gpar(fontsize = 10.5, col = "#495866"))
}

report_dir <- file.path(group_root, "reports")
figure_dir <- file.path(group_root, "figures", "step09")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.table(counts, file.path(report_dir, "step09_phospho_residue_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
draw_residue_counts(counts, file.path(figure_dir, "figure09_phospho_residue_counts.png"))
draw_synthetic_logo(file.path(figure_dir, "figure09_synthetic_sequence_logo.png"))
print(counts, row.names = FALSE)
