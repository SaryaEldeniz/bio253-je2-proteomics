# Run from the repository root or the Sarya_George project directory.
group_root <- if (file.exists("scripts/05_figures/09_render_residue_venn.R")) "." else "Sarya_George"
data_dir <- Sys.getenv("JE2_DATA_DIR")
data_config <- file.path(group_root, "config", "data_dir.txt")
if (!nzchar(data_dir) && file.exists(data_config)) {
  data_dir <- trimws(readLines(data_config, n = 1L, warn = FALSE)[1L])
}
if (is.na(data_dir) || !nzchar(data_dir)) data_dir <- file.path(group_root, "..", "data")

media <- c("TSB", "SN")
residues <- c("S", "T", "Y", "H")
genotypes <- c("stk", "stp")
contrasts <- as.vector(outer(genotypes, media, function(g, m) paste0(g, "_vs_anc_at_", m)))
phospho <- read.delim(file.path(data_dir, "derived", "dea_phospho_STYH.tsv.gz"))
tested <- subset(phospho, contrast %in% contrasts & !is_contaminant &
                   is.finite(FDR) & is.finite(diff))
stopifnot(setequal(tested$contrast, contrasts),
          !anyNA(tested$site_id), !anyDuplicated(tested[c("contrast", "site_id")]),
          all(tested$residue %in% residues))
hits <- subset(tested, estimate_type == "observed" & FDR < 0.10 & abs(diff) > 1)

count_exclusive_hits <- function(medium, residue) {
  sets <- lapply(genotypes, function(genotype) {
    contrast_name <- paste0(genotype, "_vs_anc_at_", medium)
    hits$site_id[hits$contrast == contrast_name & hits$residue == residue]
  })

  data.frame(
    medium = medium, residue = residue,
    stk_only_observed = length(setdiff(sets[[1]], sets[[2]])),
    stp_only_observed = length(setdiff(sets[[2]], sets[[1]]))
  )
}

counts <- do.call(rbind, lapply(media, function(medium) {
  do.call(rbind, lapply(residues, function(residue) count_exclusive_hits(medium, residue)))
}))

draw_residue_venn <- function(counts, output_path) {
  shown <- subset(counts, stk_only_observed + stp_only_observed > 0L)
  panels_by_medium <- split(shown, factor(shown$medium, levels = media), drop = TRUE)
  stopifnot(nrow(shown) > 0L)

  panel_counts <- vapply(panels_by_medium, nrow, integer(1))
  max_panels <- max(panel_counts)
  panel_layout <- matrix(0L, nrow = length(panel_counts) + 1L, ncol = 2L * max_panels)
  next_panel <- 1L
  for (i in seq_along(panel_counts)) {
    columns <- max_panels - panel_counts[i] + seq_len(2L * panel_counts[i])
    panel_layout[i, columns] <- rep(seq.int(next_panel, length.out = panel_counts[i]), each = 2L)
    next_panel <- next_panel + panel_counts[i]
  }
  panel_layout[nrow(panel_layout), ] <- next_panel

  ink <- "#152D40"
  observed_colour <- "#0B5278"
  png(output_path, width = 1800, height = 1120, res = 175, bg = "white")
  on.exit(invisible(dev.off()), add = TRUE)
  layout(panel_layout, heights = c(rep(1, length(panel_counts)), 0.23))
  par(oma = c(0.3, 0.5, 4.4, 0.5), family = "sans")

  for (medium_panels in panels_by_medium) {
    for (i in seq_len(nrow(medium_panels))) {
      panel <- medium_panels[i, ]
      par(mar = c(0.6, 0.3, 1.9, 0.3))
      plot.new()
      plot.window(xlim = c(-1.82, 1.82), ylim = c(-1.12, 1.46), asp = 1)
      title(main = paste0(panel$medium, "  /  p", panel$residue), line = 0.3, col.main = ink)
      symbols(c(-0.54, 0.54), c(0, 0), circles = c(1.00, 1.00), inches = FALSE,
              add = TRUE, bg = adjustcolor("#CBD8E2", alpha.f = 0.35), fg = "#718394")
      text(c(-0.82, 0.82), 1.18, c("stk vs anc", "stp vs anc"), col = ink, cex = 0.86)

      region_x <- c(-1.02, 1.02)
      text(region_x, 0.12, c(panel$stk_only_observed, panel$stp_only_observed),
           col = observed_colour, font = 2, cex = 1.6)
      text(region_x, -0.22, c("Only stk", "Only stp"), col = ink, cex = 0.81)
    }
  }

  par(mar = rep(0, 4))
  plot.new()
  text(0.5, 0.72, "Source observed only. Shared sites are omitted.",
       col = observed_colour, font = 2)
  text(0.5, 0.21, "FDR < 0.10 and |log2 fold change| > 1. Circle areas are schematic.",
       cex = 0.89, col = "#495866")
  mtext("Exclusive significant STYH phosphosites by residue", side = 3, outer = TRUE,
        line = 2.2, font = 2, cex = 1.3, col = ink)
  mtext("Each set compares a mutant against the ancestor within the same medium.",
        side = 3, outer = TRUE, line = 0.65, cex = 0.91, col = "#495866")
}

report_dir <- file.path(group_root, "reports")
figure_dir <- file.path(group_root, "figures", "step09")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.table(counts, file.path(report_dir, "step09_phospho_residue_venn_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
draw_residue_venn(counts, file.path(figure_dir, "figure09_phospho_residue_venn.png"))
