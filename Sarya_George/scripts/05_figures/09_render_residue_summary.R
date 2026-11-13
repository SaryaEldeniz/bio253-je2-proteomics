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
    bar_counts <- t(as.matrix(panel[c("observed", "lod_imputed")]))
    centers <- barplot(bar_counts,
                       names.arg = paste0("p", panel$residue), col = colours,
                       border = NA, ylim = c(0, axis_max), axes = FALSE, cex.names = 0.95)
    axis(2, at = seq(0, axis_max, 10), las = 1, tck = -0.02, cex.axis = 0.84, col = ink)
    abline(h = 0, col = ink, xpd = FALSE)
    observed <- panel$observed > 0L
    imputed <- panel$lod_imputed > 0L
    text(centers[observed], panel$observed[observed], labels = panel$observed[observed],
         pos = 3, offset = 0.20, col = ink, cex = 0.99, font = 2)
    text(centers[imputed], colSums(bar_counts)[imputed], labels = panel$lod_imputed[imputed],
         pos = 3, offset = 0.29, col = ink, cex = 0.99, font = 2)
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

report_dir <- file.path(group_root, "reports")
figure_dir <- file.path(group_root, "figures", "step09")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.table(counts, file.path(report_dir, "step09_phospho_residue_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
draw_residue_counts(counts, file.path(figure_dir, "figure09_phospho_residue_counts.png"))
print(counts, row.names = FALSE)
