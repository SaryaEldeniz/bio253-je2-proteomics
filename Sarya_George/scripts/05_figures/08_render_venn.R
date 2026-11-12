# Run from the repository root or from the Sarya_George project directory.
group_root <- if (file.exists("scripts/05_figures/08_render_venn.R")) "." else "Sarya_George"
data_dir <- Sys.getenv("JE2_DATA_DIR")
data_config <- file.path(group_root, "config", "data_dir.txt")
if (!nzchar(data_dir) && file.exists(data_config)) {
  data_dir <- trimws(readLines(data_config, n = 1L, warn = FALSE)[1L])
}
if (is.na(data_dir) || !nzchar(data_dir)) data_dir <- file.path(group_root, "..", "data")

phospho <- read.delim(file.path(data_dir, "derived", "dea_phospho_STYH.tsv.gz"))
fdr_cutoff <- 0.10
log2fc_cutoff <- 1

count_overlap <- function(medium) {
  significant_sites <- function(genotype) {
    contrast_name <- paste0(genotype, "_vs_anc_at_", medium)
    estimates <- subset(phospho, contrast == contrast_name & !is_contaminant)
    stopifnot(nrow(estimates) > 0L, !anyNA(estimates$site_id), !anyDuplicated(estimates$site_id))
    subset(estimates, is.finite(FDR) & is.finite(diff) & FDR < fdr_cutoff & abs(diff) > log2fc_cutoff)
  }

  stk <- significant_sites("stk")
  stp <- significant_sites("stp")
  shared <- intersect(stk$site_id, stp$site_id)
  stk_diff <- stk$diff[match(shared, stk$site_id)]
  stp_diff <- stp$diff[match(shared, stp$site_id)]
  same <- sum(sign(stk_diff) == sign(stp_diff))

  data.frame(
    medium = medium,
    stk_only = nrow(stk) - length(shared),
    both = length(shared),
    stp_only = nrow(stp) - length(shared),
    same_direction = same,
    opposite_direction = length(shared) - same
  )
}

draw_venn <- function(result, output_path) {
  colours <- c("#0072B2", "#D55E00")
  ink <- "#152D40"
  png(output_path, width = 1900, height = 860, res = 170)
  on.exit(invisible(dev.off()), add = TRUE)
  par(mfrow = c(1, 2), mar = c(0.8, 0.8, 2, 0.8), oma = c(2.6, 1, 3.5, 1),
      family = "sans", lwd = 1.5)

  for (i in seq_len(nrow(result))) {
    plot.new()
    plot.window(xlim = c(-2.25, 2.25), ylim = c(-1.35, 1.6), asp = 1)
    title(main = result$medium[i], line = 0.6, col.main = ink)
    symbols(c(-0.65, 0.65), c(0, 0), circles = c(1.12, 1.12), inches = FALSE,
            add = TRUE, bg = adjustcolor(colours, alpha.f = 0.24), fg = colours)
    text(c(-0.86, 0.86), 1.29, c("stk vs anc", "stp vs anc"), col = colours, font = 2)

    region_x <- c(-1.16, 0, 1.16)
    counts <- unlist(result[i, c("stk_only", "both", "stp_only")], use.names = FALSE)
    text(region_x, 0.12, counts, cex = 1.7, font = 2, col = ink)
    text(region_x, -0.23, c("Only stk", "Both", "Only stp"), cex = 0.83, col = ink)
  }

  mtext("Overlap of significant STYH phosphosites", side = 3, outer = TRUE,
        line = 1.8, font = 2, cex = 1.2, col = ink)
  mtext("Mutant-versus-ancestor comparisons within each medium", side = 3, outer = TRUE,
        line = 0.35, cex = 0.86, col = "#495866")
  mtext("FDR < 0.10 and |log2 fold change| > 1. Circle areas are schematic.",
        side = 1, outer = TRUE, line = 1, cex = 0.86, col = "#495866")
}

result <- rbind(count_overlap("TSB"), count_overlap("SN"))
print(result, row.names = FALSE)

report_dir <- file.path(group_root, "reports")
figure_dir <- file.path(group_root, "figures", "step08")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.table(result, file.path(report_dir, "step08_phospho_venn_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
draw_venn(result, file.path(figure_dir, "figure08_phospho_venn.png"))
