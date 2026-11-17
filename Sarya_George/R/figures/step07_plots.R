# Compact volcano panels with arrows for significant source-observed estimates.
plot_step07_volcano <- function(phospho, contrast_names, titles, output_path,
                                fdr_cutoff = 0.10, log2fc_cutoff = 1) {
  plot_data <- subset(
    phospho,
    contrast %in% contrast_names & !is_contaminant & is.finite(diff) & is.finite(FDR) &
      FDR >= 0 & FDR <= 1
  )
  stopifnot(all(contrast_names %in% plot_data$contrast), length(titles) == length(contrast_names))

  # Preserve the source FDR and use a finite plotting coordinate for zeros.
  positive_fdr <- plot_data$FDR[plot_data$FDR > 0]
  display_floor <- if (length(positive_fdr)) min(positive_fdr) / 10 else .Machine$double.xmin
  display_floor <- max(display_floor, .Machine$double.xmin)
  plot_data$log_fdr <- -log10(pmax(plot_data$FDR, display_floor))
  plot_data$significant <- with(plot_data, FDR < fdr_cutoff & abs(diff) > log2fc_cutoff)
  plot_data$effect_class <- with(plot_data, 1L + as.integer(significant) * (1L + (diff > 0)))
  plot_data$model <- match(plot_data$estimate_type, c("observed", "lod_imputed"))
  plot_data$model[is.na(plot_data$model)] <- 3L

  x_bin_width <- 0.5
  y_bin_width <- 0.25
  plot_data$x_bin <- floor(plot_data$diff / x_bin_width) * x_bin_width + x_bin_width / 2
  plot_data$y_bin <- floor(plot_data$log_fdr / y_bin_width) * y_bin_width + y_bin_width / 2
  x_limits <- c(floor(min(plot_data$x_bin)), ceiling(max(plot_data$x_bin)))
  y_limits <- c(0, ceiling(max(plot_data$y_bin)))
  ink <- "#132F44"
  muted <- "#4F6479"
  axis_colour <- "#B0BFCB"
  guide_colour <- "#B9CBD9"
  colours <- c("#18546C", "#E2753A", "#A6C638")

  png(output_path, width = 2250, height = 820, res = 185, bg = "white")
  on.exit(invisible(dev.off()), add = TRUE)
  layout(rbind(1:3, rep(4L, 3)), heights = c(1, 0.16))
  par(oma = c(0.2, 2.6, 3.6, 0.5), family = "sans")
  x_ticks <- seq(4 * ceiling(x_limits[1] / 4), 4 * floor(x_limits[2] / 4), by = 4)
  y_ticks <- seq(0, y_limits[2], by = 3)

  for (i in seq_along(contrast_names)) {
    panel <- subset(plot_data, contrast == contrast_names[i])
    panel$n <- 1L
    bins <- aggregate(n ~ x_bin + y_bin + effect_class + model, panel, sum)
    bins <- bins[order(bins$effect_class, bins$model), ]
    stopifnot(sum(bins$n) == nrow(panel))

    par(mar = c(4.9, 2, 2.1, 1.1))
    plot(NA, xlim = x_limits, ylim = y_limits, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
    abline(v = 0, col = axis_colour)
    abline(v = c(-log2fc_cutoff, log2fc_cutoff), h = -log10(fdr_cutoff),
           col = guide_colour, lty = 3, lwd = 1)
    points(bins$x_bin, bins$y_bin, pch = 21,
           col = ink, bg = adjustcolor(colours[bins$effect_class], alpha.f = 0.82),
           cex = 0.52 + 0.12 * pmin(log2(bins$n), 4), lwd = 0.5)

    marked <- subset(bins, model == 1L & effect_class > 1L)
    if (nrow(marked)) {
      points(marked$x_bin, marked$y_bin, pch = 21, cex = 1.15, lwd = 1.2,
             col = ink, bg = colours[marked$effect_class])
      arrows(marked$x_bin, marked$y_bin + 1, marked$x_bin, marked$y_bin + 0.21,
             length = 0.08, angle = 24, col = ink, lwd = 1.5)
    }

    axis(1, at = x_ticks, tck = 0, line = 0, mgp = c(2, 0.5, 0),
         cex.axis = 0.96, col = axis_colour, col.axis = muted)
    axis(2, at = y_ticks, las = 1, tck = 0, lwd = 0, mgp = c(2, 0.4, 0),
         cex.axis = 0.96, col.axis = muted)
    title(main = titles[i], line = 0.7, cex.main = 1.26, col.main = ink)
    mtext(expression("Interaction effect (log"[2] * ")"), side = 1, line = 2.05,
          cex = 0.93, col = muted)
    mtext(sprintf("Selected hits  +%d / -%d     |     Source observed  %d",
                  sum(panel$effect_class == 3L), sum(panel$effect_class == 2L),
                  sum(marked$n)),
          side = 1, adj = 0, line = 3.65, cex = 0.76, col = muted)
  }

  par(mar = rep(0, 4))
  plot.new()
  legend(0.5, 0.73, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n", cex = 0.83,
         legend = c("Below thresholds", "Selected, negative", "Selected, positive"),
         col = ink, pt.bg = colours, pch = 21, text.col = muted)
  text(0.5, 0.15, "Arrows mark significant Source observed hits. Points show binned estimates.",
       cex = 0.75, col = muted)
  mtext("Phospho-STYH interactions for TSB versus SN", side = 3, outer = TRUE,
        line = 1.45, font = 2, cex = 1.25, col = ink)
  mtext(expression(-log[10] * "(phosphosite FDR)"), side = 2, outer = TRUE,
        line = 0.45, cex = 0.96, col = muted)
}
