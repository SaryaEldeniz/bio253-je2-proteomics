# Uses the axes, significance colours and model symbols from plot_step04_volcano.
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
  colours <- c("grey75", "#0072B2", "#D55E00")
  symbols <- c(16, 4, 1)

  png(output_path, width = 2400, height = 1090, res = 175)
  on.exit(invisible(dev.off()), add = TRUE)
  layout(rbind(1:3, rep(4L, 3)), heights = c(1, 0.18))
  par(oma = c(0.3, 2.6, 4.1, 0.7), family = "sans")

  for (i in seq_along(contrast_names)) {
    panel <- subset(plot_data, contrast == contrast_names[i])
    panel$n <- 1L
    bins <- aggregate(n ~ x_bin + y_bin + effect_class + model, panel, sum)
    bins <- bins[order(bins$effect_class, bins$model), ]
    stopifnot(sum(bins$n) == nrow(panel))

    par(mar = c(3.7, 2.5, 3.7, 0.9))
    plot(NA, xlim = x_limits, ylim = y_limits, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
    abline(v = c(-log2fc_cutoff, log2fc_cutoff), h = -log10(fdr_cutoff),
           col = "grey55", lty = 2)
    points(bins$x_bin, bins$y_bin, pch = symbols[bins$model],
           col = adjustcolor(colours[bins$effect_class], alpha.f = 0.76),
           cex = 0.66 + 0.22 * pmin(log2(bins$n), 5), lwd = 1.1)
    axis(1, tck = -0.018, cex.axis = 0.86, col.axis = "#273747")
    axis(2, las = 1, tck = -0.018, cex.axis = 0.86, col.axis = "#273747")
    box(bty = "l", col = "#273747")
    title(main = titles[i], adj = 0, line = 2.0, col.main = "#152D40")
    mtext(sprintf("%d significant of %s phosphosites", sum(panel$significant),
                  format(nrow(panel), big.mark = ",")),
          side = 3, adj = 0, line = 0.7, cex = 0.82, col = "#495866")
    mtext(expression("Interaction log"[2] * " fold change"), side = 1, line = 2.5, cex = 0.88)
  }

  par(mar = rep(0, 4))
  plot.new()
  legend(0.5, 0.72, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n", cex = 0.88,
         legend = c("Below thresholds", "Significant, negative effect", "Significant, positive effect"),
         col = colours, pch = 16)
  legend(0.5, 0.22, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n", cex = 0.88,
         legend = c("Source observed", "LOD-flagged model"), col = "#273747", pch = symbols[1:2])
  mtext("Phospho-STYH interactions for TSB versus SN", side = 3, outer = TRUE,
        line = 1.9, adj = 0, font = 2, cex = 1.24, col = "#152D40")
  mtext(expression(-log[10](FDR)), side = 2, outer = TRUE, line = 0.4, cex = 0.95)
}
