plot_experimental_design <- function(layer_samples) {
  plot_data <- layer_samples |>
    mutate(layer = factor(layer, levels = c("total", "phospho_STYH"),
                          labels = c("Total proteome", "STYH phospho")),
           x = match(as.character(medium), MEDIA) + 0.18 * (replicate - 2),
           y = match(as.character(genotype), GENOTYPES))
  ggplot(plot_data, aes(x, y, colour = genotype)) +
    geom_point(aes(shape = reinjected), size = 3.4) +
    geom_text(aes(label = replicate), vjust = 1.9, size = 2.8, show.legend = FALSE) +
    facet_wrap(~layer, nrow = 1) +
    scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                       breaks = c(FALSE, TRUE), labels = c("No", "Yes")) +
    scale_x_continuous(breaks = 1:3, limits = c(0.5, 3.5),
                       labels = c("TSB\nControl", "SNpulse\n15 min PASN", "SN\n~15 h PASN")) +
    scale_y_reverse(breaks = 1:3, limits = c(3.45, 0.55),
                    labels = c("Ancestor", "stk1::Tn", "stp1::Tn")) +
    labs(title = "JE2 experimental design",
         subtitle = "3 replicate labels per group; 27 matching sample labels in each measurement layer",
         x = NULL, y = NULL, shape = "Re-injected acquisition",
         caption = paste("Points represent replicate labels; triangles mark the two re-injected phospho acquisitions.",
                         "Re-injections are not additional experimental replicates. PASN is cell-free PAO1 spent medium.", sep = "\n")) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.caption = element_text(hjust = 0),
          strip.background = element_rect(fill = "grey95", colour = "grey60"))
}
