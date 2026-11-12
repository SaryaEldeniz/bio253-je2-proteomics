# Step 4 figure helpers. The calling script supplies tables and saves figures.
# These functions do not read inputs, refit models, or fill missing measurements.

plot_step04_volcano <- function(x, contrast_names, titles = contrast_names) {
  required <- c("contrast", "diff", "FDR", "estimate_type", "selected_hit")
  if (!all(required %in% names(x))) stop("Volcano input is missing: ", paste(setdiff(required, names(x)), collapse = ", "))
  if (!length(contrast_names) || !all(contrast_names %in% x$contrast)) stop("Every requested contrast must occur in the volcano table.")
  if (length(titles) != length(contrast_names)) stop("Provide one title per requested contrast.")
  titles <- setNames(as.character(titles), contrast_names)
  plot_data <- x |>
    filter(contrast %in% contrast_names, is.finite(diff), is.finite(FDR), FDR >= 0, FDR <= 1)
  if (!nrow(plot_data)) stop("No finite effects and valid FDR values to plot.")
  # The original FDR remains unchanged. A zero is placed one log unit above the
  # smallest positive FDR for display, avoiding infinite plotting coordinates.
  positive_fdr <- plot_data$FDR[plot_data$FDR > 0]
  display_floor <- if (length(positive_fdr)) max(.Machine$double.xmin, min(positive_fdr) / 10) else .Machine$double.xmin
  has_zero <- any(plot_data$FDR == 0)
  plot_data <- plot_data |>
    mutate(panel = factor(contrast, levels = contrast_names, labels = unname(titles)),
           log_fdr = -log10(pmax(FDR, display_floor)),
           selection = factor(case_when(coalesce(selected_hit, FALSE) & diff < 0 ~ "Selected, negative effect",
                                        coalesce(selected_hit, FALSE) & diff > 0 ~ "Selected, positive effect",
                                        TRUE ~ "Not selected"),
                              levels = c("Not selected", "Selected, negative effect", "Selected, positive effect")),
           estimate_display = factor(case_when(estimate_type == "observed" ~ "Source observed",
                                                estimate_type == "lod_imputed" ~ "LOD-flagged model",
                                                TRUE ~ "Other / unknown"),
                                     levels = c("Source observed", "LOD-flagged model", "Other / unknown")))
  zero_note <- if (has_zero) " FDR zeros use a finite display floor; source values are unchanged." else ""
  ggplot(plot_data, aes(diff, log_fdr, colour = selection, shape = estimate_display)) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_point(size = 1.25, alpha = 0.65, stroke = 0.5) +
    facet_wrap(~panel, nrow = 1, drop = FALSE) +
    scale_colour_manual(values = c("Not selected" = "grey75", "Selected, negative effect" = "#0072B2",
                                   "Selected, positive effect" = "#D55E00"), drop = FALSE) +
    scale_shape_manual(values = c("Source observed" = 16, "LOD-flagged model" = 4, "Other / unknown" = 1), drop = TRUE) +
    labs(x = "Native contrast effect (log2, A minus B)", y = "-log10(FDR)", colour = NULL, shape = "Source model",
         caption = paste0("Selected: FDR < ", FDR_CUTOFF, " and |effect| > ", LFC_CUTOFF,
                          ". Model flags are global; inspect detections in the groups used by each contrast.", zero_note)) +
    theme(legend.position = "bottom", legend.box = "vertical", plot.caption = element_text(hjust = 0),
          panel.grid.minor = element_blank(), strip.text = element_text(size = 10))
}

plot_step04_hit_counts <- function(summary) {
  required <- c("contrast", "contrast_label", "direction", "n_hits_observed", "n_hits_lod_imputed", "n_hits_other")
  if (!all(required %in% names(summary))) stop("Hit-count input is missing: ", paste(setdiff(required, names(summary)), collapse = ", "))
  if (!nrow(summary)) stop("No contrast hit counts to plot.")
  # Each direction is explicitly attached to its own contrast; a common 'up/down'
  # label would obscure the different biological meaning of interaction signs.
  short_labels <- gsub(": PASN response difference", ": interaction", summary$contrast_label, fixed = TRUE)
  direction_labels <- stringr::str_wrap(summary$direction, width = 30)
  ordered_labels <- paste(short_labels, direction_labels, sep = "\n")
  plot_data <- summary |>
    mutate(bar_id = ordered_labels,
           bar_id = factor(bar_id, levels = rev(unique(ordered_labels)))) |>
    pivot_longer(c(n_hits_observed, n_hits_lod_imputed, n_hits_other),
                 names_to = "estimate", values_to = "n") |>
    mutate(estimate = factor(estimate, levels = c("n_hits_observed", "n_hits_lod_imputed", "n_hits_other"),
                             labels = c("Source observed", "LOD-flagged model", "Other / unknown")))
  ggplot(plot_data, aes(n, bar_id, fill = estimate)) +
    geom_col(width = 0.72) +
    geom_text(data = summary |> mutate(bar_id = factor(ordered_labels, levels = rev(unique(ordered_labels)))),
              mapping = aes(x = n_hits, y = bar_id, label = n_hits), inherit.aes = FALSE,
              hjust = -0.15, size = 3) +
    scale_fill_manual(values = c("Source observed" = "#0072B2", "LOD-flagged model" = "#E69F00", "Other / unknown" = "grey60"), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.09))) +
    labs(x = "Selected ordinary protein entries", y = NULL, fill = "Source model",
         caption = paste0("Course selection rule: FDR < ", FDR_CUTOFF, " and |effect| > ", LFC_CUTOFF,
                          ". The global model flag does not identify which contrast groups needed imputation.")) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))
}

plot_step04_samples <- function(abundances, protein_keys, labels = NULL) {
  required <- c("protein_key", "genotype", "medium", "replicate", "normalized_abundance")
  if (!all(required %in% names(abundances))) stop("Sample input is missing: ", paste(setdiff(required, names(abundances)), collapse = ", "))
  if (!length(protein_keys) || anyNA(protein_keys)) stop("Supply explicit protein keys.")
  protein_keys <- unique(protein_keys)
  x <- abundances |> filter(protein_key %in% protein_keys, medium %in% c("TSB", "SN"))
  if (!nrow(x)) stop("Requested proteins have no TSB/SN rows in the abundance table.")
  if (length(protein_keys) > 24) stop("Use at most 24 explicit proteins per sample figure.")
  if (is.null(labels)) {
    label_map <- x |> distinct(protein_key)
    if ("locus_tag" %in% names(x)) {
      label_map <- x |> distinct(protein_key, locus_tag) |>
        mutate(label = coalesce(locus_tag, protein_key)) |> select(protein_key, label)
    } else label_map <- mutate(label_map, label = protein_key)
  } else if (is.character(labels) && !is.null(names(labels))) {
    label_map <- tibble(protein_key = names(labels), label = unname(labels))
  } else if (is.data.frame(labels) && all(c("protein_key", "label") %in% names(labels))) {
    label_map <- as_tibble(labels) |> select(protein_key, label)
  } else stop("labels must be a named character vector or a table with protein_key and label.")
  if (anyDuplicated(label_map$protein_key)) stop("Each protein key needs exactly one label.")
  offsets <- c(anc = -0.25, stk = 0, stp = 0.25)
  x <- x |>
    left_join(label_map, by = "protein_key", relationship = "many-to-one") |>
    mutate(label = coalesce(label, protein_key),
           # Include the key only if supplied labels would otherwise merge facets.
           panel_label = label,
           genotype = factor(genotype, levels = GENOTYPES),
           medium = factor(medium, levels = c("TSB", "SN")),
           x_position = as.numeric(medium) + unname(offsets[as.character(genotype)]))
  label_collision <- x |> distinct(protein_key, label) |> count(label) |> filter(n > 1) |> pull(label)
  x <- x |> mutate(panel_label = if_else(label %in% label_collision, paste0(label, "\n", protein_key), label))
  observed_counts <- x |> group_by(protein_key, genotype, medium) |>
    summarise(n = sum(is.finite(normalized_abundance)), .groups = "drop")
  counts <- expand_grid(protein_key = unique(x$protein_key), genotype = GENOTYPES, medium = c("TSB", "SN")) |>
    left_join(observed_counts |> mutate(genotype = as.character(genotype), medium = as.character(medium)),
              by = c("protein_key", "genotype", "medium")) |>
    left_join(x |> distinct(protein_key, panel_label), by = "protein_key") |>
    mutate(n = coalesce(n, 0L), text = paste0(n, "/3"),
           genotype = factor(genotype, levels = GENOTYPES),
           x_position = match(medium, c("TSB", "SN")) + unname(offsets[as.character(genotype)]))
  ggplot(x, aes(x_position, normalized_abundance, colour = genotype)) +
    geom_point(position = position_jitter(width = 0.035, height = 0, seed = 1), size = 1.9, alpha = 0.9, na.rm = TRUE) +
    geom_text(data = counts, aes(y = Inf, label = text), vjust = 1.3, size = 2.8, show.legend = FALSE) +
    facet_wrap(~panel_label, scales = "free_y", ncol = min(3L, n_distinct(x$protein_key))) +
    scale_colour_manual(values = GENOTYPE_COLOURS, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 2), labels = c("TSB", "SN"), limits = c(0.55, 2.45)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
    labs(x = "Medium", y = "Normalized abundance (log2, supplied vsn values)", colour = "Genotype",
         caption = "Each point is one sample measurement. Numbers show measured values out of 3; missing values remain NA.\nPanels have separate y-scales; points are not paired across media.") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0), strip.text = element_text(size = 9))
}
