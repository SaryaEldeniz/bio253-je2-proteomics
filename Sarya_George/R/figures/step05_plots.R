# Step 5 figures. These helpers display supplied enrichment results; they do not
# run enrichment tests, choose backgrounds or alter source p-values/FDRs.

.step05_empty_plot <- function(method, fdr_column, fdr_cutoff) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No terms pass the selected FDR threshold", size = 4.5) +
    labs(title = paste(method, "functional enrichment"),
         subtitle = paste0(fdr_column, " < ", fdr_cutoff),
         caption = "An empty selection does not establish absence of a biological response.") +
    theme_void() + theme(plot.caption = element_text(hjust = 0))
}

.step05_plot_input <- function(x, contrasts, policy, fdr_column, required) {
  required <- unique(c("term_description", "collection", "contrast", "contrast_label", fdr_column, required))
  if (!all(required %in% names(x))) stop("Enrichment plot input is missing: ", paste(setdiff(required, names(x)), collapse = ", "))
  if (!is.null(contrasts)) x <- filter(x, contrast %in% contrasts)
  if ("policy" %in% names(x)) {
    if (!is.null(policy)) x <- filter(x, .data$policy %in% .env$policy)
    if (n_distinct(x$policy) > 1L) stop("Choose one policy before combining enrichment rows in a plot.")
  }
  if ("namespace" %in% names(x) && n_distinct(x$namespace) > 1L) {
    stop("Choose one identifier namespace before plotting.")
  }
  term <- if ("term_id" %in% names(x)) x$term_id else if ("term" %in% names(x)) x$term else x$term_description
  x |>
    mutate(.q = .data[[fdr_column]], .term_key = paste(collection, .env$term, sep = "::"),
           .term_label = paste0("[", collection, "] ", term_description))
}

.step05_q_display <- function(q) {
  # Only the plotting transformation is floored; returned result tables keep q.
  -log10(pmax(q, .Machine$double.xmin))
}

.step05_direction_key <- function(x, width = 145L) {
  if (!all(c("positive_direction", "negative_direction") %in% names(x))) {
    return("NES direction follows the declared signed ranking; see the comparison-direction table.")
  }
  directions <- x |> distinct(contrast, contrast_label, positive_direction, negative_direction)
  paste(str_wrap(paste0(directions$contrast_label, ": NES > 0 = ", directions$positive_direction,
                        "; NES < 0 = ", directions$negative_direction, "."), width), collapse = "\n")
}

plot_step05_ora <- function(x, contrasts = NULL, policy = NULL, n_per_panel = 6L,
                             fdr_column = "FDR_study", fdr_cutoff = 0.1) {
  if (length(n_per_panel) != 1 || !is.finite(n_per_panel) || n_per_panel < 1) stop("n_per_panel must be positive.")
  x <- .step05_plot_input(x, contrasts, policy, fdr_column,
                          c("direction", "fold_enrichment", "overlap")) |>
    filter(is.finite(.q), .q >= 0, .q < fdr_cutoff,
           is.finite(fold_enrichment), fold_enrichment > 0, overlap > 0) |>
    group_by(contrast, direction) |>
    arrange(.q, desc(overlap), desc(fold_enrichment), .term_key, .by_group = TRUE) |>
    slice_head(n = as.integer(n_per_panel)) |> ungroup()
  if (!nrow(x)) return(.step05_empty_plot("ORA", fdr_column, fdr_cutoff))
  x <- x |>
    mutate(.panel = paste(str_wrap(contrast_label, 32), str_wrap(direction, 38), sep = "\n"),
           .row_key = paste(.panel, .term_key, sep = "___"),
           .logq = .step05_q_display(.q))
  row_levels <- rev(unique(x$.row_key))
  row_labels <- setNames(str_wrap(x$.term_label, 38), x$.row_key)
  x$.row_key <- factor(x$.row_key, levels = row_levels)
  ggplot(x, aes(fold_enrichment, .row_key, size = overlap, colour = .logq)) +
    geom_vline(xintercept = 1, colour = "grey65", linetype = "dashed", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    facet_wrap(~.panel, scales = "free_y", ncol = min(2L, n_distinct(x$.panel))) +
    scale_y_discrete(labels = row_labels) +
    scale_colour_viridis_c(option = "C", name = paste0("-log10(", fdr_column, ")")) +
    scale_size_continuous(range = c(2.5, 7), name = "Overlap\nproteins") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(title = "Overrepresented functional categories",
         subtitle = paste0("Top ", n_per_panel, " terms per contrast/direction passing ", fdr_column, " < ", fdr_cutoff),
         x = "Fold enrichment relative to the declared measured background", y = NULL,
         caption = "Point size counts overlapping identifiers. Related terms may share the same proteins. Enrichment does not measure pathway activity.") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9), axis.text.y = element_text(size = 9),
          plot.caption = element_text(hjust = 0), plot.margin = margin(8, 12, 8, 8))
}

plot_step05_gsea <- function(x, contrasts = NULL, policy = NULL, n_per_panel = 6L,
                              fdr_column = "FDR_study", fdr_cutoff = 0.1) {
  if (length(n_per_panel) != 1 || !is.finite(n_per_panel) || n_per_panel < 1) stop("n_per_panel must be positive.")
  x <- .step05_plot_input(x, contrasts, policy, fdr_column, c("NES", "size")) |>
    filter(is.finite(.q), .q >= 0, .q < fdr_cutoff, is.finite(NES), size > 0) |>
    group_by(contrast) |>
    arrange(.q, desc(abs(NES)), .term_key, .by_group = TRUE) |>
    slice_head(n = as.integer(n_per_panel)) |> ungroup()
  if (!nrow(x)) return(.step05_empty_plot("GSEA", fdr_column, fdr_cutoff))
  panel <- str_wrap(x$contrast_label, 35)
  if (all(c("positive_direction", "negative_direction") %in% names(x))) {
    panel <- paste(panel, str_wrap(paste0("NES > 0: ", x$positive_direction), 42),
                   str_wrap(paste0("NES < 0: ", x$negative_direction), 42), sep = "\n")
  }
  x <- x |>
    mutate(.panel = .env$panel, .row_key = paste(.panel, .term_key, sep = "___"),
           .logq = .step05_q_display(.q))
  row_levels <- rev(unique(x$.row_key))
  row_labels <- setNames(str_wrap(x$.term_label, 38), x$.row_key)
  x$.row_key <- factor(x$.row_key, levels = row_levels)
  ggplot(x, aes(NES, .row_key, size = size, colour = .logq)) +
    geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    facet_wrap(~.panel, scales = "free_y", ncol = min(2L, n_distinct(x$.panel))) +
    scale_y_discrete(labels = row_labels) +
    scale_colour_viridis_c(option = "C", name = paste0("-log10(", fdr_column, ")")) +
    scale_size_continuous(range = c(2.5, 7), name = "Set size\nin ranked list") +
    labs(title = "Ranked functional enrichment",
         subtitle = paste0("Top ", n_per_panel, " terms per contrast passing ", fdr_column, " < ", fdr_cutoff),
         x = "Normalized enrichment score (NES)", y = NULL,
         caption = "NES direction follows the declared signed ranking statistic. Enrichment does not measure pathway activity.") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9), axis.text.y = element_text(size = 9),
          plot.caption = element_text(hjust = 0), plot.margin = margin(8, 12, 8, 8))
}

plot_step05_nes_heatmap <- function(x, contrasts = NULL, policy = NULL, max_terms = 20L,
                                     fdr_column = "FDR_study", fdr_cutoff = 0.1) {
  if (length(max_terms) != 1 || !is.finite(max_terms) || max_terms < 1) stop("max_terms must be positive.")
  x <- .step05_plot_input(x, contrasts, policy, fdr_column, "NES")
  if (anyDuplicated(x[c("contrast", ".term_key")])) stop("Heatmap needs one row per term/contrast; choose one policy, namespace and collection representation.")
  eligible <- x |> filter(is.finite(.q), .q >= 0, .q < fdr_cutoff, is.finite(NES))
  if (!nrow(eligible)) return(.step05_empty_plot("GSEA", fdr_column, fdr_cutoff))
  selected <- eligible |>
    group_by(.term_key) |>
    summarise(.best_q = min(.q), .largest_abs = max(abs(NES)),
              .term_label = first(.term_label), .groups = "drop") |>
    arrange(.best_q, desc(.largest_abs), .term_key) |> slice_head(n = as.integer(max_terms))
  if (!nrow(selected)) return(.step05_empty_plot("GSEA", fdr_column, fdr_cutoff))
  wanted_contrasts <- if (is.null(contrasts)) unique(x$contrast) else contrasts
  labels <- x |> distinct(contrast, contrast_label)
  grid <- expand_grid(.term_key = selected$.term_key, contrast = wanted_contrasts) |>
    left_join(select(x, .term_key, contrast, NES, .q), by = c(".term_key", "contrast"), relationship = "one-to-one") |>
    left_join(labels, by = "contrast", relationship = "many-to-one") |>
    mutate(.term_key = factor(.term_key, levels = rev(selected$.term_key)),
           contrast = factor(contrast, levels = wanted_contrasts),
           .passes = is.finite(.q) & .q >= 0 & .q < fdr_cutoff & is.finite(NES))
  colour_limit <- max(abs(grid$NES[is.finite(grid$NES)]))
  if (!is.finite(colour_limit) || colour_limit == 0) colour_limit <- 1
  row_labels <- setNames(str_wrap(selected$.term_label, 42), selected$.term_key)
  column_labels <- setNames(str_wrap(labels$contrast_label, 24), labels$contrast)
  ggplot(grid, aes(contrast, .term_key, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = filter(grid, .passes), shape = 16, colour = "black", size = 1.5) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                          limits = c(-colour_limit, colour_limit), na.value = "grey85", name = "NES") +
    scale_y_discrete(labels = row_labels) + scale_x_discrete(labels = column_labels, drop = FALSE) +
    labs(title = "Functional enrichment across comparisons",
         subtitle = paste0("Terms passing ", fdr_column, " < ", fdr_cutoff, " in at least one comparison; tested nonsignificant cells are retained"),
         x = NULL, y = NULL,
         caption = paste0("Black dot: ", fdr_column, " < ", fdr_cutoff,
                          ". Grey: unavailable/untested NES. Enrichment does not measure pathway activity.\n",
                          .step05_direction_key(x))) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9), plot.caption = element_text(hjust = 0),
          plot.margin = margin(8, 12, 8, 8))
}
