# Poster SVG/PNG exports of the existing Step 4 plots. No models are rerun.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

export_step04_poster_figures <- function(figures = NULL) {
  step_dir <- file.path(REPO_ROOT, "steps", "04_total_proteome")
  figure_dir <- file.path(step_dir, "figures")
  report_dir <- file.path(step_dir, "reports")
  archive_dir <- file.path(REPO_ROOT, "90_archive", "figures", "step04_before_poster_styling")
  result_dir <- file.path(REPO_ROOT, "results", "step04", "poster_export")
  for (path in c(figure_dir, report_dir, archive_dir, result_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  plan <- step04_contrast_plan()
  map <- read_identifier_map()
  data_dir <- file.path(REPO_ROOT, "results", "step04", "total")
  all_results <- read_tsv(file.path(data_dir, "annotated_all.tsv.gz"), show_col_types = FALSE)
  summary <- read_tsv(file.path(data_dir, "hit_summary.tsv"), show_col_types = FALSE)
  abundances <- read_abundances("total", drop_contaminants = FALSE)
  total <- filter(all_results, eligible_cds)
  environmental <- filter(plan, contrast_type == "environmental")
  interactions <- filter(plan, contrast_type == "interaction")

  # Reconstruct the original plot objects with precisely the same helpers,
  # source estimates, feature selection, factor order and seeded sample jitter.
  original <- list()
  original$figure02_total_PASN_volcanoes <-
    plot_step04_volcano(total, environmental$contrast, environmental$contrast_label) +
    labs(title = "Total proteome: long PASN response within each genotype",
         x = "Source log2(TSB / SN): negative = higher in SN; positive = higher in TSB")
  original$figure03_total_interaction_volcanoes <-
    plot_step04_volcano(total, interactions$contrast, c("stk vs anc", "stp vs anc", "stp vs stk")) +
    labs(title = "Total proteome: genotype differences in PASN response",
         subtitle = "Define R = log2(SN / TSB). Panels show R_anc - R_stk, R_anc - R_stp, and R_stk - R_stp.",
         x = "Interaction effect (log2); inspect within-genotype responses to interpret its sign")
  original$figure04_total_hit_counts <- plot_step04_hit_counts(summary)

  candidates <- read_tsv(file.path(REPO_ROOT, "config", "interesting_proteins.tsv"), show_col_types = FALSE)
  lookup <- lookup_je2(candidates$query)
  candidate_keys <- unique(lookup$protein_key[lookup$lookup_status == "matched"])
  labels <- map |> transmute(protein_key, label = paste(display_name, locus_tag, sep = "\n"))
  original$figure05_candidate_measurements <- plot_step04_samples(abundances, candidate_keys, labels)
  examples <- total |> filter(contrast_type == "interaction", conservative_hit) |>
    group_by(contrast) |> arrange(FDR, desc(abs(diff)), feature_id, .by_group = TRUE) |>
    slice_head(n = 1) |> ungroup()
  if (nrow(examples)) {
    original$figure06_interaction_examples <-
      plot_step04_samples(abundances, unique(examples$protein_key), labels)
  }
  stopifnot(length(original) == 5L)
  titles <- c(
    figure02_total_PASN_volcanoes = "Total proteome: long PASN response within each genotype",
    figure03_total_interaction_volcanoes = "Total proteome: genotype differences in PASN response",
    figure04_total_hit_counts = "Total proteome: selected protein entries by contrast",
    figure05_candidate_measurements = "Total proteome: candidate protein measurements",
    figure06_interaction_examples = "Total proteome: examples of genotype-dependent PASN responses"
  )
  dimensions <- tibble(
    figure = names(titles),
    width = c(12, 12, 12.5, 10.5, 10.5),
    height = c(5.6, 5.6, 11, 9, 4.5)
  )

  sample_figures <- c("figure05_candidate_measurements", "figure06_interaction_examples")
  sample_means <- lapply(original[sample_figures], function(p) {
    p$data |> filter(is.finite(normalized_abundance)) |>
      group_by(protein_key, panel_label, genotype, medium, x_position) |>
      summarise(n_observed = n(), mean_log2_abundance = mean(normalized_abundance),
                .groups = "drop")
  })
  style_for_poster <- function(p, title, name) {
    styled <- p + labs(title = title) +
      # Keep the original 11-point base so theme-derived geom defaults (for
      # example bar-border widths) also remain exactly unchanged.
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 10, margin = margin(b = 8)),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 10.5, hjust = 0),
            legend.position = "bottom", legend.box = "horizontal",
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 10),
            legend.spacing.x = grid::unit(12, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            legend.key.width = grid::unit(14, "pt"),
            legend.key.height = grid::unit(14, "pt"),
            panel.spacing = grid::unit(12, "pt"),
            plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 6)),
            plot.margin = margin(7, 9, 7, 7)) +
      # A legend override has no effect on the points drawn in the plot panels.
      guides(colour = guide_legend(order = 1,
               override.aes = list(shape = 16, alpha = 1, size = 3.8)),
             shape = guide_legend(order = 2), fill = guide_legend(order = 1))
    if (name %in% sample_figures) {
      # Summarise the supplied log2 measurements, retaining missing groups as
      # absent. Short horizontal ticks show means separately in each medium;
      # they do not pair individual samples or connect independent groups.
      styled <- styled +
        geom_segment(data = sample_means[[name]],
          aes(x = x_position - 0.07, xend = x_position + 0.07,
              y = mean_log2_abundance, yend = mean_log2_abundance, colour = genotype),
          inherit.aes = FALSE, linewidth = 0.6, show.legend = FALSE) +
        labs(caption = paste0(
          "Points: individual measurements; horizontal bars: mean of observed log2 values. Numbers: measured values out of 3.\n",
          "Missing values remain NA. Panels have separate y-scales; points are not paired across media.")) +
        theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
              panel.grid.minor = element_blank())
    }
    styled
  }
  axis_signature <- function(build) {
    lapply(build$layout$panel_params, function(panel) {
      axis <- function(x) list(
        limits = x$limits, continuous_range = x$continuous_range,
        breaks = x$breaks, minor_breaks = x$minor_breaks,
        labels = x$get_labels(),
        transformation = x$scale$get_transformation()$name
      )
      list(x_range = panel$x.range, y_range = panel$y.range,
           x = axis(panel$x), y = axis(panel$y))
    })
  }
  checks <- list()
  archive_rows <- list()
  if (is.null(figures)) figures <- names(original)
  if (!length(figures) || !all(figures %in% names(original))) stop("Unknown Step 4 figure requested.")
  for (name in figures) {
    baseline <- original[[name]]
    poster <- style_for_poster(baseline, unname(titles[name]), name)
    before <- ggplot_build(baseline)
    after <- ggplot_build(poster)
    checks[[name]] <- tibble(
      figure = name,
      # Mean bars are an explicitly requested additional layer. All pre-existing
      # point/text/bar layers must remain identical, including seeded jitter.
      layer_data_unchanged = isTRUE(all.equal(before$data,
        after$data[seq_along(before$data)], check.attributes = TRUE)),
      facet_order_unchanged = identical(before$layout$layout, after$layout$layout),
      axes_scales_breaks_labels_unchanged = identical(axis_signature(before), axis_signature(after)),
      x_axis_title_unchanged = identical(baseline$labels$x, poster$labels$x),
      y_axis_title_unchanged = identical(baseline$labels$y, poster$labels$y)
    )
    if (!all(unlist(checks[[name]][, -1]))) {
      print(checks[[name]], width = Inf)
      print(all.equal(before$data, after$data, check.attributes = TRUE))
      stop("Poster styling changed plot data or scales: ", name)
    }
    if (name %in% sample_figures) {
      mean_layer <- after$data[[length(after$data)]]
      stopifnot(length(after$data) == length(before$data) + 1L,
                nrow(mean_layer) == nrow(sample_means[[name]]),
                isTRUE(all.equal(mean_layer$y, sample_means[[name]]$mean_log2_abundance)),
                identical(mean_layer$y, mean_layer$yend))
    }
    png_path <- file.path(figure_dir, paste0(name, ".png"))
    archive_path <- file.path(archive_dir, paste0(name, ".png"))
    if (file.exists(png_path) && !file.exists(archive_path)) {
      stopifnot(file.copy(png_path, archive_path, overwrite = FALSE),
                unname(tools::md5sum(png_path)) == unname(tools::md5sum(archive_path)))
    }
    archive_rows[[name]] <- tibble(figure = name,
      original_png = file.path("90_archive", "figures", "step04_before_poster_styling", paste0(name, ".png")),
      original_png_md5 = unname(tools::md5sum(archive_path)))
    size <- filter(dimensions, figure == name)
    ggsave(file.path(figure_dir, paste0(name, ".svg")), poster,
           width = size$width, height = size$height, units = "in", device = grDevices::svg, bg = "white")
    ggsave(png_path, poster, width = size$width, height = size$height,
           units = "in", dpi = 300, bg = "white")
    # SVGs contain paths/text/vector symbols, never an embedded PNG screenshot.
    svg <- readLines(file.path(figure_dir, paste0(name, ".svg")), warn = FALSE)
    stopifnot(any(grepl("<svg", svg, fixed = TRUE)), !any(grepl("<image", svg, fixed = TRUE)))
    saveRDS(list(original = list(data = before$data, axes = axis_signature(before)),
                 poster = list(data = after$data, axes = axis_signature(after))),
            file.path(result_dir, paste0(name, "_invariance.rds")))
  }
  checks <- bind_rows(checks)
  # A partial rerun updates only the requested files and retains other reports.
  update_report <- function(rows, path) {
    if (file.exists(path)) {
      rows <- bind_rows(read_tsv(path, show_col_types = FALSE) |>
                          filter(!figure %in% figures), rows)
    }
    write_tsv(arrange(rows, figure), path, quote = "needed")
  }
  update_report(checks, file.path(report_dir, "poster_style_invariance_checks.tsv"))
  update_report(bind_rows(archive_rows), file.path(report_dir, "poster_original_png_archive.tsv"))
  # Facet labels contain newlines: quote them so this summary stays a valid TSV.
  # Both sample summaries are available, so write a fresh table each time.
  write_tsv(bind_rows(sample_means, .id = "figure"),
            file.path(report_dir, "poster_sample_means.tsv"), quote = "needed")
  svg_paths <- file.path(figure_dir, paste0(names(original), ".svg"))
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(file.path(figure_dir, "step04_poster_SVGs.zip"), files = svg_paths, root = figure_dir)
  }
  message(length(figures), " Step 4 SVGs and matching PNG previews exported. All data/axis checks passed.")
  print(checks, width = Inf)
  invisible(checks)
}

step04_poster_checks <- export_step04_poster_figures(
  getOption("step04.poster.figures", NULL))
