# Recreate the Step 7-9 phosphosite comparisons in our earlier volcano style.
# Run from Sarya_George or the repository root. Source effects/FDRs are not refit.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

recreate_phospho_volcanoes <- function() {
  fdr_cutoff <- 0.10
  effect_cutoff <- 1
  localization_cutoff <- 0.75
  plan <- step04_contrast_plan() |>
    filter(contrast_type %in% c("interaction", "mutant"))
  interactions <- filter(plan, contrast_type == "interaction")$contrast
  mutants <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB",
               "stk_vs_anc_at_SN", "stp_vs_anc_at_SN")
  contrasts <- c(interactions, mutants)
  panel_titles <- setNames(c("stk vs ancestor", "stp vs ancestor", "stp vs stk",
                            "TSB: stk vs ancestor", "TSB: stp vs ancestor",
                            "SN: stk vs ancestor", "SN: stp vs ancestor"), contrasts)
  # Match the requested residue, not just a protein or a nearby phosphosite.
  arrow_targets <- tibble::tribble(
    ~site_id, ~arrow_label,
    "lcl_CP020619.1_prot_ARG47142.1_2776_S194", "IsaA S194",
    "lcl_CP020619.1_prot_ARG46254.1_1820_T537", "Pyk T537",
    "lcl_CP020619.1_prot_ARG45539.1_1070_T552", "QoxB T552"
  )
  figure_dir <- file.path(REPO_ROOT, "steps/07_phosphosite_results/figures/volcanoes")
  supplementary_figure_dir <- file.path(REPO_ROOT, "steps/08_mutant_venn/figures/volcanoes")
  report_dir <- file.path(REPO_ROOT, "steps/07_phosphosite_results/reports/volcano_reproduction")
  result_dir <- file.path(REPO_ROOT, "results", "step07_volcano_reproduction")
  for (path in c(figure_dir, supplementary_figure_dir, report_dir, result_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  # Keep the same population as the existing Step 7-9 scripts: no extra CDS or
  # UniProt filter. A point represents a phosphosite estimate, not a protein.
  source_dea <- read_dea("phospho_STYH", contrasts, drop_contaminants = FALSE)
  excluded <- source_dea |>
    filter(is_contaminant | !is.finite(diff) | !is.finite(FDR) | FDR < 0 | FDR > 1)
  dea <- source_dea |>
    filter(!is_contaminant, is.finite(diff), is.finite(FDR), FDR >= 0, FDR <= 1)
  stopifnot(!anyDuplicated(dea[c("site_id", "contrast")]))

  # Confidence and replicate counts are annotations; they do not silently change
  # the source-reproduction plots or the multiple-testing correction.
  site_info <- read_sites("STYH", drop_contaminants = FALSE) |>
    distinct(site_id, best_loc_prob)
  stopifnot(!anyDuplicated(site_info$site_id))
  abundances <- read_abundances("phospho_STYH", drop_contaminants = FALSE)
  stopifnot(!anyDuplicated(abundances[c("site_id", "sample")]))
  support <- step04_group_support(abundances, plan, "site_id")
  total <- read_dea("total", contrasts) |>
    transmute(protein_key, contrast, parent_diff = diff, parent_FDR = FDR,
              parent_estimate_type = estimate_type)
  stopifnot(!anyDuplicated(total[c("protein_key", "contrast")]))
  id_map <- read_identifier_map() |>
    select(protein_key, display_name, uniprot_acc, entry_class)
  stopifnot(!anyDuplicated(id_map$protein_key))

  annotated <- dea |>
    left_join(site_info, by = "site_id", relationship = "many-to-one") |>
    left_join(support$support, by = c("site_id" = "feature_id", "contrast"),
              relationship = "one-to-one") |>
    left_join(total, by = c("protein_key", "contrast"), relationship = "many-to-one") |>
    left_join(id_map, by = "protein_key", relationship = "many-to-one") |>
    left_join(select(plan, contrast, positive_direction, negative_direction, formula),
              by = "contrast", relationship = "many-to-one") |>
    mutate(selected_hit = FDR < fdr_cutoff & abs(diff) > effect_cutoff,
           source_observed = estimate_type == "observed",
           localization_pass = coalesce(best_loc_prob >= localization_cutoff, FALSE),
           eligible_je2 = coalesce(entry_class == "cds", FALSE),
           conservative_quality_pass = eligible_je2 & source_observed & localization_pass &
             coalesce(all_relevant_groups_2of3, FALSE),
           conservative_hit = selected_hit & conservative_quality_pass,
           parent_testable = is.finite(parent_diff) & is.finite(parent_FDR),
           parent_selected = parent_testable & parent_FDR < fdr_cutoff &
             abs(parent_diff) > effect_cutoff,
           parent_relationship = case_when(
             !parent_testable ~ "Parent estimate unavailable",
             !parent_selected ~ "Parent does not meet both thresholds",
             sign(diff) == sign(parent_diff) ~ "Parent selected, same sign",
             TRUE ~ "Parent selected, opposite sign"),
           direction = if_else(diff >= 0, positive_direction, negative_direction),
           site_label = paste0(coalesce(display_name, locus_tag), " ", residue, position))
  stopifnot(nrow(annotated) == nrow(dea), identical(annotated$diff, dea$diff),
            identical(annotated$FDR, dea$FDR), !anyNA(annotated$best_loc_prob),
            !anyNA(annotated$n_relevant_groups),
            all(support$groups$n_samples == 3L))

  counts <- annotated |>
    group_by(contrast) |>
    summarise(sites_plotted_all = n(), sites_plotted_observed = sum(source_observed),
              hits_all = sum(selected_hit), hits_observed = sum(selected_hit & source_observed),
              hits_lod = sum(selected_hit & estimate_type == "lod_imputed"),
              sites_plotted_je2_observed = sum(eligible_je2 & source_observed),
              hits_je2_observed = sum(selected_hit & eligible_je2 & source_observed),
              all_positive = sum(selected_hit & diff > 0),
              all_negative = sum(selected_hit & diff < 0),
              observed_positive = sum(selected_hit & source_observed & diff > 0),
              observed_negative = sum(selected_hit & source_observed & diff < 0),
              observed_hits_localization_pass = sum(selected_hit & source_observed & localization_pass),
              conservative_hits = sum(conservative_hit), .groups = "drop") |>
    arrange(match(contrast, contrasts))
  stopifnot(all(counts$hits_all == counts$all_positive + counts$all_negative),
            all(counts$hits_all == counts$hits_observed + counts$hits_lod))
  # Compare independent sums to the saved Step 7/8 report counts.
  old_interactions <- read_tsv(file.path(REPO_ROOT, "steps/07_phosphosite_results/reports/step07_phospho_parent_counts.tsv"),
                               show_col_types = FALSE)
  old_venn <- read_tsv(file.path(REPO_ROOT, "steps/08_mutant_venn/reports/step08_phospho_venn_counts.tsv"),
                      show_col_types = FALSE)
  stopifnot(identical(as.integer(counts$hits_all[1:3]),
                      as.integer(old_interactions$significant_sites)),
            identical(as.integer(counts$hits_all[4:7]),
                      as.integer(c(old_venn$stk_only[1] + old_venn$both[1],
                                   old_venn$stp_only[1] + old_venn$both[1],
                                   old_venn$stk_only[2] + old_venn$both[2],
                                   old_venn$stp_only[2] + old_venn$both[2]))))

  # All point coordinates and per-sample measurements stay in ignored results/.
  write_tsv(annotated, file.path(result_dir, "all_annotated_coordinates.tsv.gz"))
  write_tsv(excluded, file.path(result_dir, "excluded_source_rows.tsv.gz"))
  write_tsv(support$groups, file.path(result_dir, "site_group_detection_counts.tsv.gz"))
  selected_ids <- unique(annotated$site_id[annotated$selected_hit])
  write_tsv(filter(abundances, site_id %in% selected_ids, medium %in% c("TSB", "SN")),
            file.path(result_dir, "selected_site_sample_values.tsv.gz"))
  write_tsv(counts, file.path(report_dir, "volcano_counts.tsv"))
  hit_columns <- c("contrast", "site_id", "site_label", "protein_key", "locus_tag",
                   "display_name", "uniprot_acc", "residue", "position", "diff", "FDR",
                   "direction", "estimate_type", "entry_class", "eligible_je2", "best_loc_prob", "localization_pass",
                   "n_observed_min_relevant", "all_relevant_groups_2of3",
                   "group_detection_pattern", "conservative_hit", "parent_diff", "parent_FDR",
                   "parent_estimate_type", "parent_relationship")
  write_tsv(annotated |> filter(selected_hit) |> select(all_of(hit_columns)) |>
              arrange(match(contrast, contrasts), FDR), file.path(report_dir, "selected_sites_annotated.tsv"))
  arrow_audit <- annotated |>
    inner_join(arrow_targets, by = "site_id", relationship = "many-to-one") |>
    select(contrast, site_id, arrow_label, diff, FDR, estimate_type, selected_hit)
  stopifnot(nrow(arrow_audit) == 21L, all(arrow_audit$estimate_type == "observed"),
            !any(arrow_audit$selected_hit[arrow_audit$contrast %in% interactions]),
            sum(arrow_audit$selected_hit) == 3L)
  write_tsv(filter(arrow_audit, contrast %in% interactions),
            file.path(report_dir, "labelled_interaction_sites.tsv"))
  write_tsv(filter(arrow_audit, contrast %in% mutants),
            file.path(report_dir, "labelled_mutant_sites.tsv"))

  # Parametric segment/rectangle intersection, vectorised over label candidates.
  # This prevents one label from covering another site's leader line.
  segment_crosses_label <- function(x0, y0, x1, y1, cx, cy, hw, hh) {
    axis_interval <- function(start, end, centre, half_size) {
      # Broadcast a fixed segment over multiple candidate rectangle centres.
      delta <- end - start + 0 * centre
      parallel <- abs(delta) < 1e-12
      inside <- abs(start - centre) <= half_size
      low <- (centre - half_size - start) / delta
      high <- (centre + half_size - start) / delta
      list(enter = ifelse(parallel, ifelse(inside, -Inf, Inf), pmin(low, high)),
           leave = ifelse(parallel, ifelse(inside, Inf, -Inf), pmax(low, high)))
    }
    x <- axis_interval(x0, x1, cx, hw)
    y <- axis_interval(y0, y1, cy, hh)
    pmax(0, x$enter, y$enter) <= pmin(1, x$leave, y$leave)
  }
  stopifnot(identical(segment_crosses_label(-2, 0, 2, 0, c(0, 0, 4), c(0, 3, 0), .5, .5),
                      c(TRUE, FALSE, FALSE)),
            identical(segment_crosses_label(0, -2, 0, 2, c(0, 3, 0), c(0, 0, 4), .5, .5),
                      c(TRUE, FALSE, FALSE)))

  # Find nearby empty regions for three compact, two-line labels. The padded
  # rectangles are deliberately larger than the rendered boxes. No plotted
  # point may lie in one, and labels may not overlap one another.
  place_site_labels <- function(points, targets, x_limit, y_limit) {
    x_span <- 2 * x_limit
    half_width <- .065 * x_span
    half_height <- .065 * y_limit
    candidates <- expand_grid(
      label_x = seq(-x_limit + half_width, x_limit - half_width, length.out = 61),
      label_y = seq(half_height + .02 * y_limit,
                    y_limit - half_height, length.out = 51))
    candidates$n_points_under_label <- vapply(seq_len(nrow(candidates)), function(i) {
      sum(abs(points$diff - candidates$label_x[i]) <= half_width + .008 * x_span &
            abs(points$y - candidates$label_y[i]) <= half_height + .01 * y_limit)
    }, integer(1))
    candidates <- filter(candidates, n_points_under_label == 0L)
    if (!nrow(candidates)) stop("No empty region is available for site labels.")
    targets <- arrange(targets, desc(y), diff)
    placed <- list()
    for (j in seq_len(nrow(targets))) {
      options <- candidates
      for (previous in placed) {
        options <- filter(options,
          abs(label_x - previous$label_x) > 2 * half_width + .018 * x_span |
            abs(label_y - previous$label_y) > 2 * half_height + .018 * y_limit)
        new_arrow_hits_old_label <- segment_crosses_label(
          options$label_x, options$label_y, targets$diff[j], targets$y[j],
          previous$label_x, previous$label_y, half_width, half_height)
        old_arrow_hits_new_label <- segment_crosses_label(
          previous$label_x, previous$label_y, previous$diff, previous$y,
          options$label_x, options$label_y, half_width, half_height)
        options <- options[!new_arrow_hits_old_label & !old_arrow_hits_new_label, ]
      }
      if (!nrow(options)) stop("No nonoverlapping site-label placement is available.")
      distance <- ((options$label_x - targets$diff[j]) / x_span)^2 +
        ((options$label_y - targets$y[j]) / y_limit)^2
      chosen <- options[which.min(distance), ]
      placed[[j]] <- bind_cols(targets[j, ], chosen) |>
        mutate(label_half_width = half_width, label_half_height = half_height,
               arrow_label_plot = sub(" ", "\n", arrow_label, fixed = TRUE))
    }
    labels <- bind_rows(placed)
    # Begin each arrow at the label edge rather than drawing through its text.
    vx <- labels$diff - labels$label_x
    vy <- labels$y - labels$label_y
    edge <- pmin(.80 * half_width / abs(vx), .85 * half_height / abs(vy))
    labels$arrow_x <- labels$label_x + edge * vx
    labels$arrow_y <- labels$label_y + edge * vy
    stopifnot(all(labels$n_points_under_label == 0L))
    for (a in seq_len(nrow(labels))) {
      others <- labels[-a, ]
      stopifnot(!any(segment_crosses_label(labels$arrow_x[a], labels$arrow_y[a],
        labels$diff[a], labels$y[a], others$label_x, others$label_y,
        half_width, half_height)))
    }
    labels
  }
  label_layouts <- list()

  # Exact unbinned coordinates; classic theme and blue/orange colours match
  # steps/10_integration_poster/figures/slide_summary/volcano_1.png and the earlier total-proteome trio.
  make_plot <- function(which_contrasts, observed_only = FALSE, ncol = 3, individual = FALSE,
                        je2_only = FALSE) {
    d <- filter(annotated, contrast %in% which_contrasts)
    if (observed_only) d <- filter(d, source_observed)
    if (je2_only) d <- filter(d, eligible_je2)
    positive_fdr <- d$FDR[d$FDR > 0]
    floor_fdr <- if (length(positive_fdr)) max(.Machine$double.xmin, min(positive_fdr) / 10) else .Machine$double.xmin
    d <- d |>
      mutate(y = -log10(pmax(FDR, floor_fdr)),
             effect_class = factor(case_when(selected_hit & diff < 0 ~ "Negative effect",
                                             selected_hit & diff > 0 ~ "Positive effect",
                                             TRUE ~ "Below thresholds"),
                                   levels = c("Below thresholds", "Negative effect", "Positive effect")),
             source_model = factor(if_else(source_observed, "Source observed", "LOD-imputed"),
                                   levels = c("Source observed", "LOD-imputed")))
    panels <- d |> group_by(contrast) |>
      summarise(positive = sum(selected_hit & diff > 0), negative = sum(selected_hit & diff < 0),
                n_sites = n(), observed_hits = sum(selected_hit & source_observed), .groups = "drop")
    panel_labels <- setNames(paste0(panel_titles[panels$contrast], "\nPositive ", panels$positive,
                                    "; negative ", panels$negative), panels$contrast)
    d$panel <- factor(d$contrast, levels = which_contrasts,
                      labels = unname(panel_labels[which_contrasts]))
    # Draw selected and observed points last so they remain visible among crosses.
    d <- arrange(d, selected_hit, source_observed)
    x_limit <- max(2, ceiling(max(abs(d$diff))))
    y_limit <- max(2, ceiling(max(d$y)))
    is_interaction <- all(which_contrasts %in% interactions)
    title_text <- if (is_interaction) "Phosphosites: genotype differences in PASN response" else {
      "Phosphosites: mutant versus ancestor in TSB and SN"
    }
    subtitle_text <- if (is_interaction) {
      "PASN response = log\u2082(SN/TSB); interaction = response(second genotype) \u2212 response(first genotype)."
    } else "Within each medium: site abundance difference = log\u2082(mutant/ancestor)."
    p <- ggplot(d, aes(diff, y, colour = effect_class, shape = source_model)) +
      geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), linetype = "dashed",
                 colour = "#8A97A6", linewidth = .35) +
      geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed",
                 colour = "#8A97A6", linewidth = .35) +
      geom_point(size = .95, alpha = .75, stroke = .45) +
      scale_colour_manual(values = c("Below thresholds" = "#B4BEC9", "Negative effect" = "#0072B2",
                                     "Positive effect" = "#D55E00"), drop = FALSE) +
      scale_shape_manual(values = c("Source observed" = 16, "LOD-imputed" = 4), drop = TRUE) +
      scale_x_continuous(limits = c(-x_limit, x_limit), breaks = pretty(c(-x_limit, x_limit), 5)) +
      scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, .03))) +
      labs(title = title_text, subtitle = subtitle_text,
           x = if (is_interaction) "Interaction effect (log2)" else {
             "log2(mutant / ancestor): negative = higher in ancestor; positive = higher in mutant"
           }, y = "-log10(phosphosite FDR)", colour = NULL, shape = NULL,
           caption = "Selected: FDR < 0.10 and |effect| > 1. Site effects are not adjusted for parent-protein abundance.") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 10, margin = margin(b = 8)),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 10.5, hjust = 0),
            legend.position = "bottom", legend.box = "horizontal",
            legend.spacing.x = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 9), legend.key.width = grid::unit(12, "pt"),
            panel.spacing = grid::unit(12, "pt"),
            plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 6)),
            plot.margin = margin(7, 9, 7, 7)) +
      guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 1,
                                                               size = 3.8)),
             shape = if (observed_only) "none" else guide_legend(order = 2, override.aes = list(colour = "#526273", alpha = 1)))
    targets <- inner_join(d, arrow_targets, by = "site_id", relationship = "many-to-one")
    stopifnot(nrow(targets) == 3L * length(which_contrasts))
    targets <- bind_rows(lapply(which_contrasts, function(contrast_name) {
      place_site_labels(filter(d, contrast == contrast_name),
                        filter(targets, contrast == contrast_name), x_limit, y_limit)
    }))
    # Arrows terminate at unchanged source coordinates in both figure sets.
      p <- p +
        geom_segment(data = targets,
          aes(x = arrow_x, y = arrow_y, xend = diff, yend = y),
          inherit.aes = FALSE, colour = "#34495B", linewidth = .30,
          arrow = grid::arrow(length = grid::unit(1.35, "mm"), angle = 22, type = "open"),
          show.legend = FALSE) +
        geom_label(data = targets, aes(x = label_x, y = label_y, label = arrow_label_plot),
          inherit.aes = FALSE, colour = "#263A49", fill = "white", size = 3.0,
          linewidth = 0, lineheight = .9, label.padding = grid::unit(.10, "lines"),
          label.r = grid::unit(0, "lines"), show.legend = FALSE)
    if (individual) {
      p <- p + labs(title = unname(panel_titles[which_contrasts]),
                    subtitle = paste0("Positive ", panels$positive, "; negative ", panels$negative,
                                      " | ", panels$n_sites, " sites"),
                    caption = stringr::str_wrap(
                      "Selected: FDR < 0.10 and |effect| > 1. Site effects are not adjusted for parent-protein abundance.",
                      width = 80)) +
        theme(axis.title.x = element_text(size = 10), legend.text = element_text(size = 8))
    } else p <- p + facet_wrap(~panel, ncol = ncol)
    attr(p, "site_label_placement") <- targets
    p
  }
  save_plot <- function(p, name, width, height) {
    plot_dir <- if (grepl("^(phospho_mutants_|observed_(stk|stp)_vs_anc_at_)", name)) {
      supplementary_figure_dir
    } else figure_dir
    label_layouts[[name]] <<- attr(p, "site_label_placement") |>
      mutate(figure = name, .before = 1) |>
      select(figure, contrast, site_id, arrow_label, diff, y, label_x, label_y,
             label_half_width, label_half_height, n_points_under_label)
    ggsave(file.path(plot_dir, paste0(name, ".png")), p, width = width, height = height,
            dpi = 300, bg = "white")
    ggsave(file.path(plot_dir, paste0(name, ".pdf")), p, width = width, height = height,
            device = grDevices::cairo_pdf, bg = "white")
    ggsave(file.path(plot_dir, paste0(name, ".svg")), p, width = width, height = height,
            device = grDevices::svg, bg = "white")
  }
  for (observed_only in c(FALSE, TRUE)) {
    scope_name <- if (observed_only) "observed_only" else "all_models"
    save_plot(make_plot(interactions, observed_only, ncol = 3),
               paste0("phospho_interactions_", scope_name), 12.2, 4.9)
    save_plot(make_plot(mutants, observed_only, ncol = 2),
               paste0("phospho_mutants_TSB_SN_", scope_name), 10.5, 7.4)
  }
  save_plot(make_plot(interactions, observed_only = TRUE, ncol = 3, je2_only = TRUE),
               "phospho_interactions_JE2_observed_only", 12.2, 4.9)
  save_plot(make_plot(mutants, observed_only = TRUE, ncol = 2, je2_only = TRUE),
               "phospho_mutants_TSB_SN_JE2_observed_only", 10.5, 7.4)
  # Separate JE2 observed-only panels are ready to arrange on a poster.
  for (contrast in contrasts) {
    save_plot(make_plot(contrast, observed_only = TRUE, individual = TRUE, je2_only = TRUE),
               paste0("observed_", contrast), 5.5, 4.3)
  }
  write_tsv(bind_rows(label_layouts), file.path(report_dir, "site_label_placement_checks.tsv"))

  paths <- c(derived_path("dea_phospho_STYH.tsv.gz"), derived_path("sites_STYH.tsv.gz"),
             derived_path("abundances_phospho_STYH.tsv.gz"), derived_path("dea_total.tsv.gz"))
  write_tsv(tibble(input = basename(paths), md5 = unname(tools::md5sum(paths))),
             file.path(report_dir, "input_checksums.tsv"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  markdown_counts <- c("| Contrast | All source hits | Observed source hits | LOD-imputed hits | JE2 observed hits | Quality subset hits |",
                       "| --- | ---: | ---: | ---: | ---: | ---: |",
                       vapply(seq_len(nrow(counts)), function(i) {
                         z <- counts[i, ]
                         sprintf("| `%s` | %d | %d | %d | %d | %d |", z$contrast, z$hits_all,
                                  z$hits_observed, z$hits_lod, z$hits_je2_observed, z$conservative_hits)
                       }, character(1)))
  writeLines(c("# Phosphosite volcano reproduction", "",
    "The October STYH DEA is redrawn in the compact classic-theme blue/orange style of the earlier total-proteome volcanoes. Each point is one phosphosite estimate, not one protein. Effects and FDRs are copied unchanged; no new tests or imputation were performed.", "",
    "## Comparisons and selection", "",
    "The three Step 7 panels test genotype-by-medium interactions. The four Step 8/9 panels compare each mutant with the ancestor separately in TSB and SN. SN is long PASN exposure, not the 15-minute pulse.", "",
    "Selection is FDR < 0.10 and absolute log2 effect > 1. Source-contaminant and nonfinite entries are excluded. No new UniProt/CDS eligibility filter is imposed; it would change the population used in the original Step 7-9 figures.", "",
    "All-model panels reproduce the existing Step 7/8 hit counts. Filled circles are source-observed estimates; crosses are LOD-imputed estimates. Observed-only panels exclude all LOD-flagged rows. The observed label is global to the source model, not a promise of three detections in every contrast group.", "",
    "The separately named JE2_observed_only panels additionally retain only entries classified as JE2 CDS in the existing identifier map. This removes external entries that escaped the source contaminant flag. For example, P02662_S61 is bovine alpha-S1-casein: it is an observed selected stp-vs-stk interaction in the source reproduction, but is excluded from the JE2-specific plot and quality subset. Consequently that panel has three JE2 observed hits rather than four source-observed hits. The original all-model and observed-only reproductions remain available for comparison.", "",
    markdown_counts, "",
    "## Localization, replication and parent protein", "",
    "The annotated selected-site table supplies the best localization probability over all samples, the normalized-abundance detection count in every relevant group, and the matched total-protein effect/FDR. The illustrative quality subset requires JE2 CDS, source-observed, best localization probability >= 0.75 and >=2/3 measured values in each relevant group (four groups for an interaction, two for a simple contrast). It is an additional annotation, not a silent filter on the reproduction plots. Best localization probability does not guarantee high localization confidence in every sample.", "",
    "Glutamate dehydrogenase T269 (B7H15_04940, Q2FZU5) is a source-observed selected stp-vs-stk interaction, but only 1/3 kinase-mutant SN values were measured. Its best localization probability is 0.8957. It remains in the JE2 observed plot and is flagged outside the conservative quality subset. A large interaction estimate alone is not sufficient evidence of a well-replicated phosphorylation response.", "",
    "Site and parent-protein estimates were fitted separately. A nonsignificant parent effect does not establish unchanged protein abundance. These plots do not measure phosphorylation occupancy or establish direct kinase substrates. A protein-adjusted significance test has not been performed.", "",
    "Per-sample values for every selected site are saved locally in results/step07_volcano_reproduction/selected_site_sample_values.tsv.gz for inspection. No raw measurement tables are added to reports/.", "",
    "## Reading the sign", "",
    "In mutant-versus-ancestor panels, positive means higher site abundance in the mutant; negative means higher in the ancestor. For an interaction A vs B at TSB vs SN, the effect is (A_TSB - A_SN) - (B_TSB - B_SN). Defining R = SN - TSB, a positive effect means R_B > R_A; a negative effect means R_A > R_B. This does not by itself establish an increase or decrease within either genotype.", "",
    "All panels within a figure share symmetric x limits and a common y scale. Observed-only and all-model figures each use their own data-derived limits. There is no binning, jitter, coordinate clipping or FDR recomputation. Any source FDR of zero receives a finite coordinate only for display.", "",
    "## Outputs and rerun", "",
    "Both the three main interaction panels and the four supplementary mutant-versus-ancestor panels label IsaA S194, Pyk T537 and QoxB T552 with small arrows terminating at the original coordinates. All nine labelled interaction tests are nonsignificant; three of the twelve labelled simple comparisons meet selection. See labelled_interaction_sites.tsv and labelled_mutant_sites.tsv. Labels are placed in nearby empty regions; conservative padded rectangles are checked against every data point, the other labels and the other arrow paths. These placement checks are saved in site_label_placement_checks.tsv. Only the coloured legend markers are enlarged; plotted points retain their original size. Both figure sets use a single explanatory subtitle and omit the technical source-filter footer. The source/JE2 selection policies remain documented above and in the file names.", "",
    "The figures are in ../../figures/volcanoes/ and ../../../08_mutant_venn/figures/volcanoes/ as PNG, PDF and SVG. The individual panels use JE2 source-observed estimates and are available for poster assembly. See volcano_counts.tsv, selected_sites_annotated.tsv, input_checksums.tsv and session_info.txt beside this file.", "",
    "From Sarya_George.Rproj run:", "", "```r",
    'source("steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R")', "```", "",
    "From the repository root prefix the path with Sarya_George/. The script reuses R/setup.R and the existing site/DEA readers; course data stay in the configured local data directory."),
    file.path(report_dir, "README.md"))
  print(counts, width = Inf)
  message("Phosphosite volcanoes saved in: ", figure_dir)
  invisible(list(counts = counts, annotated = annotated))
}

phospho_volcano_results <- recreate_phospho_volcanoes()
