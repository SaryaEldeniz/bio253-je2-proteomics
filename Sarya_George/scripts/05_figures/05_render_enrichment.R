# Render Step 5 enrichment figures from saved tables; no enrichment is rerun.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

render_step05_enrichment <- function() {
  result_dir <- file.path(REPO_ROOT, "results", "step05")
  paths <- file.path(result_dir, c("ora_all.tsv.gz", "gsea_all.tsv.gz"))
  if (any(!file.exists(paths))) stop("Run the Step 5 enrichment analysis before rendering its figures.")
  ora <- read_tsv(paths[1], show_col_types = FALSE)
  gsea <- read_tsv(paths[2], show_col_types = FALSE)
  required <- c("policy", "namespace", "contrast", "contrast_label", "pathway", "collection",
                "term_id", "term_description", "direction", "FDR_study", "FDR_local")
  for (x in list(ora, gsea)) {
    if (!all(required %in% names(x))) stop("Saved enrichment table is missing: ", paste(setdiff(required, names(x)), collapse = ", "))
  }
  # Results in other policies/namespaces are sensitivity analyses, not pooled rows.
  ora <- filter(ora, policy == "primary", namespace == "UniProt")
  gsea <- filter(gsea, policy == "primary", namespace == "UniProt")
  if (!nrow(ora) && !nrow(gsea)) stop("No primary UniProt results were found; check the saved policy/namespace fields.")
  figure_dir <- file.path(REPO_ROOT, "figures", "step05")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ancestor <- "TSB_vs_SN_at_anc"
  interactions <- c("interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
                    "interaction_stp_vs_stk_at_TSB_vs_SN")
  focus <- c(ancestor, interactions)
  all_contrasts <- step04_contrast_plan()$contrast
  manifest <- list()
  save_pair <- function(plot, stem, width, height, method, collection, scope, candidates) {
    png <- file.path(figure_dir, paste0(stem, ".png"))
    pdf <- file.path(figure_dir, paste0(stem, ".pdf"))
    ggsave(png, plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
    ggsave(pdf, plot, width = width, height = height, units = "in", device = grDevices::pdf,
           useDingbats = FALSE, bg = "white")
    manifest[[length(manifest) + 1L]] <<- tibble(
      figure = stem, method, collection, scope, policy = "primary", namespace = "UniProt",
      selection = "FDR_study < 0.1; display caps stated in figure", significant_candidates = candidates,
      png = paste0("figures/step05/", stem, ".png"), pdf = paste0("figures/step05/", stem, ".pdf"))
  }
  passing <- function(x) filter(x, is.finite(FDR_study), FDR_study >= 0, FDR_study < 0.1)
  ora_height <- function(x, n = 6L) {
    panels <- passing(x) |> distinct(contrast, direction) |> nrow()
    if (!panels) return(4.5)
    max(5, ceiling(panels / 2) * (1.9 + 0.55 * n))
  }
  gsea_height <- function(x, n = 6L) {
    panels <- passing(x) |> distinct(contrast) |> nrow()
    if (!panels) return(4.5)
    max(5, ceiling(panels / 2) * (1.7 + 0.5 * n))
  }
  heat_height <- function(x, max_terms = 20L) {
    terms <- passing(x) |> distinct(collection, term_id) |> nrow()
    directions <- if (all(c("positive_direction", "negative_direction") %in% names(x))) n_distinct(x$contrast) else 0L
    max(4.5, 2.8 + 0.38 * min(terms, max_terms) + 0.23 * directions)
  }
  for (collection_name in c("KEGG", "GO_BP")) {
    ora_collection <- filter(ora, collection == collection_name)
    gsea_collection <- filter(gsea, collection == collection_name)
    if (!nrow(ora_collection) && !nrow(gsea_collection)) next
    for (scope in c("ancestor", "interactions")) {
      contrasts <- if (scope == "ancestor") ancestor else interactions
      ora_view <- filter(ora_collection, contrast %in% contrasts)
      gsea_view <- filter(gsea_collection, contrast %in% contrasts)
      ora_n <- if (scope == "ancestor") 6L else 4L
      gsea_n <- if (scope == "ancestor") 10L else 5L
      save_pair(plot_step05_ora(ora_view, policy = "primary", n_per_panel = ora_n),
                paste0("step05_", collection_name, "_", scope, "_ORA"),
                width = 15, height = ora_height(ora_view, ora_n), method = "ORA",
                collection = collection_name, scope = scope, candidates = nrow(passing(ora_view)))
      save_pair(plot_step05_gsea(gsea_view, policy = "primary", n_per_panel = gsea_n),
                paste0("step05_", collection_name, "_", scope, "_GSEA"),
                width = 15, height = gsea_height(gsea_view, gsea_n), method = "GSEA",
                collection = collection_name, scope = scope, candidates = nrow(passing(gsea_view)))
    }
    focus_view <- filter(gsea_collection, contrast %in% focus)
    save_pair(plot_step05_nes_heatmap(focus_view, contrasts = focus, policy = "primary", max_terms = 20L),
              paste0("step05_", collection_name, "_focus_NES_heatmap"),
              width = 14, height = heat_height(focus_view), method = "GSEA NES heatmap",
              collection = collection_name, scope = "ancestor plus three interactions",
              candidates = nrow(passing(focus_view)))
  }
  kegg_all <- filter(gsea, collection == "KEGG", contrast %in% all_contrasts)
  save_pair(plot_step05_nes_heatmap(kegg_all, contrasts = all_contrasts,
                                   policy = "primary", max_terms = 20L),
            "step05_KEGG_all10_NES_heatmap", width = 18, height = heat_height(kegg_all),
            method = "GSEA NES heatmap", collection = "KEGG", scope = "all ten contrasts",
            candidates = nrow(passing(kegg_all)))

  # Two transparent examples chosen by study FDR, then contrast/pathway as ties.
  # Their running curves use the exact cached ranks and gene sets used by GSEA.
  curves <- gsea |> filter(collection == "KEGG", contrast %in% focus,
                           is.finite(FDR_study), FDR_study >= 0, FDR_study < 0.1,
                           is.finite(NES)) |>
    arrange(FDR_study, contrast, pathway) |> slice_head(n = 2L)
  if (nrow(curves) && !requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea is required to reproduce the selected running-score curves.")
  }
  for (i in seq_len(nrow(curves))) {
    row <- curves[i, ]
    cache_path <- file.path(result_dir, "cache", paste0("primary__", row$contrast, ".rds"))
    if (!file.exists(cache_path)) stop("Missing exact GSEA cache for a running-score curve: ", cache_path)
    cache <- readRDS(cache_path)
    if (is.null(cache$ranks) || is.null(cache$pathways[[row$pathway]])) {
      stop("GSEA cache lacks the exact ranks or requested pathway: ", row$pathway)
    }
    # Match the canonical tie order used by the fgsea analysis wrapper.
    ranks <- cache$ranks[order(-cache$ranks, names(cache$ranks), method = "radix")]
    curve <- fgsea::plotEnrichment(cache$pathways[[row$pathway]], ranks, gseaParam = 1) +
      labs(title = str_wrap(paste0("KEGG: ", row$term_description), 65),
           subtitle = str_wrap(paste0(row$contrast_label, "; NES = ", sprintf("%.2f", row$NES),
                                     "; study FDR = ", format(row$FDR_study, digits = 3)), 90),
           x = "Position in the declared signed-statistic ranking", y = "Running enrichment score",
           caption = str_wrap(paste0(row$direction, ". Exact cached ranks and set membership; example chosen by lowest study FDR, not by a biological-priority claim."), 110)) +
      theme(plot.caption = element_text(hjust = 0), panel.grid.minor = element_blank())
    save_pair(curve, paste0("step05_KEGG_running_curve_", i), width = 11, height = 5.7,
              method = "GSEA running-score example", collection = "KEGG", scope = row$contrast,
              candidates = 1L)
  }
  manifest <- bind_rows(manifest)
  write_tsv(manifest, file.path(REPO_ROOT, "reports", "step05_figure_manifest.tsv"))
  message("STEP 5 FIGURES COMPLETE: ", nrow(manifest), " PNG/PDF pairs saved.")
  invisible(manifest)
}

render_step05_enrichment()
source(file.path(REPO_ROOT,"scripts","05_figures","05_render_module_evidence.R"))
