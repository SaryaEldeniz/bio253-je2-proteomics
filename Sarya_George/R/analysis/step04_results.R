# Step 4: annotate existing DEA results; do not refit or recalculate their FDRs.

step04_contrast_plan <- function() {
  plan <- project_i_contrasts() |> filter(priority != "optional")
  bind_rows(lapply(seq_len(nrow(plan)), function(i) {
    name <- plan$contrast[i]
    interaction <- str_match(name, "^interaction_(.+)_vs_(.+)_at_TSB_vs_SN$")
    environmental <- str_match(name, "^TSB_vs_SN_at_(.+)$")
    mutant <- str_match(name, "^(stk|stp)_vs_anc_at_(TSB|SN)$")
    if (!is.na(interaction[1, 1])) {
      a <- interaction[1, 2]; b <- interaction[1, 3]
      type <- "interaction"
      label <- paste0(a, " vs ", b, ": PASN response difference")
      positive <- paste0("SN-TSB change more positive in ", b)
      negative <- paste0("SN-TSB change more positive in ", a)
    } else if (!is.na(environmental[1, 1])) {
      type <- "environmental"
      label <- paste0(environmental[1, 2], ": TSB vs SN")
      positive <- "Higher in TSB"; negative <- "Higher in SN"
    } else if (!is.na(mutant[1, 1])) {
      type <- "mutant"
      label <- paste0(mutant[1, 2], " vs anc in ", mutant[1, 3])
      positive <- paste0("Higher in ", mutant[1, 2]); negative <- "Higher in anc"
    } else stop("Unrecognized Step 4 contrast: ", name)
    mutate(plan[i, ], contrast_type = type, contrast_label = label,
           positive_direction = positive, negative_direction = negative,
           formula = contrast_formula(name))
  }))
}

step04_group_support <- function(abundances, plan, feature_col) {
  group_counts <- abundances |>
    group_by(feature_id = .data[[feature_col]], group) |>
    summarise(n_samples = n_distinct(sample), n_observed = sum(is.finite(normalized_abundance)),
              mean_observed = if (any(is.finite(normalized_abundance))) {
                mean(normalized_abundance[is.finite(normalized_abundance)])
              } else NA_real_, .groups = "drop") |>
    mutate(group = as.character(group))
  # One coefficient row per relevant group. NA measurements remain NA.
  coefficients <- bind_rows(lapply(plan$contrast, function(name) {
    values <- contrast_weights(name)
    tibble(contrast = name, group = names(values), coefficient = unname(values)) |>
      filter(coefficient != 0)
  }))
  support <- inner_join(coefficients, group_counts, by = "group", relationship = "many-to-many") |>
    group_by(contrast, feature_id) |>
    summarise(n_relevant_groups = n(), n_observed_min_relevant = min(n_observed),
              n_observed_max_relevant = max(n_observed),
              any_relevant_group_unobserved = any(n_observed == 0),
              all_relevant_groups_2of3 = all(n_observed >= 2),
              group_detection_pattern = paste0(group, "=", n_observed, "/", n_samples, collapse = "; "),
              observed_group_contrast = if (all(n_observed > 0)) sum(coefficient * mean_observed) else NA_real_,
              .groups = "drop")
  list(groups = group_counts, coefficients = coefficients, support = support)
}

step04_classify <- function(x, fdr_cutoff = FDR_CUTOFF, lfc_cutoff = LFC_CUTOFF) {
  x |>
    mutate(eligible_cds = !is.na(entry_class) & entry_class == "cds",
           testable_for_selection = is.finite(diff) & is.finite(FDR),
           selected_hit = eligible_cds & testable_for_selection & FDR < fdr_cutoff & abs(diff) > lfc_cutoff,
           direction = case_when(!is.finite(diff) ~ "No finite effect estimate",
                                 diff > 0 ~ positive_direction, diff < 0 ~ negative_direction,
                                 TRUE ~ "Zero estimated difference"),
           quantitative_support = coalesce(estimate_type == "observed", FALSE) & coalesce(all_relevant_groups_2of3, FALSE),
           conservative_hit = selected_hit & quantitative_support,
           source_model_flag = case_when(estimate_type == "observed" ~ "Source model: observed",
                                         estimate_type == "lod_imputed" ~ "Source model: LOD-flagged",
                                         TRUE ~ "Source model: other/unknown"),
           selection_note = case_when(!eligible_cds ~ paste0("Excluded entry class: ", coalesce(entry_class, "unmapped")),
                                      !testable_for_selection ~ "Nonfinite effect or FDR",
                                      any_relevant_group_unobserved ~ "At least one relevant group has no measured values",
                                      !all_relevant_groups_2of3 ~ "Some relevant groups have fewer than 2/3 measured values",
                                      estimate_type == "lod_imputed" ~ "LOD flag in source model; relevant groups have >=2/3 values",
                                      estimate_type == "observed" ~ "Observed source model; >=2/3 values in each relevant group",
                                      TRUE ~ "Other/unknown source model; inspect support"))
}

annotate_step04 <- function(dea, map, plan, support, dataset, sites = NULL) {
  fields <- c("protein_key", "entry_class", "display_name", "uniprot_acc", "ortholog_quality",
              "eligible_uniprot_annotation", "n_JE2_good_all", "n_JE2_good_total_DEA",
              "string_id", "n_annotation_terms", "up_protein_name", "up_keywords", "up_go_ids", "up_kegg")
  feature_col <- if (dataset == "total") "protein_key" else "site_id"
  x <- dea |>
    filter(contrast %in% plan$contrast) |>
    mutate(feature_id = .data[[feature_col]], dataset = dataset) |>
    left_join(select(map, all_of(fields)), by = "protein_key", relationship = "many-to-one") |>
    left_join(plan, by = "contrast", relationship = "many-to-one") |>
    left_join(support, by = c("contrast", "feature_id"), relationship = "one-to-one")
  if (!is.null(sites)) {
    confidence <- sites |> distinct(site_id, best_loc_prob)
    if (anyDuplicated(confidence$site_id)) stop("Site localization metadata are not unique per site.")
    x <- left_join(x, confidence, by = "site_id", relationship = "many-to-one")
  }
  x |> step04_classify() |>
    relocate(dataset, contrast, contrast_label, contrast_type, feature_id, protein_key,
             locus_tag, display_name, uniprot_acc, entry_class, diff, FDR, estimate_type,
             direction, selected_hit, conservative_hit, n_observed_min_relevant,
             any_relevant_group_unobserved, group_detection_pattern)
}

step04_hit_summary <- function(x, plan) {
  directions <- bind_rows(transmute(plan, contrast, contrast_label, direction = positive_direction),
                           transmute(plan, contrast, contrast_label, direction = negative_direction))
  hits <- x |> filter(selected_hit) |>
    group_by(contrast, direction) |>
    summarise(n_hits = n(), n_hits_observed = sum(estimate_type == "observed"),
              n_hits_lod_imputed = sum(estimate_type == "lod_imputed"),
              n_hits_other = sum(!estimate_type %in% c("observed", "lod_imputed")),
              n_quantitative_support = sum(quantitative_support),
              n_zero_group = sum(any_relevant_group_unobserved), .groups = "drop")
  directions |> left_join(hits, by = c("contrast", "direction"), relationship = "one-to-one") |>
    mutate(across(starts_with("n_"), ~coalesce(.x, 0L))) |>
    arrange(match(contrast, plan$contrast))
}

# A beginner-friendly reader for the local Step 4 results.
read_comparison <- function(contrast, dataset = "total", hits_only = FALSE,
                             conservative_only = FALSE) {
  dataset <- match.arg(dataset, c("total", "phospho_STYH"))
  path <- file.path(REPO_ROOT, "results", "step04", dataset, "annotated_all.tsv.gz")
  if (!file.exists(path)) stop('First run source("scripts/04_analysis/04_assemble_comparisons.R").')
  x <- read_tsv(path, show_col_types = FALSE)
  if (length(contrast) != 1 || !contrast %in% x$contrast) stop("Choose one contrast from step04_contrast_plan().")
  x <- filter(x, .data$contrast == .env$contrast, eligible_cds)
  if (hits_only) x <- filter(x, selected_hit)
  # This switch explicitly selects conservative HITS, not every well-supported row.
  if (conservative_only) x <- filter(x, conservative_hit)
  arrange(x, FDR, desc(abs(diff)))
}
