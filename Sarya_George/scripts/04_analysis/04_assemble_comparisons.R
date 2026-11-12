# Step 4: reproduce annotated comparison tables and total-proteome figures.
# Existing effect estimates, p-values and FDRs are preserved, not refitted.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
plan <- step04_contrast_plan()
map <- read_identifier_map()
out_root <- file.path(REPO_ROOT, "results", "step04")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(REPO_ROOT, "figures", "step04")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(plan, file.path(out_root, "comparison_plan.tsv"))
all_tables <- list(); all_summaries <- list(); all_checks <- list(); all_abundances <- list()
site_metadata <- read_sites("STYH", drop_contaminants = FALSE)

for (dataset in c("total", "phospho_STYH")) {
  message("Assembling ", dataset, " ...")
  out_dir <- file.path(out_root, dataset)
  dir.create(file.path(out_dir, "contrasts"), recursive = TRUE, showWarnings = FALSE)
  source_dea <- read_dea(dataset, plan$contrast, drop_contaminants = FALSE)
  abundances <- read_abundances(dataset, drop_contaminants = FALSE)
  feature_col <- if (dataset == "total") "protein_key" else "site_id"
  support <- step04_group_support(abundances, plan, feature_col)
  x <- annotate_step04(source_dea, map, plan, support$support, dataset,
                       sites = if (dataset == "total") NULL else site_metadata)
  x$analysis_status <- if (dataset == "total") "Protein abundance results" else {
    "Preliminary phosphosite source selection; localization and parent-protein checks deferred to Step 7"
  }
  # Interaction effects are algebraically related; their FDRs remain separate source tests.
  effects <- x |> filter(contrast_type == "interaction") |>
    select(feature_id, contrast, diff) |> pivot_wider(names_from = contrast, values_from = diff)
  identity_error <- with(effects, abs(interaction_stp_vs_stk_at_TSB_vs_SN -
                                       (interaction_stp_vs_anc_at_TSB_vs_SN - interaction_stk_vs_anc_at_TSB_vs_SN)))
  expected_groups <- tibble(contrast = plan$contrast,
                           expected = if_else(plan$contrast_type == "interaction", 4L, 2L))
  counts <- x |> left_join(expected_groups, by = "contrast")
  checks <- tibble(dataset = dataset,
                  check = c("Annotation preserves row count", "Source effects unchanged", "Source FDR unchanged",
                            "Unique feature-contrast keys", "All rows have class metadata",
                            "Every relevant group has a support-count row", "Relevant groups have three sample labels",
                            "Source estimate labels are explicit", "Selected entries are ordinary CDS",
                            "Third interaction effect equals second minus first"),
                  passed = c(nrow(x) == nrow(source_dea), identical(x$diff, source_dea$diff),
                             identical(x$FDR, source_dea$FDR), !anyDuplicated(x[c("feature_id", "contrast")]),
                             !anyNA(x$entry_class), all(counts$n_relevant_groups == counts$expected),
                             all(support$groups$n_samples == 3), !anyNA(x$estimate_type),
                             all(x$entry_class[x$selected_hit] == "cds"),
                             all(identity_error[is.finite(identity_error)] < 1e-8)))
  write_tsv(checks, file.path(out_dir, "validation_checks.tsv"))
  if (!all(checks$passed)) { print(filter(checks, !passed)); stop("Step 4 checks failed for ", dataset) }
  summary <- step04_hit_summary(x, plan)
  inventory <- x |> group_by(contrast, contrast_label, contrast_type) |>
    summarise(source_entries = n(), ordinary_CDS_rows = sum(eligible_cds),
              finite_CDS_rows = sum(eligible_cds & testable_for_selection),
              selected_rows = sum(selected_hit), conservative_selected_rows = sum(conservative_hit),
              source_LOD_flag_with_all_relevant_3of3 = sum(eligible_cds & estimate_type == "lod_imputed" & n_observed_min_relevant == 3),
              source_observed_with_sparse_relevant_group = sum(eligible_cds & estimate_type == "observed" & !all_relevant_groups_2of3),
              .groups = "drop")
  write_tsv(x, file.path(out_dir, "annotated_all.tsv.gz"))
  write_tsv(filter(x, eligible_cds), file.path(out_dir, "ordinary_CDS_all.tsv.gz"))
  write_tsv(filter(x, !eligible_cds), file.path(out_dir, "excluded_entries_audit.tsv"))
  write_tsv(support$groups, file.path(out_dir, "group_observation_counts.tsv"))
  write_tsv(summary, file.path(out_dir, "hit_summary.tsv"))
  write_tsv(inventory, file.path(out_dir, "comparison_inventory.tsv"))
  for (name in plan$contrast) {
    tab <- filter(x, contrast == name, eligible_cds) |> arrange(FDR, desc(abs(diff)))
    prefix <- file.path(out_dir, "contrasts", name)
    write_tsv(tab, paste0(prefix, "_all.tsv"))
    write_tsv(filter(tab, selected_hit), paste0(prefix, "_hits.tsv"))
    write_tsv(filter(tab, conservative_hit), paste0(prefix, "_conservative_hits.tsv"))
    # Explicit direction lists remain in JE2 IDs. Do not silently collapse UniProt collisions.
    write_tsv(filter(tab, selected_hit, diff > 0), paste0(prefix, "_positive_hits.tsv"))
    write_tsv(filter(tab, selected_hit, diff < 0), paste0(prefix, "_negative_hits.tsv"))
  }
  all_tables[[dataset]] <- x; all_summaries[[dataset]] <- summary
  all_checks[[dataset]] <- checks; all_abundances[[dataset]] <- abundances
  print(inventory |> select(contrast, ordinary_CDS_rows, selected_rows, conservative_selected_rows), width = Inf)
}

source(file.path(REPO_ROOT,"scripts","05_figures","04_render_figures.R"))
write_tsv(bind_rows(all_checks), file.path(out_root,"validation_checks.tsv"))

# Summary tables for the report. All complete data-derived tables remain local.
inventory <- read_tsv(file.path(out_root, "total", "comparison_inventory.tsv"), show_col_types = FALSE)
report_table <- inventory |> left_join(select(plan, contrast, question), by = "contrast") |>
  select(contrast, ordinary_CDS_rows, selected_rows, conservative_selected_rows)
table_lines <- apply(report_table, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
phospho_inventory <- read_tsv(file.path(out_root, "phospho_STYH", "comparison_inventory.tsv"), show_col_types = FALSE)
phospho_lines <- apply(select(phospho_inventory, contrast, ordinary_CDS_rows, selected_rows, conservative_selected_rows),
                       1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
flag_example <- filter(inventory, contrast == "TSB_vs_SN_at_anc")
report <- c("# Step 4: annotated comparisons and total-proteome figures", "",
            paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
            sprintf("%d/%d actual-input checks passed across the two layers.", sum(bind_rows(all_checks)$passed), nrow(bind_rows(all_checks))), "",
            "The supplied DEA was annotated, not refitted. Source diff, p-values and FDR are preserved. Step 2 remains skipped at the user's request.", "",
            "## Total-proteome selection", "",
            "Ordinary protein entries use metadata entry_class == cds. Selected means source FDR < 0.1 AND |source diff| > 1.",
            "Conservative selected means the same source selection, an observed source-model label, and >=2/3 measured values in every relevant group.", "",
            "| Contrast | Ordinary protein entries | Selected | Conservative selected |",
            "|---|---:|---:|---:|", table_lines, "",
            "These are protein entries in the JE2 measurement space, not deduplicated UniProt IDs. FDRs were not recalculated after class filtering or across contrasts.", "",
            "## Read the directions", "",
            "Negative TSB_vs_SN = higher in SN. Positive TSB_vs_SN = higher in TSB.",
            "Define R_g = log2(SN_g) - log2(TSB_g). The three interactions are R_anc-R_stk, R_anc-R_stp and R_stk-R_stp.",
            "A positive interaction means the reference genotype has a more positive SN-TSB change; it does not simply mean higher abundance in that genotype.",
            "The third interaction effect equals the second minus the first. Its p-value, FDR and standard error are separate source results.", "",
            "## Source model flags versus actual measurements", "",
            paste("In the ancestor TSB-vs-SN ordinary-protein table,", flag_example$source_LOD_flag_with_all_relevant_3of3,
                  "LOD-flagged entries nevertheless have 3/3 measurements in both relevant groups."),
            "The source label is global/model-wide in this export. Use group_detection_pattern, n_observed_min_relevant and any_relevant_group_unobserved to examine each comparison.",
            "Do not call every LOD-flagged result a presence/absence effect in the displayed groups. Conversely an observed model label does not guarantee dense replicate support.", "",
            "## Figures", "",
            "![Within-genotype PASN volcanoes](../figures/step04/figure02_total_PASN_volcanoes.png)", "",
            "![Genotype-by-medium interaction volcanoes](../figures/step04/figure03_total_interaction_volcanoes.png)", "",
            "![Directional selected counts](../figures/step04/figure04_total_hit_counts.png)", "",
            "![Predefined candidate measurements](../figures/step04/figure05_candidate_measurements.png)", "",
            "The predefined candidates are the editable Step 3 list, not a list chosen to maximize significance. Figure 6, when present, shows the lowest-FDR conservative selected entry per interaction, with effect size and key as tie-breakers. These examples are illustrative and not literature-validated biological priorities.", "",
            "## Preliminary phosphosite tables", "",
            "The same ten contrasts have annotated site-level tables. The counts below are sites whose parent is an ordinary CDS, not counts of proteins.",
            "These source selections are provisional. No localization threshold or parent-protein adjustment has been applied; those are Step 7. Do not use these as final regulated phosphosite lists or sequence-logo inputs.", "",
            "| Contrast | CDS-parent site rows | Preliminary selected | Conservative source/support subset |",
            "|---|---:|---:|---:|", phospho_lines, "",
            "## Use the tables in RStudio", "", "```r",
            'source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")',
            'ancestor <- read_comparison("TSB_vs_SN_at_anc")',
            'hits <- read_comparison("interaction_stk_vs_anc_at_TSB_vs_SN", hits_only = TRUE)',
            'View(hits)', "```", "",
            "Full per-contrast tables, directional lists, conservative subsets and exclusions are under results/step04. The explanatory guide is docs/step04_interpretation.md.", "",
            "## Limits and next step", "",
            "These data compare responses to cell-free mixed spent medium. The documented Stk1 reduction in stp complicates pure phosphatase-loss interpretations.",
            "No pathway enrichment, network-module inference, final phosphosite regulation, Venn or motif analysis was performed here. Next is Step 5 functional enrichment using the annotated ordinary-protein tables and appropriate measured backgrounds.")
writeLines(report, file.path(REPO_ROOT, "reports", "step04_comparisons.md"))
input_paths <- c(derived_path("dea_total.tsv.gz"), derived_path("dea_phospho_STYH.tsv.gz"),
                 derived_path("abundances_total.tsv.gz"), derived_path("abundances_phospho_STYH.tsv.gz"),
                 derived_path("sites_STYH.tsv.gz"), file.path(REPO_ROOT, "results", "step03", "id_map.tsv"))
writeLines(c("# Step 4 inputs", "", "| File | MD5 |", "|---|---|",
             sprintf("| %s | `%s` |", basename(input_paths), unname(tools::md5sum(input_paths)))),
           file.path(REPO_ROOT, "reports", "step04_inputs.md"))
writeLines(capture.output(sessionInfo()), file.path(REPO_ROOT, "reports", "step04_session.txt"))
message("STEP 4 COMPLETE: tables, total-proteome figures and report saved.")
