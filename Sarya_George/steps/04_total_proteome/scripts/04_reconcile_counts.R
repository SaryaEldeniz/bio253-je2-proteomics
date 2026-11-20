# Reconcile the lecturer screenshot with Step 4 selections without changing them.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
x <- read_tsv(file.path(REPO_ROOT, "results", "step04", "total", "annotated_all.tsv.gz"), show_col_types = FALSE)
names <- c("interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_stk_at_TSB_vs_SN")
x <- filter(x, contrast %in% names)
counts <- x |> group_by(contrast) |>
  summarise(export_entries = n(), source_FDR_only = sum(FDR < .1, na.rm = TRUE),
            source_FDR_and_effect = sum(FDR < .1 & abs(diff) > 1, na.rm = TRUE),
            ordinary_CDS_FDR_only = sum(eligible_cds & FDR < .1, na.rm = TRUE),
            selected = sum(selected_hit),
            selected_source_observed = sum(selected_hit & estimate_type == "observed"),
            conservative = sum(conservative_hit), .groups = "drop")
excluded <- x |> filter(FDR < .1, abs(diff) > 1, !eligible_cds) |>
  count(contrast, entry_class, name = "excluded_hits")
support <- x |> filter(selected_hit) |>
  mutate(reason = case_when(estimate_type != "observed" & !all_relevant_groups_2of3 ~ "LOD source flag AND sparse relevant group",
                            estimate_type != "observed" ~ "LOD source flag only",
                            !all_relevant_groups_2of3 ~ "Sparse relevant group only",
                            TRUE ~ "Retained conservative hit")) |>
  count(contrast, reason, name = "entries")
out <- file.path(REPO_ROOT, "results", "step04", "reconciliation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv(counts, file.path(out, "counts.tsv"))
write_tsv(excluded, file.path(out, "excluded_classes.tsv"))
write_tsv(support, file.path(out, "conservative_breakdown.tsv"))
write_tsv(filter(x, FDR < .1, abs(diff) > 1, !eligible_cds), file.path(out, "excluded_hit_identities.tsv"))
write_tsv(filter(x, selected_hit, !conservative_hit), file.path(out, "nonconservative_hit_identities.tsv"))
print(counts, width = Inf); print(excluded, width = Inf); print(support, width = Inf)

# The screenshot's 2348-entry universe may identify the older August export.
older_path <- file.path(DATA_ROOT, "DEAs", "DEAs", "DEA_20260802_WUtotal_JE2_2f_vsn",
                        "Results_WU_total_JE2_2f", "DE_WUtotal_JE2_2f.xlsx")
if (file.exists(older_path)) {
  older <- parse_dea_excel(older_path) |> filter(contrast %in% names)
  older_counts <- older |> group_by(contrast) |>
    summarise(version = "August 2026 total-proteome export", export_entries = n(),
              FDR_only = sum(FDR < .1, na.rm = TRUE),
              FDR_and_effect = sum(FDR < .1 & abs(diff) > 1, na.rm = TRUE),
              model = paste(unique(modelName), collapse = ";"), .groups = "drop")
  write_tsv(older_counts, file.path(out, "older_export_counts.tsv"))
  print(older_counts, width = Inf)
}
