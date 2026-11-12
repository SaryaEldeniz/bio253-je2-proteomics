source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

render_step04_figures <- function() {
plan <- step04_contrast_plan(); map <- read_identifier_map()
out_root <- file.path(REPO_ROOT, "results", "step04")
figure_dir <- file.path(REPO_ROOT, "figures", "step04")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
all_tables <- list(total=readr::read_tsv(file.path(out_root,"total","annotated_all.tsv.gz"),show_col_types=FALSE))
all_summaries <- list(total=readr::read_tsv(file.path(out_root,"total","hit_summary.tsv"),show_col_types=FALSE))
all_abundances <- list(total=read_abundances("total",drop_contaminants=FALSE))
# Figures focus on the ordinary total proteome. Phosphosite reliability is Step 7.
total <- filter(all_tables$total, eligible_cds)
environmental <- filter(plan, contrast_type == "environmental")
interactions <- filter(plan, contrast_type == "interaction")
volcano_env <- plot_step04_volcano(total, environmental$contrast, environmental$contrast_label) +
  labs(title = "Total proteome: long PASN response within each genotype",
       x = "Source log2(TSB / SN): negative = higher in SN; positive = higher in TSB")
volcano_interactions <- plot_step04_volcano(total, interactions$contrast,
                                          c("stk vs anc", "stp vs anc", "stp vs stk")) +
  labs(title = "Total proteome: genotype differences in PASN response",
       subtitle = "Define R = log2(SN / TSB). Panels show R_anc - R_stk, R_anc - R_stp, and R_stk - R_stp.",
       x = "Interaction effect (log2); inspect within-genotype responses to interpret its sign")
ggsave(file.path(figure_dir, "figure02_total_PASN_volcanoes.png"), volcano_env,
       width = 12, height = 5.6, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure03_total_interaction_volcanoes.png"), volcano_interactions,
       width = 12, height = 5.6, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure04_total_hit_counts.png"), plot_step04_hit_counts(all_summaries$total),
       width = 12.5, height = 11, units = "in", dpi = 300, bg = "white")

# Show predefined candidates and an observed-supported, predeclared ranking example
# from each focal interaction. Plot measurements before making biological claims.
candidate_config <- read_tsv(file.path(REPO_ROOT, "config", "interesting_proteins.tsv"), show_col_types = FALSE)
candidate_lookup <- lookup_je2(candidate_config$query)
candidate_keys <- unique(candidate_lookup$protein_key[candidate_lookup$lookup_status == "matched"])
label_map <- map |> transmute(protein_key, label = paste(display_name, locus_tag, sep = "\n"))
candidate_plot <- plot_step04_samples(all_abundances$total, candidate_keys, label_map)
ggsave(file.path(figure_dir, "figure05_candidate_measurements.png"), candidate_plot,
       width = 10.5, height = 9, units = "in", dpi = 300, bg = "white")
top_examples <- total |> filter(contrast_type == "interaction", conservative_hit) |>
  group_by(contrast) |> arrange(FDR, desc(abs(diff)), feature_id, .by_group = TRUE) |>
  slice_head(n = 1) |> ungroup()
write_tsv(top_examples, file.path(out_root, "total", "interaction_plot_examples.tsv"))
if (nrow(top_examples)) {
  ggsave(file.path(figure_dir, "figure06_interaction_examples.png"),
         plot_step04_samples(all_abundances$total, unique(top_examples$protein_key), label_map),
         width = 10.5, height = 4.5, units = "in", dpi = 300, bg = "white")
}
candidate_results <- all_tables$total |> filter(protein_key %in% candidate_keys)
write_tsv(candidate_results, file.path(out_root, "total", "candidate_comparisons.tsv"))

}
render_step04_figures()
