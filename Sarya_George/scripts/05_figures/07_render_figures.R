# Run from the repository root or from the Sarya_George project directory.
group_root <- if (file.exists("R/analysis/step07_parent_comparison.R")) "." else "Sarya_George"
source(file.path(group_root, "R/analysis/step07_parent_comparison.R"))
source(file.path(group_root, "R/figures/step07_plots.R"))

data_dir <- Sys.getenv("JE2_DATA_DIR")
data_config <- file.path(group_root, "config", "data_dir.txt")
if (!nzchar(data_dir) && file.exists(data_config)) {
  data_dir <- trimws(readLines(data_config, n = 1L, warn = FALSE)[1L])
}
if (is.na(data_dir) || !nzchar(data_dir)) data_dir <- file.path(group_root, "..", "data")
phospho <- read.delim(file.path(data_dir, "derived", "dea_phospho_STYH.tsv.gz"))
total <- read.delim(file.path(data_dir, "derived", "dea_total.tsv.gz"))

contrasts <- c(
  "interaction_stk_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_stk_at_TSB_vs_SN"
)
titles <- c("stk vs anc", "stp vs anc", "stp vs stk")
results <- summarize_step07_parent_comparison(phospho, total, contrasts)

interaction_table <- results[c(
  "contrast", "significant_sites", "same_direction", "opposite_direction",
  "parent_below_thresholds", "parent_unavailable"
)]
interaction_table$contrast <- titles
print(interaction_table, row.names = FALSE)

report_dir <- file.path(group_root, "reports")
figure_dir <- file.path(group_root, "figures", "step07")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.table(interaction_table, file.path(report_dir, "step07_phospho_parent_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
plot_step07_volcano(phospho, contrasts, titles,
                    file.path(figure_dir, "figure07_phospho_interaction_volcanoes.png"))
