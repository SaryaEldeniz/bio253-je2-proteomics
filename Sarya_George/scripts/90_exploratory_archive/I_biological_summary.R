# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
targets <- tibble::tribble(
  ~label, ~locus_tag,
  "Stp1", "B7H15_06350", "Stk1", "B7H15_06355",
  "AhpC", "B7H15_02095", "SodM", "B7H15_08530",
  "KatA", "B7H15_07015", "Hla", "B7H15_06050",
  "CrtN", "B7H15_14250", "TcyA", "B7H15_13375",
  "TcyB", "B7H15_13370", "TcyC", "B7H15_13365"
)
ab <- read_abundances("total") |> inner_join(targets, by = "locus_tag")
p <- plot_protein(targets$locus_tag, ab, max_proteins = 12)
# Sample values are checked before interpreting any model difference.
ggsave("results/exploratory/I_marker_sample_check.png",
       p, width = 14, height = 10, dpi = 120)
summary <- ab |> group_by(label, locus_tag, genotype, medium) |>
  summarise(n_detected = sum(!is.na(normalized_abundance)),
            mean_log2 = if(all(is.na(normalized_abundance))) NA_real_ else
              mean(normalized_abundance, na.rm = TRUE), .groups = "drop")
write_tsv(summary, "results/exploratory/I_marker_abundance_summary.tsv")
dea <- read_dea("total", c("TSB_vs_SN_at_anc", "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp")) |>
  inner_join(targets, by = "locus_tag") |>
  select(label, locus_tag, product, contrast, diff, FDR, estimate_type)
write_tsv(dea, "results/exploratory/I_marker_dea.tsv")
print(filter(dea, contrast == "TSB_vs_SN_at_anc"), n = 20, width = 180)
cat("\nREGULATOR SAMPLE SUMMARY\n")
print(filter(summary, label %in% c("Stk1", "Stp1")), n = 20, width = 160)
