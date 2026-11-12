# Step 1: validate the design, document contrasts, and draw Figure 1.
# Run from the project root. No differential analysis is performed here.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

out_dir <- file.path(REPO_ROOT, "results", "step01")
for (path in c(out_dir, file.path(REPO_ROOT, "figures"), file.path(REPO_ROOT, "reports"))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
checks <- tibble(check = character(), passed = logical(), detail = character())
record_check <- function(name, passed, detail) {
  checks <<- add_row(checks, check = name, passed = isTRUE(passed), detail = detail)
}

# 1. Read the course inputs. Keep contaminants for the inventory so counts are
# honest; biological filtering belongs in later modules.
samples <- read_samples()
abundance <- list(total = read_abundances("total", drop_contaminants = FALSE),
                  phospho_STYH = read_abundances("phospho_STYH", drop_contaminants = FALSE))
dea <- list(total = read_dea("total", drop_contaminants = FALSE),
            phospho_STYH = read_dea("phospho_STYH", drop_contaminants = FALSE))
sites <- read_sites("STYH", drop_contaminants = FALSE)
metadata_columns <- c("sample", "genotype", "medium", "replicate", "group")
expected <- expand_grid(genotype = GENOTYPES, medium = MEDIA, replicate = 1:3) |>
  mutate(sample = paste0("JE2", genotype, "_", if_else(medium == "SNpulse", "pulseSN", medium), "_", replicate),
         group = paste(genotype, medium, sep = "_")) |>
  select(all_of(metadata_columns))
actual <- samples |>
  mutate(across(c(genotype, medium, group), as.character)) |>
  select(all_of(metadata_columns))
record_check("Sample sheet: 27 unique nonmissing IDs", nrow(samples) == 27 &&
               !anyNA(samples$sample) && !anyDuplicated(samples$sample), paste(nrow(samples), "rows"))
record_check("Sample sheet matches the full 3 x 3 x 3 design",
             nrow(anti_join(expected, actual, by = metadata_columns)) == 0 &&
               nrow(anti_join(actual, expected, by = metadata_columns)) == 0,
             "Expected genotypes anc/stk/stp; media TSB/SNpulse/SN; replicates 1-3")
for (field in c("sample", "raw_file_total", "sample_name_in_excel")) {
  parsed <- parse_sample_name(samples[[field]])
  record_check(paste("Names agree after parsing:", field),
               identical(parsed$sample, samples$sample),
               "Canonical pulse sample IDs retain pulseSN; medium/group use SNpulse")
}

# 2. Cross-check layer metadata and feature coverage. Feature counts may differ
# between abundance and DEA exports; report the differences instead of guessing.
layer_samples <- list(); coverage <- list(); differences <- list(); contrasts <- list()
for (layer in names(abundance)) {
  ab <- abundance[[layer]]; de <- dea[[layer]]
  feature <- if (layer == "total") "protein_key" else "site_id"
  labels <- ab |>
    mutate(across(c(genotype, medium, group), as.character)) |>
    select(all_of(c(metadata_columns, "reinjected"))) |> distinct()
  record_check(paste(layer, "sample metadata"), nrow(labels) == 27 &&
                 nrow(anti_join(select(labels, all_of(metadata_columns)), expected,
                                by = metadata_columns)) == 0,
               "Every measurement row has the expected canonical labels")
  record_check(paste(layer, "unique feature x sample keys"),
               !anyDuplicated(ab[c(feature, "sample")]), "No repeated feature/sample rows")
  record_check(paste(layer, "complete feature x sample grid"),
               nrow(ab) == n_distinct(ab[[feature]]) * 27,
               "Complete row grid includes NA intensities; it does not imply complete detection")
  record_check(paste(layer, "unique feature x contrast keys"),
               !anyDuplicated(de[c(feature, "contrast")]), "No repeated feature/contrast rows")
  record_check(paste(layer, "nonmissing contrast names"),
               !anyNA(de$contrast) && all(nzchar(de$contrast)), "Every DEA row names its comparison")
  record_check(paste(layer, "nonmissing feature IDs"),
               !anyNA(ab[[feature]]) && !anyNA(de[[feature]]), "Feature IDs required for joins")
  expected_reinjected <- if (layer == "total") character() else c("JE2stk_pulseSN_1", "JE2stk_pulseSN_2")
  record_check(paste(layer, "layer-specific reinjections"), !anyNA(labels$reinjected) &&
                 setequal(labels$sample[labels$reinjected], expected_reinjected),
               paste(length(expected_reinjected), "flagged sample labels expected"))
  layer_samples[[layer]] <- mutate(labels, layer = layer, .before = 1)
  ab_keys <- distinct(ab, feature_id = .data[[feature]], is_contaminant)
  de_keys <- distinct(de, feature_id = .data[[feature]], is_contaminant)
  ab_only <- anti_join(ab_keys, de_keys, by = "feature_id")
  de_only <- anti_join(de_keys, ab_keys, by = "feature_id")
  record_check(paste(layer, "DEA features occur in abundance table"), nrow(de_only) == 0,
               paste(nrow(ab_only), "abundance-only features; identities saved locally"))
  differences[[layer]] <- bind_rows(mutate(ab_only, layer = layer, location = "abundance_only"),
                                    mutate(de_only, layer = layer, location = "DEA_only"))
  coverage[[layer]] <- tibble(layer, feature_unit = if (layer == "total") "protein entries" else "sites",
                              abundance_features = nrow(ab_keys), DEA_features = nrow(de_keys),
                              abundance_only = nrow(ab_only), DEA_only = nrow(de_only),
                              abundance_contaminant_features = sum(ab_keys$is_contaminant),
                              DEA_contaminant_features = sum(de_keys$is_contaminant),
                              samples = n_distinct(ab$sample), contrasts = n_distinct(de$contrast))
  contrasts[[layer]] <- de |> count(contrast, name = "tested_rows") |> mutate(layer = layer)
}
layer_samples <- bind_rows(layer_samples)
coverage <- bind_rows(coverage)
contrast_inventory <- bind_rows(contrasts)
record_check("Both DEA layers share the same 33 contrast names",
             length(unique(dea$total$contrast)) == 33 &&
               setequal(dea$total$contrast, dea$phospho_STYH$contrast),
             "Condition-specific and pooled contrasts are inventoried separately")
site_labels <- sites |>
  mutate(across(c(genotype, medium, group), as.character)) |>
  select(all_of(c(metadata_columns, "reinjected"))) |> distinct()
phospho_labels <- layer_samples |> filter(layer == "phospho_STYH") |> select(-layer)
record_check("FragPipe sites and phospho abundance metadata agree",
             nrow(site_labels) == 27 &&
               nrow(anti_join(site_labels, phospho_labels, by = names(site_labels))) == 0,
             "Includes the same two phospho-only reinjection flags")

# 3. Make a sign-explicit comparison plan. Coefficients concern group means.
plan <- project_i_contrasts() |>
  mutate(formula = vapply(contrast, contrast_formula, character(1)),
         interpretation = vapply(contrast, explain_contrast, character(1)))
for (layer in names(dea)) {
  record_check(paste(layer, "all planned contrasts exist"),
               all(plan$contrast %in% dea[[layer]]$contrast), paste(nrow(plan), "planned contrasts"))
}
weights <- bind_rows(lapply(plan$contrast, function(name) {
  values <- contrast_weights(name)
  tibble(contrast = name, group = names(values), coefficient = unname(values))
}))
record_check("Planned contrast coefficients sum to zero",
             all(tapply(weights$coefficient, weights$contrast, sum) == 0),
             "Simple comparisons use two cells; interactions use four cells")
coefficient_checks <- weights |>
  group_by(contrast) |>
  summarise(n_nonzero = sum(coefficient != 0), valid_values = all(coefficient %in% c(-1, 0, 1)), .groups = "drop")
record_check("Planned contrasts have the expected nonzero coefficients",
             all(coefficient_checks$valid_values) &&
               all(coefficient_checks$n_nonzero == if_else(str_starts(coefficient_checks$contrast, "interaction_"), 4L, 2L)),
             "Exactly two +/-1 coefficients for simple comparisons; four for interactions")

# Detailed metadata outputs are local-only according to .gitignore.
group_counts <- layer_samples |>
  count(layer, genotype, medium, group, name = "n_replicates")
record_check("Every group has three sample labels in each layer",
             nrow(group_counts) == 18 && all(group_counts$n_replicates == 3), "9 groups x 2 layers")
write_tsv(layer_samples, file.path(out_dir, "sample_manifest.tsv"))
write_tsv(group_counts, file.path(out_dir, "group_counts.tsv"))
write_tsv(coverage, file.path(out_dir, "feature_coverage.tsv"))
write_tsv(bind_rows(differences), file.path(out_dir, "feature_coverage_differences.tsv"))
write_tsv(contrast_inventory, file.path(out_dir, "all_contrasts.tsv"))
write_tsv(plan, file.path(out_dir, "planned_contrasts.tsv"))
write_tsv(weights, file.path(out_dir, "contrast_coefficients.tsv"))
write_tsv(checks, file.path(out_dir, "validation_checks.tsv"))
if (!all(checks$passed)) {
  print(filter(checks, !passed), width = Inf)
  stop("Design validation failed. Inspect results/step01/validation_checks.tsv.")
}

# 4. Save the figure and a readable report, including limits of these checks.
figure <- plot_experimental_design(layer_samples)
ggsave(file.path(REPO_ROOT, "figures", "figure01_design.png"), figure,
       width = 10.5, height = 4.7, units = "in", dpi = 300, bg = "white")
coverage_lines <- apply(coverage, 1, function(row) {
  sprintf("| %s | %s | %s | %s | %s | %s |", row[["layer"]], row[["abundance_features"]],
          row[["DEA_features"]], row[["abundance_only"]], row[["samples"]], row[["contrasts"]])
})
plan_lines <- paste0("| `", plan$contrast, "` | ", plan$priority, " | ", plan$question, " |")
alias_count <- sum(samples$sample != samples$sample_name_in_excel)
report <- c("# Step 1: experimental design and comparison plan", "",
            paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
            sprintf("**Result: %d/%d validation checks passed.**", sum(checks$passed), nrow(checks)), "",
            "## Verified design", "",
            "- 27 unique sample labels in the sample sheet and in each measurement layer.",
            "- 3 genotypes (anc, stk, stp), 3 media (TSB, SNpulse, SN), 3 replicate labels per group.",
            "- TSB: control; SNpulse: 15-minute PASN exposure; SN: approximately 15-hour PASN exposure.",
            "- PASN is 30% cell-free, sterile-filtered PAO1 spent medium; this is not a live co-culture.",
            sprintf("- %d shortened Excel sample labels map back to their canonical JE2 sample labels.", alias_count),
            "- Total layer: no reinjection flags. Phospho layer: JE2stk_pulseSN_1 and _2 flagged.",
            "- These flags do not increase the replicate count. Matching replicate numbers do not establish pairing.", "",
            "![Experimental design](../figures/figure01_design.png)", "",
            "## Dataset inventory before contaminant filtering", "",
            "| Layer | Abundance features | DEA features | Abundance-only | Samples | Contrasts |",
            "|---|---:|---:|---:|---:|---:|", coverage_lines, "",
            "Abundance and DEA feature counts differ. The omitted identities are saved locally; Step 1 does not infer why they were omitted.",
            "A complete feature-by-sample row grid includes missing (NA) values and does not mean every feature was detected.", "",
            "## Direction and planned comparisons", "",
            "`A_vs_B` means `log2(A) - log2(B)`. Negative `TSB_vs_SN` means higher in SN.",
            "Interaction: `(stk_TSB - anc_TSB) - (stk_SN - anc_SN)` for the kinase-versus-ancestor TSB-versus-SN contrast.",
            "The stp version follows the same formula with stp in place of stk. These are differences between group means, not paired-sample subtractions.", "",
            "| Contrast | Priority | Question |", "|---|---|---|", plan_lines, "",
            "Full coefficient tables and interpretation text are in the local results/step01 directory.",
            "Pooled contrasts are inventoried but not assigned an unverified weighting formula.", "",
            "## What remains for Step 2", "",
            "Check actual detection, missingness, replicate agreement/PCA, and Stk1/Stp1 measurements. Correct metadata do not establish good measurement quality.",
            "The course design notes reduced Stk1 in the stp strain; verify per-sample values before interpreting mutant effects.",
            "No hit selection, enrichment, site-protein adjustment, Venn diagram or motif analysis was run in Step 1.", "",
            "## Reproduce", "", "```r", 'source("scripts/01_data/01_validate_design.R")', "```", "",
            "See step01_inputs.md for input hashes and step01_session.txt for package versions.")
writeLines(report, file.path(REPO_ROOT, "reports", "step01_design.md"))

# Hashes establish which local inputs were used without copying measurements.
input_names <- c("samples.tsv", "protein_annotation.tsv.gz", "abundances_total.tsv.gz",
                 "abundances_phospho_STYH.tsv.gz", "dea_total.tsv.gz",
                 "dea_phospho_STYH.tsv.gz", "sites_STYH.tsv.gz")
input_paths <- vapply(input_names, derived_path, character(1))
input_lines <- sprintf("| derived/%s | %s | `%s` |", input_names,
                       file.info(input_paths)$size, unname(tools::md5sum(input_paths)))
writeLines(c("# Step 1 input provenance", "", "Input location is configured locally and excluded from Git.", "",
             "| Input | Bytes | MD5 |", "|---|---:|---|", input_lines),
           file.path(REPO_ROOT, "reports", "step01_inputs.md"))
writeLines(capture.output(sessionInfo()), file.path(REPO_ROOT, "reports", "step01_session.txt"))
message(sprintf("STEP 1 COMPLETE: %d/%d checks passed; figure and report saved.", sum(checks$passed), nrow(checks)))
print(coverage, width = Inf)
