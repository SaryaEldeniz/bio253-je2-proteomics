# Synthetic checks for selection, interaction signs and missingness handling.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
passed <- character()
check_result <- function(label, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", label)
  passed <<- c(passed, label)
}
plan <- step04_contrast_plan()
check_result("Ten nonoptional comparisons, including all three interactions",
             nrow(plan) == 10 && sum(plan$contrast_type == "interaction") == 3)
check_result("Positive kinase interaction names the ancestor response",
             plan$positive_direction[plan$contrast == "interaction_stk_vs_anc_at_TSB_vs_SN"] == "SN-TSB change more positive in anc")
test <- tibble(entry_class = c("cds", "cds", "cds", "contaminant", "cds", "cds"),
                diff = c(-2, 1, 2, 3, 2, NA_real_), FDR = c(.01, .01, .1, .001, .01, .01),
                estimate_type = c("observed", "observed", "observed", "observed", "lod_imputed", "observed"),
                all_relevant_groups_2of3 = TRUE, any_relevant_group_unobserved = FALSE,
                positive_direction = "Higher in TSB", negative_direction = "Higher in SN")
classified <- step04_classify(test)
check_result("Strict FDR/effect cutoffs and CDS filter are respected", identical(which(classified$selected_hit), c(1L, 5L)))
check_result("Negative environmental effect means higher in SN", classified$direction[1] == "Higher in SN")
check_result("A LOD source flag can coexist with observed relevant groups",
             classified$selected_hit[5] && !classified$conservative_hit[5] &&
               classified$selection_note[5] == "LOD flag in source model; relevant groups have >=2/3 values")
unknown <- test[1, ]
unknown$estimate_type <- NA_character_
unknown <- step04_classify(unknown)
check_result("Unknown source flags are not silently labelled observed",
             !unknown$quantitative_support && unknown$selection_note == "Other/unknown source model; inspect support")
ab <- expand_grid(protein_key = c("measured", "missing"), group = c("anc_TSB", "anc_SN"), replicate = 1:3) |>
  mutate(sample = paste(group, replicate), normalized_abundance = case_when(
    protein_key == "missing" & group == "anc_TSB" ~ NA_real_, group == "anc_TSB" ~ 2, TRUE ~ 5))
support <- step04_group_support(ab, filter(plan, contrast == "TSB_vs_SN_at_anc"), "protein_key")$support
check_result("Fully observed group means preserve TSB-minus-SN direction",
             support$observed_group_contrast[support$feature_id == "measured"] == -3)
check_result("An unmeasured group remains NA rather than becoming zero",
             is.na(support$observed_group_contrast[support$feature_id == "missing"]) &&
               support$any_relevant_group_unobserved[support$feature_id == "missing"])
writeLines(c("# Step 4 code checks", "", paste(length(passed), "checks passed."), "", paste0("- PASS: ", passed)),
           file.path(REPO_ROOT, "reports", "step04_code_checks.md"))
message(length(passed), "/", length(passed), " Step 4 code checks passed.")
