# Focused checks for the new code. No private data or testthat package required.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
results <- character()
check_code <- function(label, passed) {
  if (!isTRUE(passed)) stop("FAILED: ", label)
  results <<- c(results, label)
}
means <- setNames(rep(0, length(GROUPS)), GROUPS)
means[c("anc_TSB", "anc_SN", "stk_TSB", "stk_SN")] <- c(2, 5, 4, 7)
check_code("Negative TSB-vs-SN means higher in SN",
           sum(means * contrast_weights("TSB_vs_SN_at_anc")) == -3)
check_code("Positive stk-vs-anc means higher in stk",
           sum(means * contrast_weights("stk_vs_anc_at_TSB")) == 2)
check_code("The same medium response in both genotypes gives zero interaction",
           sum(means * contrast_weights("interaction_stk_vs_anc_at_TSB_vs_SN")) == 0)
means["stk_SN"] <- 4
check_code("A larger rise in SN in the ancestor gives a positive interaction",
           sum(means * contrast_weights("interaction_stk_vs_anc_at_TSB_vs_SN")) == 3)
means[c("stp_TSB", "stp_SN")] <- c(1, 6)
check_code("The stp interaction uses stp and ancestor in the correct order",
           sum(means * contrast_weights("interaction_stp_vs_anc_at_TSB_vs_SN")) == -2)
check_code("The stp-vs-stk interaction preserves that mutant ordering",
           sum(means * contrast_weights("interaction_stp_vs_stk_at_TSB_vs_SN")) == -5)
unsupported <- tryCatch({ contrast_weights("stk_vs_anc"); FALSE }, error = function(e) TRUE)
check_code("Pooled contrasts are not assigned assumed weights", unsupported)
parsed <- parse_sample_name(c("2anc_pulseSN_1", "JE2stk_TSB_3",
                             "20260612_053_rep_C41881_S1143466_JE2stk_pulseSN_1_Enriched"))
check_code("Shortened pulse names recover the canonical sample ID",
           identical(parsed$sample, c("JE2anc_pulseSN_1", "JE2stk_TSB_3", "JE2stk_pulseSN_1")))
check_code("Medium labels distinguish SNpulse from raw pulseSN spelling",
           identical(as.character(parsed$medium), c("SNpulse", "TSB", "SNpulse")))
check_code("Reinjection flags follow the acquisition filename",
           identical(parsed$reinjected, c(FALSE, FALSE, TRUE)))
output_dir <- file.path(REPO_ROOT, "reports")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c("# Step 1 code checks", "", paste(length(results), "focused checks passed."), "",
             paste0("- PASS: ", results), "",
             "These synthetic checks use no input measurements. The original testthat suite is separate."),
           file.path(output_dir, "step01_code_checks.md"))
message(length(results), "/", length(results), " focused code checks passed.")
