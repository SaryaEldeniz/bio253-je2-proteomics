# Network-free synthetic checks. fgsea/BiocParallel are required for real GSEA.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if (!requireNamespace("fgsea", quietly = TRUE) || !requireNamespace("BiocParallel", quietly = TRUE)) {
  stop("GSEA checks require fgsea and BiocParallel. Install the declared Step 5 dependencies; do not report these checks as passed.")
}
passed <- character()
check_gsea <- function(label, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", label)
  passed <<- c(passed, label)
}
rejected <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
ranks <- setNames(seq(5, -5, length.out = 200), sprintf("ID%03d", 1:200))
sets <- list(top = names(ranks)[1:15], bottom = names(ranks)[186:200],
             mixed = names(ranks)[seq(4, 194, by = 10)],
             too_small = names(ranks)[1:3])
set.seed(79)
previous_rng <- .Random.seed
a <- run_step05_gsea(ranks, sets, seed = 42)
check_gsea("Caller RNG state is preserved", identical(previous_rng, .Random.seed))
b <- run_step05_gsea(ranks, sets, seed = 42)
# Compare every result column, not elapsed-time metadata attached to the table.
check_gsea("The same seed gives reproducible estimates and leading edges",
           identical(lapply(a, identity), lapply(b, identity)) &&
             identical(attr(a, "gsea_method"), attr(b, "gsea_method")))
check_gsea("Effective size filtering excludes the small set", setequal(a$pathway, c("top", "bottom", "mixed")))
check_gsea("Coherent top and bottom sets have the expected enrichment signs",
           a$NES[a$pathway == "top"] > 0 && a$NES[a$pathway == "bottom"] < 0)
check_gsea("The coherent synthetic signals have strong raw support",
           all(a$pval[a$pathway %in% c("top", "bottom")] < .01))
check_gsea("Leading edges contain identifiers from their pathway",
           all(vapply(seq_len(nrow(a)), function(i) all(a$leadingEdge[[i]] %in% sets[[a$pathway[i]]]), logical(1))))
check_gsea("Raw engine adjusted values and uncertainty/status are retained",
           all(c("padj_engine", "log2err", "status", "pval_at_floor") %in% names(a)))
bad <- ranks; names(bad)[2] <- names(bad)[1]
check_gsea("Duplicate target identifiers are rejected", rejected(run_step05_gsea(bad, sets)))
bad <- ranks; bad[1] <- NA_real_
check_gsea("Nonfinite rank values are rejected", rejected(run_step05_gsea(bad, sets)))
empty <- run_step05_gsea(ranks, list(tiny = names(ranks)[1:3]), seed = 42)
check_gsea("No eligible pathways returns an explicit empty result", nrow(empty) == 0 && attr(empty, "gsea_method")$n_pathways_eligible == 0)
method <- attr(a, "gsea_method")
check_gsea("Method records serial multilevel implementation and size/seed settings",
           method$engine == "fgsea::fgseaMultilevel" && method$backend == "BiocParallel::SerialParam" && method$seed == 42 && method$min_size == 10)
dir.create(file.path(REPO_ROOT, "reports"), recursive = TRUE, showWarnings = FALSE)
writeLines(c("# Step 5 GSEA engine checks", "", paste(length(passed), "checks passed."), "",
             paste0("- PASS: ", passed), "", "Synthetic ranked sets only; no unpublished measurements or network requests.",
             paste0("fgsea version: ", method$fgsea_version)),
           file.path(REPO_ROOT, "reports", "step05_gsea_code_checks.md"))
message(length(passed), "/", length(passed), " GSEA engine checks passed.")
