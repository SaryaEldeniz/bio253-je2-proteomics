# Standalone, read-only Step 5 result checks using base R only.
# Run from Sarya_George: Rscript --vanilla tests/check_step05_results.R
# A first argument may instead name the absolute Sarya_George project folder.
# Before the enrichment run has finished, this reports missing final outputs and
# exits successfully without polling or treating an incomplete run as a failure.

arguments <- commandArgs(trailingOnly = TRUE)
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[1]) else ""
candidates <- c(arguments[1], if (nzchar(script_path)) dirname(dirname(script_path)),
                getwd(), file.path(getwd(), "Sarya_George"))
candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
roots <- candidates[file.exists(file.path(candidates, "Sarya_George.Rproj"))]
if (!length(roots)) stop("Cannot locate Sarya_George.Rproj; supply its project folder as the first argument.")
project <- normalizePath(roots[1], winslash = "/", mustWork = TRUE)
output <- file.path(project, "results", "step05")
report <- file.path(project, "reports", "step05")
required <- c(file.path(output, c("ora_all.tsv.gz", "gsea_all.tsv.gz")),
              file.path(report, c("comparison_plan.tsv", "universe_and_hit_counts.tsv",
                                  "annotation_coverage.tsv", "gsea_run_log.tsv")))
if (!all(file.exists(required))) {
  cat("Step 5 final outputs are not yet complete; result validation was not run.\n")
  cat("Missing:", paste(basename(required[!file.exists(required)]), collapse = ", "), "\n")
  quit(save = "no", status = 0L)
}

read_table <- function(path) {
  connection <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(connection))
  read.delim(connection, header = TRUE, sep = "\t", quote = "\"", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
}
numeric_equal <- function(x, y, tolerance = 1e-10) {
  if (length(x) != length(y) || !identical(is.na(x), is.na(y))) return(FALSE)
  keep <- !is.na(x)
  if (!any(keep)) return(TRUE)
  x <- x[keep]; y <- y[keep]
  all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= tolerance * pmax(abs(x), abs(y), .Machine$double.xmin))
}
split_ids <- function(value) {
  if (length(value) != 1L || is.na(value) || !nzchar(value)) return(character())
  strsplit(value, ";", fixed = TRUE)[[1]]
}
key <- function(policy, contrast) paste(policy, contrast, sep = "__")
expected_direction <- function(contrast, effect_sign) {
  positive <- effect_sign == "positive"
  if (grepl("^TSB_vs_SN_at_(anc|stk|stp)$", contrast)) {
    return(if (positive) "Higher in TSB" else "Higher in SN")
  }
  if (grepl("^(stk|stp)_vs_anc_at_(TSB|SN)$", contrast)) {
    mutant <- sub("_vs_anc_at_.*$", "", contrast)
    return(if (positive) paste("Higher in", mutant) else "Higher in anc")
  }
  if (grepl("^interaction_(stk|stp)_vs_(anc|stk)_at_TSB_vs_SN$", contrast)) {
    pair <- strsplit(gsub("^interaction_|_at_TSB_vs_SN$", "", contrast), "_vs_", fixed = TRUE)[[1]]
    # g1-vs-g2 at TSB-vs-SN equals R_g2 - R_g1, R_g = SN_g - TSB_g.
    return(paste("SN-TSB change more positive in", if (positive) pair[2] else pair[1]))
  }
  stop("Unexpected comparison in exported output: ", contrast)
}

ora <- read_table(required[1]); gsea <- read_table(required[2])
plan <- read_table(required[3]); audit <- read_table(required[4])
coverage <- read_table(required[5]); run_log <- read_table(required[6])
policies <- c("primary", "replicated", "observed_replicated", "local_JE2")
collections <- c("KEGG", "GO_BP", "GO_MF", "GO_CC", "UniProt_keywords")
expected_runs <- as.vector(outer(policies, plan$contrast, key))
stopifnot(nrow(plan) == 10L, !anyDuplicated(plan$contrast), length(expected_runs) == 40L,
          nrow(audit) == 40L, !anyDuplicated(key(audit$policy, audit$contrast)),
          setequal(key(audit$policy, audit$contrast), expected_runs),
          nrow(run_log) == 40L, !anyDuplicated(key(run_log$policy, run_log$contrast)),
          setequal(key(run_log$policy, run_log$contrast), expected_runs),
          nrow(coverage) == 200L,
          !anyDuplicated(paste(coverage$policy, coverage$contrast, coverage$collection)),
          setequal(coverage$collection, collections),
          all(coverage$annotated_proteins <= coverage$background_size),
          all(coverage$coverage_fraction >= 0 & coverage$coverage_fraction <= 1),
          setequal(unique(key(ora$policy, ora$contrast)), expected_runs),
          setequal(unique(key(gsea$policy, gsea$contrast)), expected_runs),
          !anyDuplicated(paste(ora$policy, ora$contrast, ora$pathway, ora$effect_sign)),
          !anyDuplicated(paste(gsea$policy, gsea$contrast, gsea$pathway)),
          setequal(ora$effect_sign, c("positive", "negative")),
          all(ora$collection %in% collections), all(gsea$collection %in% collections),
          all(ora$p_value >= 0 & ora$p_value <= 1),
          all(is.na(gsea$pval) | (gsea$pval >= 0 & gsea$pval <= 1)),
          all(ora$term_size >= 10 & ora$term_size <= 500),
          all(gsea$size >= 10 & gsea$size <= 500),
          all(ora$term_size < ora$background_size),
          all(gsea$size < gsea$background_size),
          all(ora$overlap >= 0 & ora$overlap <= ora$foreground_size & ora$overlap <= ora$term_size))

check_adjustment <- function(table, p_column) {
  p <- table[[p_column]]
  for (policy in policies) {
    family <- which(table$policy == policy)
    finite <- family[is.finite(p[family])]
    expected <- p.adjust(p[finite], method = "BH", n = length(family))
    stopifnot(numeric_equal(table$FDR_study[finite], expected),
              all(is.na(table$FDR_study[setdiff(family, finite)])))
    for (contrast in plan$contrast) {
      local <- family[table$contrast[family] == contrast]
      finite <- local[is.finite(p[local])]
      expected <- p.adjust(p[finite], method = "BH", n = length(local))
      stopifnot(numeric_equal(table$FDR_local[finite], expected),
                all(is.na(table$FDR_local[setdiff(local, finite)])))
    }
  }
}
check_adjustment(ora, "p_value"); check_adjustment(gsea, "pval")
cat("Passed: 40-run coverage, unique tests, effective set sizes and independently recomputed BH families.\n")

cache_by_run <- list(); ranks_by_run <- list()
for (i in seq_len(nrow(audit))) {
  policy <- audit$policy[i]; contrast <- audit$contrast[i]; run_key <- key(policy, contrast)
  cache_file <- file.path(output, "cache", paste0(run_key, ".rds"))
  rank_file <- file.path(output, "contrasts", paste0(run_key, "__rank_audit.tsv.gz"))
  stopifnot(file.exists(cache_file), file.exists(rank_file))
  cache <- readRDS(cache_file); ranks <- read_table(rank_file)
  cache_by_run[[run_key]] <- cache; ranks_by_run[[run_key]] <- ranks
  stopifnot(nrow(ranks) == audit$universe[i], !anyDuplicated(ranks$target_id),
            !anyDuplicated(ranks$protein_key), all(is.finite(ranks$statistic)),
            all(is.finite(ranks$diff)), all(is.finite(ranks$FDR)),
            setequal(names(cache$ranks), ranks$target_id),
            numeric_equal(unname(cache$ranks[ranks$target_id]), ranks$statistic),
            all(vapply(cache$pathways, function(ids) all(ids %in% ranks$target_id), logical(1))))
  if (policy != "local_JE2") stopifnot(identical(ranks$target_id, ranks$uniprot_acc))
  if (policy %in% c("replicated", "observed_replicated")) stopifnot(all(ranks$all_relevant_groups_2of3))
  if (policy == "observed_replicated") stopifnot(all(ranks$estimate_type == "observed"))
  if (policy == "local_JE2") stopifnot(identical(ranks$target_id, ranks$protein_key))
  foreground <- list(positive = ranks$target_id[ranks$FDR < .1 & ranks$diff > 1],
                     negative = ranks$target_id[ranks$FDR < .1 & ranks$diff < -1])
  stopifnot(length(foreground$positive) == audit$positive_hits[i],
            length(foreground$negative) == audit$negative_hits[i])
  o <- ora[ora$policy == policy & ora$contrast == contrast, , drop = FALSE]
  g <- gsea[gsea$policy == policy & gsea$contrast == contrast, , drop = FALSE]
  stopifnot(nrow(o) == 2L * nrow(g), nrow(g) == audit$tested_terms[i],
            all(o$background_size == nrow(ranks)), all(g$background_size == nrow(ranks)),
            all(g$seed == cache$seed),
            sum(!is.finite(g$pval)) == audit$failed_GSEA[i],
            all(c("positive", "negative") %in% o$effect_sign))
  for (sign in c("positive", "negative")) {
    selected <- o[o$effect_sign == sign, , drop = FALSE]
    expected <- expected_direction(contrast, sign)
    stopifnot(all(selected$direction == expected),
              all(selected$foreground_size == length(foreground[[sign]])))
  }
  for (j in seq_len(nrow(g))) {
    ids <- split_ids(g$leading_edge_ids[j])
    eligible_set <- intersect(cache$pathways[[g$pathway[j]]], ranks$target_id)
    stopifnot(length(eligible_set) == g$size[j], !anyDuplicated(ids),
              all(ids %in% eligible_set))
    cached_row <- match(g$pathway[j], cache$gsea$pathway)
    stopifnot(!is.na(cached_row), setequal(ids, cache$gsea$leadingEdge[[cached_row]]),
              numeric_equal(g$pval[j], cache$gsea$pval[cached_row]),
              numeric_equal(g$NES[j], cache$gsea$NES[cached_row]))
    if (is.finite(g$NES[j]) && g$NES[j] != 0) {
      sign <- if (g$NES[j] > 0) "positive" else "negative"
      stopifnot(g$effect_sign[j] == sign, g$direction[j] == expected_direction(contrast, sign))
    } else stopifnot(g$effect_sign[j] == "unavailable", g$direction[j] == "Unavailable")
  }
}
cat("Passed: all cached ranks, foreground counts, policy restrictions, native directions and leading-edge membership.\n")

# Independent real-data ORA cross-check uses a Fisher exact test, not the
# hypergeometric function used by the implementation. Prefer a primary discovery.
candidate <- which(ora$policy == "primary" & is.finite(ora$FDR_study) & ora$FDR_study < .1)
if (!length(candidate)) candidate <- which(is.finite(ora$p_value))
stopifnot(length(candidate) > 0L)
selected <- candidate[which.min(ora$p_value[candidate])]
row <- ora[selected, , drop = FALSE]
run_key <- key(row$policy, row$contrast)
ranks <- ranks_by_run[[run_key]]; cache <- cache_by_run[[run_key]]
foreground <- ranks$target_id[ranks$FDR < .1 &
  if (row$effect_sign == "positive") ranks$diff > 1 else ranks$diff < -1]
term <- intersect(cache$pathways[[row$pathway]], ranks$target_id)
overlap <- intersect(foreground, term)
N <- nrow(ranks); n <- length(foreground); K <- length(term); k <- length(overlap)
contingency <- matrix(c(k, n-k, K-k, N-K-n+k), nrow = 2L, byrow = TRUE,
                       dimnames = list(c("foreground", "other_tested"), c("in_term", "outside_term")))
stopifnot(all(contingency >= 0), all(row$background_size == N), row$foreground_size == n,
          row$term_size == K, row$overlap == k, setequal(split_ids(row$overlap_ids), overlap),
          numeric_equal(row$p_value, fisher.test(contingency, alternative = "greater")$p.value),
          numeric_equal(row$expected_overlap, n*K/N))
if (n > 0L) stopifnot(numeric_equal(row$fold_enrichment, (k/n)/(K/N)))
cat("Passed: independent Fisher exact test for", row$policy, row$contrast, row$pathway,
    row$effect_sign, "(overlap", k, "of", n, "; study FDR", format(row$FDR_study, digits = 5), ").\n")

members_file <- file.path(output, "contributing_proteins.tsv.gz")
if (file.exists(members_file)) {
  members <- read_table(members_file)
  stopifnot(!anyNA(members$protein_key), !anyNA(members$uniprot_acc),
            all(members$ids == members$uniprot_acc), all(members$FDR_study < .1))
  cat("Passed: contributing proteins have mapped JE2 keys and UniProt accessions.\n")
}
cat("All completed Step 5 result checks passed. No files were changed.\n")
