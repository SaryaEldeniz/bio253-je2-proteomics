# Synthetic alias/mapping checks. No network, private data or testthat required.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
passed <- character()
check_id <- function(label, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", label)
  passed <<- c(passed, label)
}
map <- tibble(protein_key = c("A", "B", "C", "D"), locus_tag = paste0("LOC", 1:4),
              protein_id = paste0("ACC", 1:4), display_name = c("Alpha", "Beta", "Gamma", "Contaminant"),
              entry_class = c("cds", "cds", "cds", "contaminant"),
              uniprot_acc = c("Q00001", "Q00001", "Q00002", NA_character_),
              ortholog_quality = c("good", "good", "check", NA_character_),
              eligible_uniprot_annotation = c(TRUE, TRUE, FALSE, FALSE),
              string_id = paste0("ORG.", c("A", "B", "C", "D")))
aliases <- tibble(protein_key = c("A", "B", "C", "D", "A", "B", "A"),
                  alias = c("Alpha", "Beta", "Gamma", "Contaminant", "Q00001", "Q00001", "Example kinase"))
check_id("Unique aliases recover the intended JE2 protein", lookup_je2("Alpha", map = map, aliases = aliases)$protein_key == "A")
check_id("Lookup is case insensitive", lookup_je2("alpha", map = map, aliases = aliases)$protein_key == "A")
ambiguous <- lookup_je2("Q00001", map = map, aliases = aliases)
check_id("Many-to-one UniProt aliases retain both JE2 entries", nrow(ambiguous) == 2 && all(ambiguous$lookup_status == "ambiguous"))
missing <- lookup_je2("absent", map = map, aliases = aliases)
check_id("Unknown aliases remain explicit unresolved rows", nrow(missing) == 1 && missing$lookup_status == "unresolved" && is.na(missing$protein_key))
check_id("Repeated queries retain their query indices", identical(lookup_je2(c("Alpha", "Alpha"), map = map, aliases = aliases)$query_index, 1:2))
check_id("Description searches work as literal substrings", lookup_je2("kinase", exact = FALSE, map = map, aliases = aliases)$protein_key == "A")
prepared <- prepare_uniprot_query(c("Alpha", "Gamma", "Q00001", "absent"), map, aliases)
check_id("Web-query preparation keeps only resolved good mappings", identical(prepared$ids, "Q00001") && nrow(prepared$mapping) == 1 && nrow(prepared$excluded) == 4)
links <- tibble(node_a = c("ORG.A", "ORG.A", "ORG.B"), node_b = c("ORG.B", "ORG.D", "ORG.C"), combined_score = c(850, 950, 500))
edges <- local_interactions("Alpha", map = map, aliases = aliases, links = links)$edges
check_id("Local partner lookup excludes contaminant endpoints", nrow(edges) == 1 && edges$protein_key_b == "B")
within <- local_interactions("Alpha", within_queries = TRUE, map = map, aliases = aliases, links = links)$edges
check_id("Within-query networks do not silently add neighbours", nrow(within) == 0)
rejected <- tryCatch({ local_interactions("Alpha", required_score = -1, map = map, aliases = aliases, links = links); FALSE }, error = function(e) TRUE)
check_id("Invalid network thresholds are rejected", rejected)
refs <- tibble(uniprot_acc = c("Q00001", "Q00002"), database = "STRING", reference_id = c("REF.AB", "REF.C"))
public_edges <- tibble(stringId_A = c("REF.AB", "REF.AB"), stringId_B = c("REF.C", "REF.UNKNOWN"))
map$in_total_DEA <- c(TRUE, TRUE, TRUE, FALSE)
labelled <- annotate_reference_edges(public_edges, refs, map)
check_id("Reference-edge labels preserve edge count and many-to-one mappings",
         nrow(labelled) == 2 && all(labelled$n_JE2_matches_A == 2))
check_id("Unknown public partners remain unmapped rather than called unmeasured",
         is.na(labelled$any_in_total_DEA_B[2]))
dir.create(file.path(REPO_ROOT, "reports"), recursive = TRUE, showWarnings = FALSE)
writeLines(c("# Step 3 identifier code checks", "", paste(length(passed), "checks passed."), "", paste0("- PASS: ", passed)),
           file.path(REPO_ROOT, "reports", "step03_identifier_code_checks.md"))
message(length(passed), "/", length(passed), " identifier checks passed.")
