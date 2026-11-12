# Offline checks: no unpublished data, network calls or testthat required.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
passed <- character()
check_api <- function(name, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", name)
  passed <<- c(passed, name)
}
expect_rejected <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
check_api("measurement tables cannot be sent as identifiers",
          expect_rejected(.api_identifiers(data.frame(id = "Q2FZ64", abundance = 10))))
check_api("empty and NA identifiers are rejected",
          expect_rejected(.api_identifiers(c("Q2FZ64", NA_character_))))
check_api("organism is explicit", expect_rejected(.api_species(character())))
raw <- tibble::tibble(queryIndex = c("0", "1", "1", "2"),
                       stringId = c("123.A", "123.B", "123.C", "999.D"),
                       ncbiTaxonId = c("123", "123", "123", "999"))
mapping <- .string_resolve_rows(c("123.A", "Q2FZ64", "Q2FZ65", "missing", "123.A"), "123", raw)
check_api("unresolved inputs are retained", any(mapping$identifier == "missing" & mapping$mapping_status == "unresolved"))
check_api("multiple returned candidates are labelled ambiguous",
          all(mapping$mapping_status[mapping$identifier == "Q2FZ64"] == "ambiguous_candidates"))
check_api("wrong-species mappings require review",
          mapping$mapping_status[mapping$identifier == "Q2FZ65"] == "species_mismatch_review_required")
check_api("duplicate input rows remain traceable", sum(mapping$identifier == "123.A") == 2L)
single <- .string_resolve_rows("Q2FZ64", "123", raw[1, ])
check_api("service best match is not called exact", single$mapping_status == "service_best_match_review_required")
check_api("mixed STRING organism namespaces are rejected",
          expect_rejected(.string_validate_ids(c("123.A", "999.B"), "123")))
fake_entry <- list(primaryAccession = "Q2FZ64", uniProtkbId = "KINASE_STAA8",
                    proteinDescription = list(recommendedName = list(fullName = list(value = "Example kinase"))),
                    organism = list(scientificName = "Example organism", taxonId = 123),
                    genes = list(list(geneName = list(value = "stk1"))), sequence = list(length = 100))
summary <- .uniprot_summary(fake_entry, "Q2FZ64")
check_api("UniProt exact accession and organism are retained", summary$accession_status == "exact_primary_accession" &&
            summary$organism_taxon == "123" && summary$sequence_length == 100L)
check_api("redirected accession is distinguished", .uniprot_summary(fake_entry, "Q2FZ65")$accession_status == "secondary_or_redirected_accession")
check_api("unresolved UniProt lookup remains explicit", .uniprot_summary(list(), "Q2FZ64")$accession_status == "unresolved")
check_api("KEGG requires explicit entries rather than a whole organism/database",
          expect_rejected(.kegg_identifiers("hsa")))
check_api("KEGG accepts organism genes and pathway entry identifiers",
          identical(.kegg_identifiers(c("example:GENE_1", "map00010", "path:map00010")),
                    c("example:GENE_1", "map00010", "path:map00010")))
kegg_fixture <- .kegg_parse_links("example:GENE_1\tpath:example00010\nexample:GENE_1\tpath:example00020\n")
check_api("KEGG headerless links retain the first data row and both columns",
          nrow(kegg_fixture) == 2L && kegg_fixture$source_id[1] == "example:GENE_1" &&
            kegg_fixture$target_id[1] == "path:example00010")
check_api("KEGG no-link response returns a typed empty table",
          nrow(.kegg_parse_links("")) == 0L && identical(names(.kegg_parse_links("")), c("source_id", "target_id")))
check_api("Malformed KEGG responses are not silently parsed",
          expect_rejected(.kegg_parse_links("<html>Service unavailable</html>")))
output_dir <- file.path(REPO_ROOT, "reports")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c("# Web-resource client offline checks", "", paste(length(passed), "checks passed."), "",
             paste0("- PASS: ", passed), "", "No live API requests were made by these checks."),
           file.path(output_dir, "step03_api_code_checks.md"))
message(length(passed), "/", length(passed), " web-resource offline checks passed.")
