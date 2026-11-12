# Small live API demonstration using public IDs, not experimental measurements.
# API failures stay visible; the local map works independently of these services.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out_dir <- file.path(REPO_ROOT, "results", "step03_api")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
selected <- read_tsv(file.path(REPO_ROOT, "config", "interesting_proteins.tsv"), show_col_types = FALSE)
query <- prepare_uniprot_query(selected$query)
write_tsv(query$excluded, file.path(out_dir, "query_exclusions.tsv"))
if (!length(query$ids)) stop("No unambiguous quality-checked UniProt queries; inspect query_exclusions.tsv.")

# Five public examples are sufficient to verify the API workflow (Stk1/Stp1/TcyABC).
accessions <- head(query$ids, 5)
entries <- setNames(lapply(accessions, uniprot_entry), accessions)
summaries <- bind_rows(lapply(entries, function(x) x$summary))
write_tsv(summaries, file.path(out_dir, "uniprot_summary.tsv"))
crossrefs <- bind_rows(lapply(names(entries), function(accession) {
  refs <- entries[[accession]]$entry$uniProtKBCrossReferences
  if (!length(refs)) return(tibble(uniprot_acc = character(), database = character(), reference_id = character()))
  tibble(uniprot_acc = accession,
         database = vapply(refs, function(x) x$database, character(1)),
         reference_id = vapply(refs, function(x) x$id, character(1)))
}))
write_tsv(crossrefs, file.path(out_dir, "uniprot_crossreferences.tsv"))
if (any(summaries$accession_status != "exact_primary_accession") ||
    n_distinct(summaries$organism_taxon) != 1 || anyNA(summaries$organism_taxon)) {
  stop("Review unresolved/redirected accessions or mixed organisms in uniprot_summary.tsv before querying STRING.")
}
species <- unique(summaries$organism_taxon)
resolved <- string_resolve(accessions, species = species)
resolution_provenance <- api_provenance(resolved)
reference_pairs <- filter(crossrefs, database == "STRING")
resolved <- resolved |>
  mutate(verified_against_UniProt_STRING_crossref = paste(identifier, stringId) %in%
           paste(reference_pairs$uniprot_acc, reference_pairs$reference_id),
         accepted_for_demo = verified_against_UniProt_STRING_crossref &
           !mapping_status %in% c("unresolved", "ambiguous_candidates", "species_mismatch_review_required"))
write_tsv(resolved, file.path(out_dir, "string_resolution.tsv"))
accepted <- unique(resolved$stringId[resolved$accepted_for_demo])
if (!length(accepted)) stop("STRING matches need review; none were verified against UniProt STRING cross-references.")

functional <- string_network(accepted, species, network_type = "functional", required_score = 700)
physical <- string_network(accepted, species, network_type = "physical", required_score = 700)
partners <- string_partners(accepted, species, network_type = "physical", required_score = 700, limit = 10)
annotations <- string_functional_annotation(accepted, species)
write_tsv(functional, file.path(out_dir, "string_functional_selected_edges.tsv"))
write_tsv(physical, file.path(out_dir, "string_physical_selected_edges.tsv"))
write_tsv(partners, file.path(out_dir, "string_public_context_partners.tsv"))
write_tsv(annotations, file.path(out_dir, "string_functional_annotations.tsv"))
write_tsv(annotate_reference_edges(functional, crossrefs), file.path(out_dir, "functional_edges_with_JE2_labels.tsv"))
write_tsv(annotate_reference_edges(physical, crossrefs), file.path(out_dir, "physical_edges_with_JE2_labels.tsv"))
write_tsv(annotate_reference_edges(partners, crossrefs), file.path(out_dir, "public_partners_with_JE2_labels.tsv"))
kegg_ids <- unique(filter(crossrefs, database == "KEGG")$reference_id)
if (!length(kegg_ids)) stop("No KEGG cross-reference returned for the demonstration proteins.")
kegg_pathways <- kegg_links(kegg_ids, target = "pathway")
kegg_example <- kegg_entry(kegg_ids[1])
write_tsv(kegg_pathways, file.path(out_dir, "kegg_gene_pathway_links.tsv"))
writeLines(kegg_example$text, file.path(out_dir, "kegg_example_gene.txt"))
all_provenance <- c(lapply(entries, api_provenance),
                    list(resolution = resolution_provenance, functional = api_provenance(functional),
                         physical = api_provenance(physical), partners = api_provenance(partners),
                         annotations = api_provenance(annotations),
                         kegg_links = api_provenance(kegg_pathways), kegg_entry = api_provenance(kegg_example)))
saveRDS(all_provenance, file.path(out_dir, "api_provenance.rds"))
summary_lines <- apply(select(summaries, requested_accession, protein_name, organism_taxon, accession_status), 1,
                       function(x) paste0("| ", paste(x, collapse = " | "), " |"))
writeLines(c("# Step 3: live public API demonstration", "",
             paste("Run completed:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
             "Requests used only public identifiers for the five selected reference proteins. No experimental measurements or sample identifiers were sent.", "",
             "| Accession | Reference protein name | Taxon | Accession check |",
             "|---|---|---|---|", summary_lines, "",
             paste("STRING organism selected from returned UniProt annotation:", species),
             paste("STRING IDs confirmed by UniProt STRING cross-references:", length(accepted)),
             paste("Selected-protein functional edges:", nrow(functional)),
             paste("Selected-protein physical edges:", nrow(physical)),
             paste("Physical public-context partner rows:", nrow(partners)),
             paste("Functional-annotation rows:", nrow(annotations)), "",
             paste("KEGG gene-to-pathway links:", nrow(kegg_pathways)),
             paste("KEGG example gene entry:", kegg_example$status), "",
             "These API outputs demonstrate reference-resource access, not new differential-expression or enrichment results.",
             "Context partners are not automatically measured or changed in the JE2 experiment; reference mappings do not transfer exact phosphosite coordinates.",
             "Cached response bodies, API URLs, response headers and retrieval times are retained locally in results/api_cache and results/step03_api.", "",
             "Official documentation: [UniProt](https://www.uniprot.org/help/api), [STRING](https://string-db.org/help/api/), [KEGG](https://www.genome.jp/kegg/rest/keggapi.html)."),
           file.path(REPO_ROOT, "reports", "step03_live_api.md"))
message("LIVE API DEMO COMPLETE: ", length(accessions), " UniProt entries, ", length(accepted), " verified STRING IDs.")
