# Step 3: annotation/alias joins. Step 2 is skipped at the user's request.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out_dir <- file.path(REPO_ROOT, "results", "step03")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- read_je2_metadata()
terms <- read_je2_terms()
annotation <- read_annotation()
total <- read_dea("total", "TSB_vs_SN_at_anc", drop_contaminants = FALSE) |>
  distinct(protein_key, is_contaminant)
phospho <- read_dea("phospho_STYH", "TSB_vs_SN_at_anc", drop_contaminants = FALSE) |>
  distinct(protein_key, is_contaminant)
stopifnot(!anyDuplicated(meta$protein_key), !anyDuplicated(annotation$protein_key))

# The union retains unmatched and excluded entries in the audit instead of dropping them.
keys <- tibble(protein_key = unique(c(meta$protein_key, annotation$protein_key,
                                     total$protein_key, phospho$protein_key)))
map <- keys |>
  left_join(meta, by = "protein_key", relationship = "one-to-one") |>
  left_join(annotation |> transmute(protein_key, course_locus = locus_tag,
                                    course_product = product, course_protein_id = protein_id),
            by = "protein_key", relationship = "one-to-one") |>
  left_join(rename(total, total_DEA_contaminant_flag = is_contaminant), by = "protein_key") |>
  left_join(rename(phospho, phospho_DEA_contaminant_flag = is_contaminant), by = "protein_key") |>
  mutate(in_total_DEA = protein_key %in% total$protein_key,
         in_phospho_DEA = protein_key %in% phospho$protein_key,
         metadata_match = !is.na(entry_class),
         entry_class = coalesce(entry_class,
                                if_else(coalesce(total_DEA_contaminant_flag, FALSE) |
                                          coalesce(phospho_DEA_contaminant_flag, FALSE), "contaminant", "unclassified")),
         locus_tag = coalesce(locus_tag, course_locus),
         product = coalesce(product, course_product),
         protein_id = coalesce(protein_id, course_protein_id),
         eligible_uniprot_annotation = entry_class == "cds" & coalesce(ortholog_quality == "good", FALSE) & !is.na(uniprot_acc))

# Resolve supplied STRING identifiers through exact aliases, not a guessed prefix.
alias_path <- file.path(DATA_ROOT, "raw", "resources", "string_aliases.txt")
string_aliases <- read_tsv(alias_path, show_col_types = FALSE)
names(string_aliases)[1] <- "string_id"
string_matches <- string_aliases |>
  filter(alias %in% map$protein_key) |>
  distinct(protein_key = alias, string_id)
string_counts <- string_matches |> count(protein_key, name = "n_string_matches")
single_matches <- string_matches |> inner_join(filter(string_counts, n_string_matches == 1), by = "protein_key")
term_counts <- terms |> distinct(protein_key, category, term) |> count(protein_key, name = "n_annotation_terms")
physical <- read_local_physical_links()
network_ids <- union(physical$node_a, physical$node_b)
map <- map |>
  left_join(select(single_matches, protein_key, string_id), by = "protein_key", relationship = "one-to-one") |>
  left_join(string_counts, by = "protein_key", relationship = "one-to-one") |>
  left_join(term_counts, by = "protein_key", relationship = "one-to-one") |>
  mutate(n_string_matches = coalesce(n_string_matches, 0L),
         n_annotation_terms = coalesce(n_annotation_terms, 0L),
         has_local_physical_edges = !is.na(string_id) & string_id %in% network_ids) |>
  left_join(KEY_PROTEINS |> select(locus_tag, curated_name = short_name), by = "locus_tag", relationship = "many-to-one") |>
  mutate(display_name = coalesce(curated_name, gene_symbol_best, up_gene_primary, locus_tag, protein_id, protein_key))

# Count target collisions both across the reference and within this total DEA universe.
target_counts <- map |> filter(eligible_uniprot_annotation) |>
  group_by(uniprot_acc) |>
  summarise(n_JE2_good_all = n_distinct(protein_key),
            n_JE2_good_total_DEA = n_distinct(protein_key[in_total_DEA]), .groups = "drop")
map <- map |> left_join(target_counts, by = "uniprot_acc", relationship = "many-to-one") |>
  select(protein_key, entry_class, metadata_match, locus_tag, protein_id, product, aa_length,
         display_name, gene_symbol_best, up_gene_primary, uniprot_acc, up_protein_name,
         ortholog_quality, blastortho_evalue, up_length_ratio, eligible_uniprot_annotation,
         n_JE2_good_all, n_JE2_good_total_DEA, string_id, n_string_matches,
         in_string, n_annotation_terms, has_local_physical_edges, in_total_DEA, in_phospho_DEA,
         total_DEA_contaminant_flag, phospho_DEA_contaminant_flag, up_keywords, up_go_ids,
         up_subcellular, up_kegg)

# Build a searchable alias vocabulary. Names such as Stk1/PknB need curated aliases
# because the provided gene-symbol fields are empty for these regulators.
base_aliases <- map |>
  select(protein_key, locus_tag, protein_id, uniprot_acc, display_name, gene_symbol_best, up_gene_primary,
         product, up_protein_name) |>
  pivot_longer(-protein_key, names_to = "alias_source", values_to = "alias")
provided_aliases <- string_matches |>
  inner_join(string_aliases, by = "string_id", relationship = "many-to-many") |>
  transmute(protein_key, alias, alias_source = paste0("STRING:", source))
synonyms <- tibble(locus_tag = c("B7H15_06355", "B7H15_06355", "B7H15_06350", "B7H15_07015"),
                   alias = c("PknB", "Stk1", "Stp1", "KatA")) |>
  inner_join(select(map, protein_key, locus_tag), by = "locus_tag", relationship = "many-to-one") |>
  transmute(protein_key, alias, alias_source = "curated course/literature name")
lookup_aliases <- bind_rows(base_aliases, provided_aliases, synonyms,
                            map |> transmute(protein_key, alias = protein_key, alias_source = "JE2 protein key"),
                            map |> transmute(protein_key, alias = string_id, alias_source = "local STRING ID")) |>
  filter(!is.na(alias), nzchar(alias)) |> distinct(protein_key, alias, alias_source)

checks <- tibble(
  check = c("One master row per protein key", "All total DEA keys retained", "All phospho DEA parent keys retained",
            "Primary annotation mappings meet stated E-value and length policy", "No unhandled local STRING ambiguity"),
  passed = c(nrow(map) == n_distinct(map$protein_key),
             all(total$protein_key %in% map$protein_key), all(phospho$protein_key %in% map$protein_key),
             all(map$blastortho_evalue[map$eligible_uniprot_annotation] <= 1e-5 &
                   between(map$up_length_ratio[map$eligible_uniprot_annotation], 0.9, 1.1)),
             all(map$n_string_matches <= 1)))
write_tsv(checks, file.path(out_dir, "mapping_checks.tsv"))
if (!all(checks$passed)) stop("Mapping checks failed. Inspect results/step03/mapping_checks.tsv.")
write_tsv(map, file.path(out_dir, "id_map.tsv"))
write_tsv(lookup_aliases, file.path(out_dir, "lookup_aliases.tsv.gz"))
write_tsv(string_matches, file.path(out_dir, "local_string_matches.tsv"))
write_tsv(terms |> semi_join(filter(map, entry_class == "cds"), by = "protein_key"),
          file.path(out_dir, "terms_cds.tsv.gz"))
write_tsv(filter(map, !eligible_uniprot_annotation | coalesce(n_JE2_good_all > 1, FALSE)),
          file.path(out_dir, "excluded_weak_or_many_to_one.tsv"))

coverage <- bind_rows(lapply(c("total", "phospho"), function(layer) {
  included <- if (layer == "total") map$in_total_DEA else map$in_phospho_DEA
  x <- map[included, ]
  tibble(layer, tested_parent_entries = nrow(x), ordinary_CDS = sum(x$entry_class == "cds"),
         good_mapped_JE2 = sum(x$eligible_uniprot_annotation),
         unique_good_UniProt = n_distinct(x$uniprot_acc[x$eligible_uniprot_annotation]),
         ordinary_CDS_STRING_alias = sum(x$entry_class == "cds" & !is.na(x$string_id)),
         ordinary_CDS_with_terms = sum(x$entry_class == "cds" & x$n_annotation_terms > 0),
         ordinary_CDS_with_physical_edges = sum(x$entry_class == "cds" & x$has_local_physical_edges))
}))
write_tsv(coverage, file.path(out_dir, "mapping_coverage.tsv"))
audit <- filter(map, in_total_DEA) |>
  count(entry_class, metadata_match, total_DEA_contaminant_flag, name = "entries")
write_tsv(audit, file.path(out_dir, "total_DEA_classification.tsv"))
selected <- read_tsv(file.path(REPO_ROOT, "config", "interesting_proteins.tsv"), show_col_types = FALSE)
examples <- lookup_je2(selected$query, map = map, aliases = lookup_aliases) |>
  left_join(selected, by = "query", relationship = "many-to-one")
write_tsv(examples, file.path(out_dir, "interesting_protein_lookup.tsv"))
write_tsv(protein_resource_links(selected$query, map, lookup_aliases), file.path(out_dir, "resource_links.tsv"))
example_edges <- local_interactions(selected$query, map = map, aliases = lookup_aliases, links = physical)
write_tsv(example_edges$edges, file.path(out_dir, "example_local_physical_partners.tsv"))

# Save a compact, shareable report; the full identifier tables remain local.
coverage_lines <- apply(coverage, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
example_lines <- apply(select(examples, label, locus_tag, uniprot_acc, ortholog_quality), 1,
                       function(x) paste0("| ", paste(x, collapse = " | "), " |"))
writeLines(c("# Step 3: identifier map and annotation", "",
             paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
             "Step 2 was skipped at the user's request because the data were already validated.", "",
             sprintf("Master map: %s entries; %s searchable alias records; %s/%s mapping checks passed.",
                     nrow(map), nrow(lookup_aliases), sum(checks$passed), nrow(checks)), "",
             "## Coverage of the current DEA feature universe", "",
             "Phospho counts below are unique parent proteins, not phosphosite counts. Both universes use the current ancestor TSB-vs-SN contrast.", "",
             "| Layer | Tested parent entries | CDS | Good JE2 maps | Unique UniProt | CDS STRING aliases | CDS with terms | CDS with physical edges |",
             "|---|---:|---:|---:|---:|---:|---:|---:|", coverage_lines, "",
             "STRING alias resolution, functional-term coverage and physical-edge coverage are separate quantities.", "",
             "## Classification correction", "",
             paste("Metadata identifies", sum(map$in_total_DEA & map$entry_class == "contaminant" & !coalesce(map$total_DEA_contaminant_flag, TRUE)),
                   "total-DEA contaminant entries not flagged by its is_contaminant column."),
             "The primary ordinary-protein analysis uses entry_class == cds. Unmatched, weak and excluded rows remain in the local audit.", "",
             "## Example public identifiers", "",
             "These are user-editable starting candidates, not newly selected significant hits.", "",
             "| Label | JE2 locus | UniProt reference ortholog | Mapping quality |",
             "|---|---|---|---|", example_lines, "",
             "## Use and limits", "",
             "UniProt assignments are reference-proteome ortholog mappings to S. aureus NCTC8325. Good means E-value <= 1e-5 and length ratio 0.9-1.1.",
             "Keep JE2 sequence coordinates for phosphosites. Preserve many-to-one mappings; ORA uses unique target IDs and ranked analysis needs an explicit collision rule.",
             "Local physical-network partners are associations from the supplied export, not evidence of direct kinase substrates.", "",
             "See docs/step03_identifiers.md and docs/web_resources.md for interactive R use and official API references.", "",
             "## Input provenance", "",
             paste0("- Mapping zip MD5: `", unname(tools::md5sum(annotation_zip_path())), "`"),
             paste0("- Local STRING aliases MD5: `", unname(tools::md5sum(alias_path)), "`")),
           file.path(REPO_ROOT, "reports", "step03_identifiers.md"))
writeLines(capture.output(sessionInfo()), file.path(REPO_ROOT, "reports", "step03_session.txt"))
message("STEP 3 LOCAL MAP COMPLETE")
print(coverage, width = Inf)
print(audit, width = Inf)
