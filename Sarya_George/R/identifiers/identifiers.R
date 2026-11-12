# Project I: local identifier lookup. No web requests happen in these functions.

annotation_zip_path <- function() {
  path <- Sys.getenv("JE2_ANNOTATION_ZIP", unset = "")
  config <- file.path(REPO_ROOT, "config", "annotation_zip.txt")
  if (!nzchar(path) && file.exists(config)) path <- trimws(readLines(config, n = 1, warn = FALSE))
  if (!nzchar(path) || !file.exists(path)) {
    stop("Set JE2_ANNOTATION_ZIP or config/annotation_zip.txt. See docs/step03_identifiers.md.")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_je2_metadata <- function() {
  read_tsv(unz(annotation_zip_path(), "JE2_protein_metadata.tsv"), show_col_types = FALSE)
}

read_je2_terms <- function() {
  read_tsv(unz(annotation_zip_path(), "JE2_protein_terms.tsv"), show_col_types = FALSE)
}

read_identifier_map <- function() {
  path <- file.path(REPO_ROOT, "results", "step03", "id_map.tsv")
  if (!file.exists(path)) stop('First run source("scripts/02_identifiers/03_build_identifier_map.R").')
  read_tsv(path, show_col_types = FALSE)
}

read_identifier_aliases <- function() {
  path <- file.path(REPO_ROOT, "results", "step03", "lookup_aliases.tsv.gz")
  if (!file.exists(path)) stop('First run source("scripts/02_identifiers/03_build_identifier_map.R").')
  read_tsv(path, show_col_types = FALSE)
}

# Ambiguous aliases produce multiple rows, and unmatched queries remain visible.
# Pass map/aliases explicitly to use this helper in tests without private data.
lookup_je2 <- function(queries, exact = TRUE, map = read_identifier_map(),
                       aliases = read_identifier_aliases()) {
  if (!is.character(queries) || !length(queries) || anyNA(queries) || any(!nzchar(trimws(queries)))) {
    stop("Provide one or more nonempty protein identifiers or names.")
  }
  bind_rows(lapply(seq_along(queries), function(i) {
    query <- trimws(queries[i])
    keep <- if (exact) toupper(aliases$alias) == toupper(query) else {
      str_detect(toupper(aliases$alias), fixed(toupper(query)))
    }
    keys <- unique(aliases$protein_key[!is.na(keep) & keep])
    hits <- map |> filter(protein_key %in% keys)
    if (!nrow(hits)) {
      hits <- map[NA_integer_, , drop = FALSE]
      status <- "unresolved"
    } else status <- if (nrow(hits) == 1) "matched" else "ambiguous"
    mutate(hits, query = query, query_index = i, lookup_status = status, .before = 1)
  }))
}

# Explicitly return exclusions alongside the public IDs suitable for a web query.
prepare_uniprot_query <- function(queries, map = read_identifier_map(), aliases = read_identifier_aliases()) {
  hits <- lookup_je2(queries, map = map, aliases = aliases) |>
    mutate(query_eligible = lookup_status == "matched" & coalesce(eligible_uniprot_annotation, FALSE))
  list(ids = unique(hits$uniprot_acc[hits$query_eligible]),
       mapping = filter(hits, query_eligible), excluded = filter(hits, !query_eligible))
}

protein_resource_links <- function(queries, map = read_identifier_map(), aliases = read_identifier_aliases()) {
  lookup_je2(queries, map = map, aliases = aliases) |>
    transmute(query, lookup_status, protein_key, locus_tag, uniprot_acc, ortholog_quality,
              uniprot_url = if_else(!is.na(uniprot_acc), paste0("https://www.uniprot.org/uniprotkb/", uniprot_acc), NA_character_),
              ncbi_protein_url = if_else(!is.na(protein_id), paste0("https://www.ncbi.nlm.nih.gov/protein/", protein_id), NA_character_))
}

# The local export is an undirected PHYSICAL network. Keep each edge only once.
read_local_physical_links <- function() {
  path <- file.path(DATA_ROOT, "raw", "resources", "string_physical_links.txt")
  if (!file.exists(path)) stop("Missing local STRING physical-links export: ", path)
  read_table(path, col_types = cols(protein1 = col_character(), protein2 = col_character(),
                                    combined_score = col_double()), progress = FALSE) |>
    transmute(node_a = pmin(protein1, protein2), node_b = pmax(protein1, protein2), combined_score) |>
    group_by(node_a, node_b) |>
    summarise(combined_score = max(combined_score), .groups = "drop")
}

local_interactions <- function(queries, required_score = 700, within_queries = FALSE,
                               map = read_identifier_map(), aliases = read_identifier_aliases(),
                               links = read_local_physical_links()) {
  if (length(required_score) != 1 || !is.finite(required_score) || required_score < 0 || required_score > 1000) {
    stop("required_score must be between 0 and 1000.")
  }
  selected <- lookup_je2(queries, map = map, aliases = aliases)
  selected_ids <- unique(selected$string_id[selected$lookup_status == "matched" & !is.na(selected$string_id)])
  nodes <- map |> filter(entry_class == "cds", !is.na(string_id)) |>
    select(string_id, protein_key, locus_tag, display_name, uniprot_acc)
  edges <- links |> filter(combined_score >= required_score,
                           node_a %in% nodes$string_id, node_b %in% nodes$string_id)
  edges <- if (within_queries) filter(edges, node_a %in% selected_ids, node_b %in% selected_ids) else {
    filter(edges, node_a %in% selected_ids | node_b %in% selected_ids)
  }
  edges <- edges |>
    left_join(rename_with(nodes, ~paste0(.x, "_a")), by = c("node_a" = "string_id_a"), relationship = "many-to-one") |>
    left_join(rename_with(nodes, ~paste0(.x, "_b")), by = c("node_b" = "string_id_b"), relationship = "many-to-one") |>
    arrange(desc(combined_score))
  list(queries = selected, edges = edges, network_type = "local physical STRING export")
}

# Add JE2 labels to online reference-organism edges without multiplying edge rows.
# Multiple JE2 orthologs remain visible as counts and semicolon-separated labels.
# Only supplied cross-references are mapped; an unknown partner stays unknown.
annotate_reference_edges <- function(edges, crossrefs, map = read_identifier_map()) {
  if (!all(c("stringId_A", "stringId_B") %in% names(edges))) {
    if (!nrow(edges)) return(edges)
    stop("Expected STRING edge columns stringId_A and stringId_B.")
  }
  collapse_ids <- function(x) {
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))
    if (length(x)) paste(x, collapse = ";") else NA_character_
  }
  nodes <- crossrefs |> filter(database == "STRING") |>
    distinct(uniprot_acc, reference_id) |>
    inner_join(filter(map, eligible_uniprot_annotation), by = "uniprot_acc", relationship = "many-to-many") |>
    group_by(reference_id) |>
    summarise(JE2_protein_keys = collapse_ids(protein_key), JE2_locus_tags = collapse_ids(locus_tag),
              JE2_names = collapse_ids(display_name), n_JE2_matches = n_distinct(protein_key),
              any_in_total_DEA = any(in_total_DEA), .groups = "drop")
  provenance <- attr(edges, "api_provenance")
  result <- edges |>
    left_join(rename_with(nodes, ~paste0(.x, "_A")), by = c("stringId_A" = "reference_id_A"), relationship = "many-to-one") |>
    left_join(rename_with(nodes, ~paste0(.x, "_B")), by = c("stringId_B" = "reference_id_B"), relationship = "many-to-one")
  attr(result, "api_provenance") <- provenance
  result
}
