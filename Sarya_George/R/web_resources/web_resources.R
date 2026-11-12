# Public-resource clients for small, explicit identifier queries.
# Namespaced dependencies: curl, jsonlite and digest (already available locally).
# Public responses and provenance are cached under ignored results/api_cache/.

STRING_API_BASE <- "https://version-12-5.string-db.org/api"
.api_last_request <- new.env(parent = emptyenv())

.api_identifiers <- function(ids, maximum = 50L) {
  if (!is.character(ids) || !length(ids) || length(ids) > maximum ||
      anyNA(ids) || any(!nzchar(ids)) || any(!grepl("^[A-Za-z0-9_.:|+-]+$", ids))) {
    stop("Supply 1-", maximum, " explicit public identifiers as a character vector; no tables or measurements.")
  }
  ids
}

.api_species <- function(species) {
  species <- as.character(species)
  if (length(species) != 1 || is.na(species) ||
      !grepl("^([0-9]+|STRG[A-Za-z0-9]+)$", species)) {
    stop("Supply one explicit STRING taxon/organism identifier.")
  }
  species
}

.api_throttle <- function(service) {
  # One request per second also leaves a conservative margin for UniProt.
  now <- as.numeric(Sys.time())
  previous <- .api_last_request[[service]]
  if (!is.null(previous)) Sys.sleep(max(0, 1.05 - (now - previous)))
  .api_last_request[[service]] <- as.numeric(Sys.time())
}

.api_request <- function(url, params = list(), service, refresh = FALSE,
                         allowed_status = 200L) {
  for (package in c("curl", "jsonlite", "digest")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install package '", package, "' before making API requests.")
    }
  }
  cache_dir <- file.path(REPO_ROOT, "results", "api_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key <- digest::digest(list(url = url, params = params[sort(names(params))]), algo = "sha256")
  cache_file <- file.path(cache_dir, paste0(service, "_", key, ".rds"))
  if (!refresh && file.exists(cache_file)) {
    result <- readRDS(cache_file)
    result$provenance$cache_hit <- TRUE
    return(result)
  }
  # Build a GET query from scalar public parameters; no body/table is uploaded.
  encode <- function(value) utils::URLencode(as.character(value), reserved = TRUE)
  query <- if (length(params)) paste(paste0(encode(names(params)), "=",
                                           vapply(params, encode, character(1))), collapse = "&") else ""
  request_url <- paste0(url, if (nzchar(query)) paste0("?", query) else "")
  if (nchar(request_url) > 7500) stop("Query too long; use a smaller explicit identifier list.")
  handle <- curl::new_handle(useragent = "JE2-Project-I-course-analysis/1.0",
                             timeout = 15, connecttimeout = 10, followlocation = TRUE,
                             accept_encoding = "")
  for (attempt in 1:3) {
    .api_throttle(service)
    response <- tryCatch(curl::curl_fetch_memory(request_url, handle = handle), error = identity)
    if (inherits(response, "error")) {
      if (attempt == 3) stop("API connection failed after 3 attempts: ", conditionMessage(response))
      Sys.sleep(2^(attempt - 1))
      next
    }
    status <- response$status_code
    header_lines <- strsplit(rawToChar(response$headers), "\r\n", fixed = TRUE)[[1]]
    header_lines <- header_lines[grepl(":", header_lines, fixed = TRUE)]
    headers <- setNames(trimws(sub("^[^:]+:", "", header_lines)),
                         tolower(sub(":.*$", "", header_lines)))
    if (status %in% allowed_status) break
    if (status %in% c(429L, 500L, 502L, 503L, 504L) && attempt < 3) {
      retry_after <- suppressWarnings(as.numeric(tail(headers[names(headers) == "retry-after"], 1)))
      if (!length(retry_after) || is.na(retry_after)) retry_after <- 2^(attempt - 1)
      Sys.sleep(max(1.05, min(retry_after, 15)))
    } else stop(service, " API returned HTTP ", status, ". Query URL: ", request_url)
  }
  result <- list(body = rawToChar(response$content), status = status,
                 provenance = list(service = service, request_url = request_url,
                                   retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                                   http_status = status, cache_file = cache_file, cache_hit = FALSE,
                                   headers = headers))
  saveRDS(result, cache_file)
  result
}

api_provenance <- function(x) {
  value <- attr(x, "api_provenance")
  if (is.null(value) && is.list(x) && "provenance" %in% names(x)) value <- x[["provenance"]]
  value
}

.api_tsv <- function(response) {
  if (!nzchar(trimws(response$body))) return(tibble::tibble())
  if (grepl("^[<{\\[]", trimws(response$body))) {
    stop("Expected a tabular API response but received markup/JSON; inspect the cached response.")
  }
  readr::read_tsv(I(response$body), col_types = readr::cols(.default = readr::col_character()),
                 show_col_types = FALSE, progress = FALSE, name_repair = "check_unique")
}

.uniprot_summary <- function(entry, requested_accession) {
  scalar <- function(value) if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1]])
  primary <- scalar(entry$primaryAccession)
  protein_name <- entry$proteinDescription$recommendedName$fullName$value
  if (is.null(protein_name)) protein_name <- entry$proteinDescription$submissionNames[[1]]$fullName$value
  gene <- if (length(entry$genes)) entry$genes[[1]]$geneName$value else NULL
  tibble::tibble(requested_accession, primary_accession = primary,
                 accession_status = if (is.na(primary)) "unresolved" else if (primary == requested_accession) {
                   "exact_primary_accession"
                 } else "secondary_or_redirected_accession",
                 entry_name = scalar(entry$uniProtkbId), protein_name = scalar(protein_name),
                 gene_name = scalar(gene), organism = scalar(entry$organism$scientificName),
                 organism_taxon = scalar(entry$organism$taxonId),
                 sequence_length = suppressWarnings(as.integer(scalar(entry$sequence$length))))
}

# Returns raw parsed JSON, a compact summary and retrieval provenance.
uniprot_entry <- function(accession, refresh = FALSE) {
  accession <- .api_identifiers(accession, maximum = 1L)
  if (!grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$", accession)) {
    stop("Use a UniProt accession, for example Q2FZ64; a JE2 locus tag is not a UniProt accession.")
  }
  response <- .api_request(paste0("https://rest.uniprot.org/uniprotkb/", accession, ".json"),
                           service = "UniProt", refresh = refresh, allowed_status = c(200L, 404L))
  entry <- if (response$status == 200L) jsonlite::fromJSON(response$body, simplifyVector = FALSE) else list()
  list(summary = .uniprot_summary(entry, accession), entry = entry, provenance = response$provenance)
}

.string_request <- function(method, identifiers, species, extra = list(), refresh = FALSE,
                            allow_unresolved = FALSE) {
  ids <- unique(.api_identifiers(identifiers))
  species <- .api_species(species)
  params <- c(list(identifiers = paste(ids, collapse = "\r"), species = species,
                   caller_identity = "JE2-Project-I-course-analysis"), extra)
  response <- .api_request(paste0(STRING_API_BASE, "/tsv/", method), params,
                           service = "STRING", refresh = refresh,
                           allowed_status = if (allow_unresolved) c(200L, 404L) else 200L)
  result <- if (response$status == 404L) tibble::tibble() else .api_tsv(response)
  attr(result, "api_provenance") <- response$provenance
  result
}

.string_resolve_rows <- function(identifiers, species, raw) {
  queries <- unique(identifiers)
  if (nrow(raw)) {
    required <- c("queryIndex", "stringId", "ncbiTaxonId")
    if (!all(required %in% names(raw))) stop("Unexpected STRING mapping columns; inspect the cached response.")
    indices <- suppressWarnings(as.integer(raw$queryIndex)) + 1L
    if (anyNA(indices) || any(indices < 1L | indices > length(queries))) stop("Invalid STRING query indices.")
    raw$identifier <- queries[indices]
  } else {
    raw <- tibble::tibble(identifier = character(), stringId = character(), ncbiTaxonId = character())
  }
  result <- tibble::tibble(input_row = seq_along(identifiers), identifier = identifiers) |>
    dplyr::left_join(raw, by = "identifier", relationship = "many-to-many") |>
    dplyr::group_by(input_row) |>
    dplyr::mutate(n_candidates_returned = sum(!is.na(stringId))) |>
    dplyr::ungroup()
  same_species <- if (grepl("^STRG", species)) {
    startsWith(result$stringId, paste0(species, "."))
  } else result$ncbiTaxonId == species
  result$mapping_status <- dplyr::case_when(
    is.na(result$stringId) ~ "unresolved",
    result$n_candidates_returned > 1L ~ "ambiguous_candidates",
    is.na(same_species) | !same_species ~ "species_mismatch_review_required",
    result$identifier == result$stringId ~ "exact_string_id",
    TRUE ~ "service_best_match_review_required"
  )
  result
}

# STRING selects a best match; this function does not call that exact validation.
# Every input survives, including unresolved and duplicate input identifiers.
string_resolve <- function(identifiers, species, refresh = FALSE) {
  identifiers <- .api_identifiers(identifiers)
  species <- .api_species(species)
  raw <- .string_request("get_string_ids", identifiers, species, list(echo_query = 1),
                         refresh = refresh, allow_unresolved = TRUE)
  result <- .string_resolve_rows(identifiers, species, raw)
  attr(result, "api_provenance") <- api_provenance(raw)
  result
}

.string_validate_ids <- function(string_ids, species) {
  ids <- unique(.api_identifiers(string_ids))
  species <- .api_species(species)
  if (!all(startsWith(ids, paste0(species, ".")))) {
    stop("Use resolved STRING IDs from the explicitly selected organism. Do not mix organism ID namespaces.")
  }
  ids
}

.string_score <- function(required_score) {
  if (length(required_score) != 1L || is.na(required_score) ||
      required_score < 0 || required_score > 1000) stop("required_score must be between 0 and 1000.")
  as.integer(required_score)
}

# Induced network: return only edges whose two endpoints were explicitly queried.
# This removes automatic single-query neighbours if the server adds any.
string_network <- function(string_ids, species, network_type = "functional",
                            required_score = 700, refresh = FALSE) {
  ids <- .string_validate_ids(string_ids, species)
  network_type <- match.arg(network_type, c("functional", "physical"))
  result <- .string_request("network", ids, species,
                            list(network_type = network_type, required_score = .string_score(required_score), add_nodes = 0),
                            refresh = refresh)
  original_n <- nrow(result)
  if (original_n) {
    if (!all(c("stringId_A", "stringId_B") %in% names(result))) stop("Unexpected STRING edge columns.")
    provenance <- api_provenance(result)
    result <- dplyr::filter(result, stringId_A %in% ids, stringId_B %in% ids)
    if ("score" %in% names(result)) result$score <- as.numeric(result$score)
    attr(result, "api_provenance") <- provenance
  }
  attr(result, "network_scope") <- list(network_type = network_type, explicit_query_ids = ids,
                                        server_edges = original_n, excluded_neighbour_edges = original_n - nrow(result))
  result
}

# Partners can include unmeasured proteins: they are public context, not new data.
string_partners <- function(string_ids, species, network_type = "functional",
                             required_score = 700, limit = 10L, refresh = FALSE) {
  ids <- .string_validate_ids(string_ids, species)
  network_type <- match.arg(network_type, c("functional", "physical"))
  if (length(limit) != 1L || is.na(limit) || limit < 1 || limit > 100) stop("limit must be between 1 and 100.")
  result <- .string_request("interaction_partners", ids, species,
                            list(network_type = network_type, required_score = .string_score(required_score), limit = as.integer(limit)),
                            refresh = refresh)
  if ("score" %in% names(result)) result$score <- as.numeric(result$score)
  attr(result, "network_scope") <- list(network_type = network_type, explicit_query_ids = ids,
                                        partners_are_public_context = TRUE)
  result
}

string_functional_annotation <- function(string_ids, species, refresh = FALSE) {
  ids <- .string_validate_ids(string_ids, species)
  .string_request("functional_annotation", ids, species, refresh = refresh)
}

# KEGG clients below are intended for this academic course project. The service
# permits academic use and requests <=3 calls/sec; our shared limiter uses <=1/sec.
.kegg_identifiers <- function(entry_ids, maximum = 10L) {
  ids <- unique(.api_identifiers(entry_ids, maximum = maximum))
  is_entry <- grepl("^([A-Za-z][A-Za-z0-9_-]*:[A-Za-z0-9_.-]+|[A-Za-z]{1,5}[0-9]{5})$", ids)
  if (!all(is_entry)) {
    stop("Supply explicit KEGG entry IDs, e.g. an organism:gene cross-reference or map00010; not a whole database name.")
  }
  ids
}

.kegg_parse_links <- function(body) {
  if (!nzchar(trimws(body))) return(tibble::tibble(source_id = character(), target_id = character()))
  lines <- strsplit(trimws(body), "\n", fixed = TRUE)[[1]]
  if (any(lengths(strsplit(lines, "\t", fixed = TRUE)) != 2L)) {
    stop("Unexpected KEGG link response: expected exactly two headerless tab-separated columns.")
  }
  readr::read_tsv(I(body), col_names = c("source_id", "target_id"),
                  col_types = readr::cols(.default = readr::col_character()),
                  show_col_types = FALSE, progress = FALSE) |>
    dplyr::distinct()
}

kegg_links <- function(entry_ids, target = "pathway", refresh = FALSE) {
  ids <- .kegg_identifiers(entry_ids)
  target <- match.arg(target, c("pathway", "module", "ko", "brite", "enzyme", "reaction", "compound", "pubmed"))
  response <- .api_request(paste0("https://rest.kegg.jp/link/", target, "/", paste(ids, collapse = "+")),
                           service = "KEGG", refresh = refresh)
  result <- .kegg_parse_links(response$body)
  attr(result, "api_provenance") <- response$provenance
  attr(result, "query_scope") <- list(entry_ids = ids, target = target,
                                      entries_without_returned_links = setdiff(ids, result$source_id))
  result
}

kegg_entry <- function(entry_id, refresh = FALSE) {
  id <- .kegg_identifiers(entry_id, maximum = 1L)
  response <- .api_request(paste0("https://rest.kegg.jp/get/", id), service = "KEGG",
                           refresh = refresh, allowed_status = c(200L, 404L))
  list(entry_id = id, status = if (response$status == 200L) "retrieved" else "not_found",
       text = if (response$status == 200L) response$body else "", provenance = response$provenance)
}
