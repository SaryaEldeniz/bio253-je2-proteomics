# Querying public protein resources from R

These clients send small, explicitly supplied public identifier lists. They do
not upload your proteomics tables, fold changes, sample names or intensities.
Anonymous access is sufficient; no API key is embedded in the project.

## Start with an exact UniProt accession

From the project root:

```r
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
stk <- uniprot_entry("Q2FZ64")
stk$summary
stk$entry$comments       # raw parsed annotation JSON
stk$entry$features       # domains/sites where annotated
stk$entry$sequence       # sequence of this reference entry, not necessarily JE2
api_provenance(stk)
```

The summary preserves the requested and returned primary accessions, organism
taxon, protein/gene name and sequence length. A secondary/redirected accession is
labelled separately from an exact primary-accession match. A 404 is represented
as unresolved. Transport errors and other HTTP failures are reported as errors.
The complete parsed JSON remains available for detailed annotation inspection.

The JE2 mapping table uses reference-strain orthologs. Check the returned organism
before selecting the STRING organism. Keep exact phosphosite coordinates attached
to the original JE2 protein sequence, not an unaligned reference sequence.

## Resolve identifiers in a declared STRING organism

```r
species <- stk$summary$organism_taxon
resolved <- string_resolve(c("Q2FZ64", "Q2FZ65"), species = species)
resolved
```

The query preserves every input, including unmapped and repeated inputs. STRING
selects a best match; the wrapper labels that `service_best_match_review_required`
rather than calling it independently validated. Check the name, organism and
reference mapping before using that match. Multiple returned candidates and
species mismatches have explicit review statuses. An already resolved STRING ID
returned unchanged within the selected organism is labelled `exact_string_id`.

`ncbiTaxonId` must match a numeric organism selection. A custom organism selection
(`STRG...`) is checked against its STRING ID prefix. Do not silently use a broad
species taxon when the response is a different strain/organism namespace; review
and explicitly select the intended STRING reference organism instead.

## Network among selected proteins versus their partners

After reviewing the mapping, explicitly select the accepted STRING IDs:

```r
# Example: after checking that rows 1 and 2 are the intended proteins/organism,
# explicitly choose those rows. Change these row numbers for your own query.
reviewed_rows <- c(1, 2)
accepted <- resolved$stringId[reviewed_rows]
stopifnot(!anyNA(accepted),
          !any(resolved$mapping_status[reviewed_rows] %in%
                 c("unresolved", "ambiguous_candidates", "species_mismatch_review_required")))

edges <- string_network(accepted, species, network_type = "functional", required_score = 700)
partners <- string_partners(accepted, species, network_type = "physical", required_score = 700, limit = 10)
terms <- string_functional_annotation(accepted, species)
api_provenance(edges)
attr(edges, "network_scope")
```

`string_network()` returns only edges whose two endpoints are in the requested
set, filtering out automatic neighbours if STRING adds any. It may return no
edges for a single protein. `string_partners()` deliberately retrieves context
proteins outside your selected set. Those partners are not automatically measured
or differentially regulated in your experiment. Join them back to your master
map/results to determine measurement status.

The default requested score threshold is 700 on STRING's 0-1000 input scale.
Returned `score` values may use a 0-1 scale. A functional association is broader
than a physical-interaction prediction. Neither constitutes a direct
kinase-substrate validation. This wrapper does not yet perform enrichment; later
ORA/GSEA needs a declared measured background and directional comparison plan.

## Caching and errors

- `results/api_cache/` contains response bodies and provenance in local RDS files.
  It is excluded by the project's existing results rule.
- The cache key includes endpoint and query parameters. Cached queries make no
  new network call; `refresh = TRUE` deliberately retrieves a new response.
- `api_provenance(result)` exposes URL, UTC retrieval time, status, response
  headers and cache location. UniProt release headers, if supplied, are retained.
- STRING calls use the explicit `https://version-12-5.string-db.org/api` release
  URL, as illustrated in the official documentation consulted on 11 November 2026.
  The local course STRING exports can be a different version and namespace.
- Requests are spaced by at least 1.05 seconds per service. Each has a 15-second
  timeout and at most three attempts; transient HTTP/network failures use bounded
  backoff. Other HTTP errors stop with the endpoint/status visible.
- Network restrictions, unavailable services or an unsupported custom organism
  remain visible failures; a cache is not fabricated as a substitute live result.

The functions use existing `curl`, `jsonlite` and `digest` packages by namespace.
Their presence was verified locally; install them on a new machine if needed.

## KEGG gene/pathway context for academic use

KEGG's public REST endpoint is for academic use by academic users. This course
project fits that setting. Its documented maximum is three requests per second;
these clients use the same more conservative one-request-per-second limiter.

Use the KEGG cross-reference returned by a reviewed UniProt entry, not a guessed
organism code or a UniProt accession passed as a KEGG gene:

```r
kegg_refs <- Filter(function(ref) identical(ref$database, "KEGG"), stk$entry$uniProtKBCrossReferences)
gene_id <- kegg_refs[[1]]$id   # first inspect whether a KEGG cross-reference exists
pathways <- kegg_links(gene_id, target = "pathway")
pathways                    # source_id and target_id, no input header was assumed
gene_record <- kegg_entry(gene_id)
cat(gene_record$text)
api_provenance(pathways)
```

`kegg_links()` accepts up to ten explicit entry IDs and returns distinct
source/target pairs. An entry without returned links is listed in the
`query_scope` attribute; no link does not mean the gene has no biological role.
Supported target collections include pathway, module, KO, BRITE and enzyme.
`kegg_entry()` retrieves one entry's text plus provenance; HTTP 404 is labelled
`not_found`. Pathway membership describes an annotation, not a demonstrated
change in pathway activity. No enrichment test is performed by these clients.

## What has a supported client here?

UniProt annotation JSON, STRING mapping/network/partners/category annotation,
and KEGG entry/link lookup have clients. BioCyc, AureoWiki, NCBI and PubMed remain useful web resources,
but this step does not pretend that their web links are working REST clients.
Additional clients can be added for a specific question after verifying the
corresponding official API, its identifiers and any access terms.

## Official documentation

- [STRING API](https://string-db.org/help/api/): mapping, species selection,
  pinned versions, request spacing, networks, partners and functional annotation.
- [UniProt programmatic access](https://www.uniprot.org/help/api).
- [UniProt query formats](https://www.uniprot.org/help/api_queries).
- [UniProt website API paper](https://doi.org/10.1093/nar/gkaf394): exact-entry
  endpoints and supported JSON/tabular formats.
- [KEGG API manual](https://www.genome.jp/kegg/rest/keggapi.html): link/get URLs,
  entry identifiers and tabular versus flat-text output.
- [KEGG API access](https://www.kegg.jp/kegg/rest/): academic-use conditions and
  request-rate limit.

Run the network-free client checks with `Rscript tests/check_web_resources.R`.
