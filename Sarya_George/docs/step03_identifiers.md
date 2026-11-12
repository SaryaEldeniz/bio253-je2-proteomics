# Step 3: annotation, aliases and public-resource lookup

Step 2 is skipped at the user's request because the data have already been validated.
This module checks identifier joins and annotation transfer; it does not repeat
PCA, sample-quality analysis or mutant abundance validation.

## Run the local map

1. Open the project and keep the existing data-directory setting.
2. Copy `config/annotation_zip.example.txt` to `config/annotation_zip.txt`; put the
   full location of `JE2-resource_and_mapping2uniprot.zip` on its first line.
   Alternatively set the `JE2_ANNOTATION_ZIP` environment variable.
3. Run `source("scripts/02_identifiers/03_build_identifier_map.R")` from the project root.
4. Read `reports/step03_identifiers.md` and the detailed local `results/step03/` tables.

The code reads `JE2_protein_metadata.tsv` and `JE2_protein_terms.tsv` directly from
the supplied zip. It joins on `protein_key`, keeps the exact JE2 locus/accession,
and resolves supplied STRING IDs through exact aliases. All source entries remain
in the master audit. Ordinary-protein analysis uses `entry_class == "cds"`.

The metadata identifies some contaminants not flagged in the DEA. The original
`read_dea(..., drop_contaminants=TRUE)` alone therefore does not define a clean
ordinary-protein enrichment universe. Use the mapped classification as well.

`eligible_uniprot_annotation` requires a CDS and a quality-checked reference
mapping. These map to the NCTC8325 reference proteome, not necessarily an identical
JE2 sequence. Weak mappings remain visible. Many-to-one UniProt targets are counted
and retained; no score is silently selected for future ranked enrichment.

## Interactive protein lookup

```r
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

# JE2 locus, public accession, or curated name/synonym:
lookup_je2(c("Stk1", "PknB", "Stp1", "B7H15_13375", "Q2FVL2"))

# Literal keyword discovery against names, aliases and descriptions:
lookup_je2("cystine", exact = FALSE)

# Browsable links; these are links, not invented REST clients:
protein_resource_links(c("Stk1", "TcyA"))

# Prepare explicit public queries and inspect any exclusions first:
query <- prepare_uniprot_query(c("Stk1", "Stp1", "TcyA"))
query$mapping
query$excluded
query$ids
```

Every input query is retained. Ambiguous aliases produce multiple rows labelled
`ambiguous`; unknown IDs produce `unresolved` rows. `prepare_uniprot_query()` sends
nothing: it returns a list of eligible public identifiers and excluded cases.

For example, `lookup_je2("Q2FVL2")` shows a good TcyA mapping and a second JE2
entry with a weaker mapping. That accession alone is ambiguous in the supplied
reference map. Query `TcyA` or the exact JE2 locus to select the intended source
protein while keeping the ambiguity visible in the audit.

Stk1/PknB/Stp1 names were added as transparent curated aliases because the supplied
gene-symbol fields are empty for these proteins. `config/interesting_proteins.tsv`
contains nine editable starting candidates. These candidates were not chosen by a
new significance analysis.

## Local interaction lookup, with no network access

```r
partners <- local_interactions(c("Stk1", "Stp1"), required_score = 700)
partners$queries
partners$edges

selected_edges <- local_interactions(c("TcyA", "TcyB", "TcyC"),
                                      within_queries = TRUE, required_score = 700)
selected_edges$edges
```

This uses `string_physical_links.txt`, removes duplicated reverse edges and
excludes non-CDS endpoints. The 700 threshold uses the file's 0-1000 score scale.
Neighbouring proteins are context, not automatically changed proteins or direct
substrates. ID alias coverage, annotation-term coverage and network-edge coverage
are reported separately.

## Live APIs

See `docs/web_resources.md` for the reusable UniProt/STRING functions. Run
`source("scripts/03_web_resources/03_query_web_resources.R")` for the small demonstration: it reads
the first five eligible accessions from the candidate list, retrieves UniProt
annotations, takes the organism from those annotations and checks STRING mappings
against UniProt's STRING cross-references before querying networks and partners.
It also retrieves KEGG gene-to-pathway links using KEGG IDs supplied by UniProt,
and the flat-text annotation of one example gene.

The demonstration also saves edges with familiar JE2 labels. The reusable helper
`annotate_reference_edges(edges, crossrefs)` adds these labels without increasing
the edge count. Multiple JE2 matches are listed and counted explicitly; an unknown
public partner stays unmapped rather than being declared unmeasured.

For one protein:

```r
stk <- uniprot_entry("Q2FZ64")
stk$summary
stk$entry$features
stk$entry$comments
```

For several proteins, follow the explicit organism/mapping review in the web
resource guide. The local custom STRING IDs and the online reference-organism IDs
are different namespaces. Do not mix their edges or gene sets without mapping.

## Files produced

- `results/step03/id_map.tsv`: one row per source protein key, with annotations,
  classification, UniProt mapping quality, exact local STRING ID and coverage flags.
- `lookup_aliases.tsv.gz`: searchable identifier/name/description vocabulary.
- `terms_cds.tsv.gz`: ordinary-protein category memberships for later enrichment.
- `excluded_weak_or_many_to_one.tsv`: entries needing exclusion or review.
- `mapping_coverage.tsv`: denominator and coverage summary for total and phospho parents.
- `interesting_protein_lookup.tsv` and `resource_links.tsv`: example candidate lookups.
- `example_local_physical_partners.tsv`: example local physical partners.
- `reports/step03_identifiers.md`: shareable aggregate report.
- `results/step03_api/` and `results/api_cache/`: local API responses and provenance.
- `reports/step03_live_api.md`: created only if the complete live demonstration succeeds.

Detailed annotation maps/results remain local by default; reproducible code,
documentation, aggregate reports and selected figures are the GitHub content.
