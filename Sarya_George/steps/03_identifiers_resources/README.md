# Step 3 · Identifiers, annotations and web resources

This step connects JE2 protein keys/locus tags with UniProt ortholog accessions, aliases, STRING identifiers and functional terms. JE2 keys retain the experimental identities; mapped UniProt IDs support external-resource lookup without silently merging distinct JE2 entries.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | `03_build_identifier_map.R` builds the map; `03_query_web_resources.R` demonstrates public-resource queries. |
| [docs/](docs/) | [Identifier guide](docs/step03_identifiers.md) and [API/web-resource guide](docs/web_resources.md). |
| [reports/](reports/) | Matching coverage, API examples, code checks and run/session notes. |
| [notebook/](notebook/) | Dated Step 3 progress. |

Required local inputs: `derived/protein_annotation.tsv.gz`, the DEA inventories used to identify measured proteins, and the configured Week 2 JE2 metadata/mapping archive. Full generated maps and aliases remain in the ignored `results/step03/`; public API caches remain in `results/api_cache/`.

From the group project:

```r
source("steps/03_identifiers_resources/scripts/03_build_identifier_map.R")
source("R/setup.R")
lookup_je2(c("Stk1", "Stp1", "TcyA"))
local_interactions(c("Stk1", "Stp1"))$edges
```

Live queries require network access. The reusable clients are in [R/web_resources/](../../R/web_resources/); local mapping helpers are in [R/identifiers/](../../R/identifiers/).

[Next: total-proteome comparisons](../04_total_proteome/README.md) · [Project overview](../../README.md)
