# Step 5 · GO/KEGG enrichment and module contributors

ORA tests selected directional hit lists; GSEA tests ranked protein estimates. Use the separately labelled observed-only results for the requested figures that exclude source LOD-imputed models. Earlier broad-universe and Top-N outputs are retained.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | Enrichment runs, public-resource checks, report/figure rendering and four-module driver extraction. |
| [docs/](docs/) | [Methods and analysis choices](docs/step05_methods.md). |
| [reports/](reports/README.md) | Original enrichment tables, annotation coverage, run settings and interpretation. |
| [reports/observed_only/](reports/observed_only/README.md) | Every significant observed-only GO/KEGG result, per-contrast tables and figure index. |
| [reports/observed_only/four_module_drivers/](reports/observed_only/four_module_drivers/README.md) | Module overview, leading-edge proteins and readable protein lists. |
| [figures/](figures/) | Original overview figures; `observed_only/` contains the all-significant GO BP/MF/CC and KEGG panels. |
| [notebook/](notebook/) | Step 5 decisions, display revisions and module-driver notes. |

Required local inputs: Step 4 comparison tables in `results/step04/`, the Step 3 map/aliases and supplied annotation terms from the Week 2/STRING resources. Enrichment caches/full tables remain in `results/step05/`. Gene-set membership and public API checks are documented separately in the methods/resource reports.

```r
source("steps/05_enrichment/scripts/05_run_observed_enrichment.R")
source("steps/05_enrichment/scripts/05_render_all_observed.R")
source("steps/05_enrichment/scripts/05_extract_module_drivers.R")
```

Observed-only GSEA panels are ordered from highest signed NES to lowest; ORA panels are ordered by fold enrichment. Enrichment/module membership supplies context for [Step 6](../06_string_network/README.md); protein selection uses the contrasts and measurement support rather than NES alone.

[Next: STRING candidate/network analysis](../06_string_network/README.md) · [Project overview](../../README.md)
