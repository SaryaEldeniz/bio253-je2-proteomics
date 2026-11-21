# Step 7 · Phosphosite results and interaction figures

This step compares phosphosite effects with parent-protein effects and annotates localization confidence, measured-replicate support and source-model flags. Site abundance and protein abundance are separate measurements; the figures do not establish protein-adjusted phosphorylation occupancy.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | Original Step 7 renderer, labelled-volcano renderer and baseline-versus-interaction figure script. |
| [reports/](reports/) | Original site/parent comparison and counts. |
| [reports/volcano_reproduction/](reports/volcano_reproduction/README.md) | Annotated selected sites, all seven contrast counts, labelled-site estimates and label-placement checks. |
| [figures/](figures/) | Original Step 7 volcano retained. |
| [figures/volcanoes/](figures/volcanoes/) | Latest labelled three-interaction PNG/PDF/SVG figures and individual panels. |
| [figures/baseline_vs_interaction/](figures/baseline_vs_interaction/) | IsaA S194, Pyk T537 and QoxB T552: sample values, baseline contrasts and separate interaction estimates. |
| [notebook/](notebook/) | Step 7 and subsequent phosphosite-figure/formatting notes. |

Required local data: `dea_phospho_STYH.tsv.gz`, `sites_STYH.tsv.gz`, `abundances_phospho_STYH.tsv.gz` and matched `dea_total.tsv.gz` in the configured `derived/` directory, plus the Step 3 identifier map. Complete plotted coordinates/sample tables remain in `results/step07_volcano_reproduction/`.

```r
source("steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R")
```

That shared renderer also writes the four supplementary mutant-versus-ancestor volcanoes to [Step 8](../08_mutant_venn/figures/volcanoes/). `all_models`, `observed_only` and `JE2_observed_only` filenames distinguish the plotted populations. Labels/arrows identify the same three sites in each panel; all nine of their interaction tests are nonsignificant.

The updated [compact all-model interaction figure](figures/volcanoes/compact_poster/phospho_interactions_all_models.svg) uses the matched poster style, one enlarged-symbol legend, and labelled IsaA S194, Pyk T537 and QoxB T552 sites. For the poster, its shared y-axis ends at 5: exactly one high LOD-imputed estimate, **Pyk S538** in the stk-versus-ancestor interaction, is omitted from this display and explicitly noted on the figure. It remains in the [full-range compact companion](figures/volcanoes/compact_poster/phospho_interactions_all_models_full_range.svg), source results and statistical counts. The earlier figure is archived in `90_archive/figures/step07_before_compact_volcanoes/`. This is a display choice, not a quality exclusion or a change to FDR.

Combined and individual SVG/PNG/PDF files are in [figures/volcanoes/compact_poster/](figures/volcanoes/compact_poster/). The regular `phospho_interactions_all_models` files also contain the updated compact poster view. See the [matched volcano guide](../10_integration_poster/reports/compact_volcanoes/README.md) for the omission record, counts, shared legend and all exports. To recreate only these interaction figures, run this **after** the older renderer above:

```r
options(poster.volcano.sections = "phospho_interactions")
source("steps/10_integration_poster/scripts/poster_a0/compact_volcanoes.R",
       encoding = "UTF-8")
options(poster.volcano.sections = NULL)
```

[Next: mutant Venn comparisons](../08_mutant_venn/README.md) · [Project overview](../../README.md)
