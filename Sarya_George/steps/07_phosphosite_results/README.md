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

[Next: mutant Venn comparisons](../08_mutant_venn/README.md) · [Project overview](../../README.md)
