# Step 4 · Total-proteome comparisons

This step annotates the supplied October DEA and organizes the planned contrasts. It preserves the source effects/FDRs and records source-model and measured-replicate support. The August screenshot and October export are reconciled in a separate report.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | Comparison assembly, count reconciliation, original figure rendering and `04_export_poster_figures.R` for the poster SVG/PNG style. |
| [docs/](docs/) | [Usage](docs/step04_usage.md) and [interpretation](docs/step04_interpretation.md). |
| [reports/](reports/) | Comparison summary, source-count reconciliation, input hashes and code checks. |
| [figures/](figures/) | Total-proteome PASN and interaction volcanoes, hit counts and per-sample candidate plots, in SVG and PNG. See the [figure guide](figures/README.md). |
| [notebook/](notebook/) | Dated Step 4 notes. |

Required local inputs: `derived/dea_total.tsv.gz`, `derived/abundances_total.tsv.gz`, the Step 3 identifier map and the STYH metadata used by the combined comparison assembly. Full annotated tables stay in `results/step04/`; shared summaries are in this step's `reports/`.

```r
source("steps/04_total_proteome/scripts/04_assemble_comparisons.R")
# Redraw saved comparisons with the original styling:
source("steps/04_total_proteome/scripts/04_render_figures.R")
# Export all five figures in the updated poster style (SVG and PNG):
source("steps/04_total_proteome/scripts/04_export_poster_figures.R")
```

The poster exporter preserves plot values, point sizes, colour mappings, facet order and every panel's axis ranges, transformations, breaks and labels. It applies bold titles, clean panels and larger legend markers. Figures 5 and 6 include readable grid lines and horizontal mean bars for each genotype within each medium. The original PNGs are retained in [the figure archive](../../90_archive/figures/step04_before_poster_styling/). Checks are recorded in [reports/poster_style_invariance_checks.tsv](reports/poster_style_invariance_checks.tsv).

For `TSB_vs_SN`, positive means higher in TSB; negative means higher in SN. Interaction signs compare genotype-specific responses and need the formulas in the interpretation guide. LOD-imputed estimates remain flagged.

[Next: ORA/GSEA](../05_enrichment/README.md) · [Project overview](../../README.md)
