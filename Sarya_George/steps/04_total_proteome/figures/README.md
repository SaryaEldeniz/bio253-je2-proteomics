# Step 4 poster figures

Use the SVG files for the poster: they contain vector paths and symbols, with no embedded raster screenshots, so they remain sharp when enlarged. Matching PNG previews are saved beside them. [Download all five SVGs](step04_poster_SVGs.zip).

| Figure | SVG | PNG preview |
| --- | --- | --- |
| PASN response within each genotype | [SVG](figure02_total_PASN_volcanoes.svg) | [PNG](figure02_total_PASN_volcanoes.png) |
| Three genotype-by-environment interactions | [SVG](figure03_total_interaction_volcanoes.svg) | [PNG](figure03_total_interaction_volcanoes.png) |
| Selected entry counts by contrast | [SVG](figure04_total_hit_counts.svg) | [PNG](figure04_total_hit_counts.png) |
| Candidate protein measurements | [SVG](figure05_candidate_measurements.svg) | [PNG](figure05_candidate_measurements.png) |
| Example interaction protein measurements | [SVG](figure06_interaction_examples.svg) | [PNG](figure06_interaction_examples.png) |

Figure 3 now uses the compact matched volcano style (21 November): one horizontal legend, bold titles and larger legend symbols, with the original plotted values, colours, point sizes and numeric scales. [Separate panels](compact_poster/) use the same axes as the combined figure. [Compact phosphosite figures and reproduction guide](../../10_integration_poster/reports/compact_volcanoes/README.md) include the matched phosphosite comparisons and shared legend.

The poster style uses bold titles and facet labels, a white background and enlarged colour-legend markers. The data points inside the panels retain their original sizes. Plot values, colour mappings, axis labels, limits, breaks, transformations and facet order are preserved, including the separate y-scales in sample panels. No statistical models were rerun and no significance criteria were changed.

Figures 5 and 6 also have light major grid lines and short horizontal bars marking each genotype's mean in each medium. These bars show the arithmetic mean of the finite, observed normalized log2 values, without replacing missing values. Groups with no observed measurements have no mean bar; the existing `n/3` labels retain the number of measurements available. The means and counts are recorded in [reports/poster_sample_means.tsv](../reports/poster_sample_means.tsv).

Open `Sarya_George.Rproj` and run:

```r
source("steps/04_total_proteome/scripts/04_export_poster_figures.R")
# Run last to retain the current compact Figure 3 and matched phosphosite plots:
source("steps/10_integration_poster/scripts/poster_a0/compact_volcanoes.R",
       encoding = "UTF-8")
```

To update only Figures 5 and 6:

```r
options(step04.poster.figures = c("figure05_candidate_measurements",
                                "figure06_interaction_examples"))
source("steps/04_total_proteome/scripts/04_export_poster_figures.R")
options(step04.poster.figures = NULL)
```

This uses the saved Step 4 tables in `results/step04/total/`, the configured total-proteome abundance file, the identifier map and `config/interesting_proteins.tsv`. The [export script](../scripts/04_export_poster_figures.R) reconstructs the existing plots with the unchanged Step 4 helpers and applies only presentation changes. It checks the built plot data and axes before writing each file.

The [invariance report](../reports/poster_style_invariance_checks.tsv) records five checks for each of the five figures. The layer check compares all original layers; the additional mean-bar layer in Figures 5 and 6 is checked separately against the observed means. The [archive manifest](../reports/poster_original_png_archive.tsv) records the paths and MD5 hashes of the original PNGs, which are retained in `90_archive/figures/step04_before_poster_styling/`. Running the original renderer restores the older appearance; rerun the poster exporter afterward if needed.
