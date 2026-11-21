# 2026-11-20: Step 4 vector exports for the poster

## Request

Provide SVG versions of all five Step 4 figures and match the recent volcano styling: bold titles, clean white panels and larger legend dots. Preserve the data, colours and axis scales inside the plots.

## Work performed

- Added `steps/04_total_proteome/scripts/04_export_poster_figures.R` as a separate, rerunnable figure entry point. Existing analysis and plotting helper files were left unchanged.
- Reconstructed the same five plot objects from saved October-export results, abundance measurements, the existing identifier map and the existing candidate selection. No models, fold changes, FDRs or selection rules were recalculated.
- Exported figures 02 through 06 as vector SVGs, matching 300-dpi PNG previews and an SVG ZIP bundle. Retained the original figure dimensions.
- Applied bold titles and facet headings, classic white panels and larger colour-legend symbols. Point sizes inside panels, colour mappings and all numeric plotting scales were retained.
- Preserved original PNGs byte-for-byte in `90_archive/figures/step04_before_poster_styling/`; recorded original MD5 hashes in `reports/poster_original_png_archive.tsv`.

## Verification

All five figures passed all five checks in `reports/poster_style_invariance_checks.tsv`: identical built layer data, facet order, axis limits/transformations/breaks/labels, x-axis title and y-axis title. The built-layer comparison also checks point positions, jitter, sizes, colours, shapes, bar coordinates and annotation values. Each SVG was checked for vector content and the absence of embedded raster images. PNG previews were visually inspected for title, panel and legend layout.

## Reproduction

Open `Sarya_George.Rproj`, then run:

```r
source("steps/04_total_proteome/scripts/04_export_poster_figures.R")
```

AI assistance: Codex helped write the presentation-only exporter, documentation and preservation checks. Source statistical estimates and scientific selections were not modified.

## Follow-up: means and grid for Figures 5 and 6

At the user's request, added short horizontal mean bars and light major grid lines only to the two sample-measurement figures. Each bar is the arithmetic mean of the observed normalized log2 abundances for one protein, genotype and medium. Missing values remain missing and groups with no finite measurements have no bar. The `n/3` labels retain the number of available measurements. No groups are connected across media and no error bars or statistical tests were added.

The exporter now supports updating specified figures without rewriting the other figure files. All original data layers, axis scales, colours and facets must pass the existing preservation checks; the new mean layer is checked separately. Derived means and measurement counts are saved to `reports/poster_sample_means.tsv`. Both SVGs and matching PNGs were updated, along with the SVG bundle. Figures 2 through 4 were retained byte-for-byte.
