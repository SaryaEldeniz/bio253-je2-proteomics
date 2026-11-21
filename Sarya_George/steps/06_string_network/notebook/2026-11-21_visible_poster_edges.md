# 2026-11-21: visible STRING associations for the poster

Request: make the connecting lines of `STRING_75_poster_halos.svg` visible from afar, with no other layout or data changes.

Changed only the 75-node edge style in `scripts/06_plot_poster_networks.R`: colour #697686 to #42566A, width range 0.08–0.36 to 0.25–0.60 mm, and alpha range 0.10–0.38 to 0.55–0.85. Width and opacity still increase monotonically with the same STRING confidence scores over 0.4–1.0. The 49/67 node render styles are unchanged.

Added `scripts/06_render_visible_poster_edges.R` to archive the previous SVG/PNG/PDF once and rerender only the 75-node live-KEGG figure. It checks identical node/edge/membership input hashes, exact equality of all 75 exported node coordinates/metadata (including halo values and labels), exact outline equality, the same 412 edges, and vector-only SVG output. The earlier files remain in `90_archive/figures/step06_before_visible_edges/`; current SVG/PNG/PDF use their usual filenames.

The older renderer's literal plus/minus character was escaped in the R source to avoid this machine's C-locale parsing issue, preserving its displayed text. Codex made the styling/export change, checked numerical/layout invariance and visually inspected the refreshed PNG. No statistical or network-confidence filtering changed, and no API requests were made.

Recreate from `Sarya_George.Rproj`:

```r
source("steps/06_string_network/scripts/06_render_visible_poster_edges.R",
       encoding = "UTF-8")
```
