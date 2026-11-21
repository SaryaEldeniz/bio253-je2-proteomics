# 2026-11-21: matched compact volcanoes with LOD estimates

Request: add LOD-imputed estimates to the four phosphosite mutant-versus-ancestor volcanoes underlying the Venn diagrams, and make these and Figure 3's total-proteome interaction volcanoes compact and stylistically consistent for the poster.

Added a rendering entry point in `scripts/poster_a0/compact_volcanoes.R`. It reads existing annotated caches and retains effects, FDRs, selection rules, point colours and plotted point sizes. It produces a combined four-panel phosphosite figure, four individual phosphosite panels, a combined three-panel total-protein interaction figure, three individual total panels and one shared legend, in SVG/PNG/PDF. No models were refitted, no FDRs were recalculated and no new imputation was performed.

The phosphosite population remains JE2 CDS. Each of the four comparisons contains 4,100 sites: 157 observed plus 3,943 source LOD-imputed estimates. Selected counts are 60, 72, 42 and 60 for stk/TSB, stp/TSB, stk/PASN and stp/PASN, respectively. These comprise observed counts 5/6/7/8 and LOD counts 55/66/35/52. The older source totals 60/76/45/63 included additional non-JE2 LOD entries, which remain excluded from this JE2-specific presentation. The total-protein population is unchanged at 2,229 eligible entries per interaction: 1,990 observed and 239 LOD-imputed.

The observed-only Venn selection is unchanged: 14 medium-specific comparison hits correspond to 13 distinct phosphosites on 12 proteins. LOD crosses are added to the volcano view for context, not silently added to those Venn sets. Source-model flags are global; an observed label does not establish 3/3 measured replicates in every comparison group. Site signals remain unadjusted for parent-protein abundance.

The shared visual style uses bold titles, short axis labels, compact spacing, enlarged legend keys and Arial text at least 24 pt at the declared physical sizes. Combined phosphosite and total figures are 529.33 × 190 mm and 529.33 × 115 mm. Individual panels are 171.11 × 115 mm, without repeated legends, to use alongside the shared legend. Their height requires enlarging the earlier poster placeholders instead of shrinking the text below the lecture's minimum.

All four phosphosite panels share x = -8 to +8 and y = 0 to 7. The total panels inherit the combined Figure 3 axes and ticks, including individual exports. Label positions are recalculated using actual Arial text measurements, with boxes kept away from points and previous labels; arrow tips retain source coordinates. The earlier canonical Figure 3 outputs are archived before replacement, and the earlier observed-only phosphosite figures remain intact.

Independent checks of the first completed export verified 16,400 phosphosite x/y coordinates and observed/LOD shapes, 6,687 total-protein coordinates, and 24 named-site label placements with unchanged source effects/FDRs and zero points under the label boxes. Ten compact SVGs parsed with no embedded raster images. Negative tick labels retain ASCII minus signs; PDF extraction and a magnified PNG crop confirmed correct rendering. The preservation/count/placement reports and current output links are in [the compact figure guide](../reports/compact_volcanoes/README.md).

Reproduce from the RStudio project:

```r
source("steps/10_integration_poster/scripts/poster_a0/compact_volcanoes.R",
       encoding = "UTF-8")
```

AI assistance: Codex prepared the R plotting/export code, audited population differences and source-coordinate preservation, checked layout and label collisions, and wrote this record. The work changes figure presentation rather than the underlying analyses.
