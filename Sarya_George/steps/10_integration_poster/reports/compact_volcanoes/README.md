# Compact matched volcano figures for the poster

These figures add source LOD-imputed estimates to the four JE2 phosphosite mutant-versus-ancestor comparisons and apply the same compact presentation to the three total-proteome interaction comparisons. Effects, FDRs and the original point colours and sizes are retained. No models were refitted, no FDRs were recalculated and no new imputation was performed.

[Download the SVG bundle](compact_volcano_SVGs.zip). SVG, PNG and PDF versions are saved together in the linked figure folders.

## Combined figures

| Figure | SVG | PNG | PDF |
| --- | --- | --- | --- |
| Four phosphosite comparisons, observed and LOD-imputed JE2 estimates | [SVG](../../../08_mutant_venn/figures/volcanoes/compact_poster/phospho_mutants_TSB_SN_JE2_all_models_compact.svg) | [PNG](../../../08_mutant_venn/figures/volcanoes/compact_poster/phospho_mutants_TSB_SN_JE2_all_models_compact.png) | [PDF](../../../08_mutant_venn/figures/volcanoes/compact_poster/phospho_mutants_TSB_SN_JE2_all_models_compact.pdf) |
| Three total-proteome interactions, observed and LOD-imputed estimates | [SVG](../../../04_total_proteome/figures/compact_poster/figure03_total_interaction_volcanoes.svg) | [PNG](../../../04_total_proteome/figures/compact_poster/figure03_total_interaction_volcanoes.png) | [PDF](../../../04_total_proteome/figures/compact_poster/figure03_total_interaction_volcanoes.pdf) |

The four phosphosite comparisons use a two-by-two arrangement. The three total-proteome comparisons use a single row. Each combined figure has one compact legend: grey, below selection thresholds; blue, selected negative effect; orange, selected positive effect; circle, observed source model; cross, LOD-imputed source model. Legend markers are enlarged without enlarging the plotted points. The three named phosphosites retain labels and arrows at their unchanged source coordinates.

## Separate panels and shared legend

Individual panels omit repeated titles, subtitles and legends, so they can be assembled with one shared legend on the poster.

| Comparison | SVG |
| --- | --- |
| Phosphosite: stk versus ancestor in TSB | [SVG](../../../08_mutant_venn/figures/volcanoes/compact_poster/JE2_all_models_stk_vs_anc_at_TSB_compact.svg) |
| Phosphosite: stp versus ancestor in TSB | [SVG](../../../08_mutant_venn/figures/volcanoes/compact_poster/JE2_all_models_stp_vs_anc_at_TSB_compact.svg) |
| Phosphosite: stk versus ancestor in PASN | [SVG](../../../08_mutant_venn/figures/volcanoes/compact_poster/JE2_all_models_stk_vs_anc_at_SN_compact.svg) |
| Phosphosite: stp versus ancestor in PASN | [SVG](../../../08_mutant_venn/figures/volcanoes/compact_poster/JE2_all_models_stp_vs_anc_at_SN_compact.svg) |
| Total-protein interaction: stk versus ancestor | [SVG](../../../04_total_proteome/figures/compact_poster/figure03_interaction_stk_vs_anc_at_TSB_vs_SN_compact.svg) |
| Total-protein interaction: stp versus ancestor | [SVG](../../../04_total_proteome/figures/compact_poster/figure03_interaction_stp_vs_anc_at_TSB_vs_SN_compact.svg) |
| Total-protein interaction: stp versus stk | [SVG](../../../04_total_proteome/figures/compact_poster/figure03_interaction_stp_vs_stk_at_TSB_vs_SN_compact.svg) |
| Shared five-item legend | [SVG](../../../08_mutant_venn/figures/volcanoes/compact_poster/shared_volcano_legend.svg) |

## Population, selection and Venn interpretation

The phosphosite figures retain **JE2 CDS entries** from the existing annotation map and include both source-model categories. Each comparison has **4,100 plotted sites: 157 observed and 3,943 LOD-imputed**. They are the same JE2 population used for the previously displayed observed-only panels, with the LOD-imputed rows added.

Selection remains FDR < 0.10 and absolute log2 effect > 1.

| Comparison | Observed selected | LOD-imputed selected | Total selected |
| --- | ---: | ---: | ---: |
| stk versus ancestor, TSB | 5 | 55 | 60 |
| stp versus ancestor, TSB | 6 | 66 | 72 |
| stk versus ancestor, PASN | 7 | 35 | 42 |
| stp versus ancestor, PASN | 8 | 52 | 60 |

The older unrestricted source totals were 60, 76, 45 and 63. Their additional 0, 4, 3 and 3 selected entries, respectively, are non-JE2 LOD-imputed entries. This compact JE2 view excludes those external entries rather than silently adding them to the poster.

**Adding LOD crosses does not change the observed-only Venn sets.** The exclusive set remains 14 medium-specific comparison hits, representing 13 distinct phosphosites on 12 proteins. The coloured crosses in these expanded volcanoes were not added to that Venn selection. Significance in only one mutant comparison does not establish biological absence in the other strain or a significant difference between the mutants.

The total-proteome figures retain the existing eligible-CDS population: 2,229 entries per interaction, comprising 1,990 observed and 239 LOD-imputed estimates. Their source estimates and thresholds are unchanged.

The source-model flag is global. An observed circle does not guarantee three measured replicates in every relevant group. LOD-imputed estimates remain visibly marked and should not be interpreted as precise, fully measured fold changes. Phosphosite effects are not adjusted for parent-protein abundance and do not measure phosphorylation occupancy.

## Reading the effects

For the four mutant comparisons, positive means higher phosphosite signal in the mutant; negative means higher in the ancestor. PASN is the prolonged SN condition, not the 15-minute pulse.

For the total-protein interaction A versus B, the native effect is `(A_TSB - A_SN) - (B_TSB - B_SN)` on log2 signal. Positive means the TSB-minus-SN response is more positive in A than B. This does not alone establish the direction of the response within either genotype.

All four phosphosite panels use common axes, including when exported individually: x from -8 to +8 and y from 0 to 7. The three total-protein panels also share the combined figure's numeric limits and ticks when exported individually. The original total-proteome combined axes are preserved.

## Poster sizes and text

| Output | Intended physical size |
| --- | --- |
| Combined phosphosite figure | 529.33 × 190 mm |
| Combined total-protein figure | 529.33 × 115 mm |
| Each individual panel | 171.11 × 115 mm |
| Shared legend | Cropped to its content |

At those placed sizes, text is Arial, with 32-pt combined titles, 26-pt panel headings and a minimum of 24 pt for axes, legends, subtitles and phosphosite labels. SVG files declare their physical dimensions in millimetres; retain those widths when importing. The individual panels need **115 mm height**, so enlarge the poster foundation's earlier 66-mm placeholder frames rather than shrinking these figures and their labels. The separate panels plus one shared legend give the most flexibility when adapting the poster layout.

The earlier canonical Figure 3 SVG/PNG/PDF was preserved in `90_archive/figures/step04_before_compact_volcanoes/` before updating its regular figure location. The earlier observed-only phosphosite figures remain available separately.

## Reproduce and inspect

From `Sarya_George.Rproj`:

```r
source("steps/10_integration_poster/scripts/poster_a0/compact_volcanoes.R",
       encoding = "UTF-8")
```

The [rendering script](../../scripts/poster_a0/compact_volcanoes.R) reads the existing annotated Step 4 and Step 7 result caches. It does not execute the older full rendering scripts or rewrite the statistical analyses. Individual measurements and full coordinate caches remain local under ignored `results/`.

- [Population and selected-hit counts](phosphosite_population_and_hits.tsv).
- [Source-point and numeric-axis preservation checks](preservation_checks.tsv).
- [Named-site coordinates and label-placement checks](site_label_placement.tsv).

Label bounds are measured from actual 24-pt Arial text before placement. The boxes are placed in empty regions and checked for point overlap. All three named sites retain their original effects/FDRs. SVGs are vector artwork with no embedded raster screenshots. Visual checks cover the combined figures and individual panels; the minus signs on negative axis ticks were independently confirmed in the cached plot labels, PDF text and PNG rendering.

## Added: three all-model phosphosite interactions

These Step 7 interaction figures use the same compact poster style as the figures above. They retain the original all-source population of **4,235 sites per comparison: 170 observed and 4,065 LOD-imputed**. This population includes 135 non-JE2 entries and is deliberately different from the JE2-only 4,100-site population of the four mutant comparisons above. No models, effects, FDRs or selection rules were changed.

| Figure | SVG | PNG | PDF |
| --- | --- | --- | --- |
| Three phosphosite interactions, compact poster view | [SVG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models.svg) | [PNG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models.png) | [PDF](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models.pdf) |
| Three phosphosite interactions, full range | [SVG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models_full_range.svg) | [PNG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models_full_range.png) | [PDF](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models_full_range.pdf) |
| stk versus ancestor, individual compact panel | [SVG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stk_vs_anc_at_TSB_vs_SN_compact.svg) | [PNG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stk_vs_anc_at_TSB_vs_SN_compact.png) | [PDF](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stk_vs_anc_at_TSB_vs_SN_compact.pdf) |
| stp versus ancestor, individual compact panel | [SVG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stp_vs_anc_at_TSB_vs_SN_compact.svg) | [PNG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stp_vs_anc_at_TSB_vs_SN_compact.png) | [PDF](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stp_vs_anc_at_TSB_vs_SN_compact.pdf) |
| stp versus stk, individual compact panel | [SVG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stp_vs_stk_at_TSB_vs_SN_compact.svg) | [PNG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stp_vs_stk_at_TSB_vs_SN_compact.png) | [PDF](../../../07_phosphosite_results/figures/volcanoes/compact_poster/all_models_interaction_stp_vs_stk_at_TSB_vs_SN_compact.pdf) |

The compact combined figure is 529.33 x 115 mm, and each individual panel is 171.11 x 115 mm. Text, point sizes, colours, enlarged legend markers and label handling match the compact figures above. Combined and individual compact panels share x = -12 to +12 and y = 0 to 5. The full-range combined view retains all points with y = 0 to 9. The canonical Step 7 `figures/volcanoes/phospho_interactions_all_models` SVG/PNG/PDF was updated to the compact view after preserving the original files in `90_archive/figures/step07_before_compact_volcanoes/`.

**One explicitly disclosed display omission:** Pyk S538 in the stk-versus-ancestor interaction is the only source point above y = 5. Its effect is +10.6980053 and FDR 9.9559974e-9 (-log10(FDR) = 8.0019152). It is LOD-imputed, with ancestor TSB detected in 0/3 replicates and each of the other relevant groups in 3/3. This is a different phosphosite from the annotated Pyk T537. Its result remains in the source tables, counts and full-range figure. The compact combined subtitle and affected individual panel disclose its omission.

| Interaction | Source sites | Source selected: observed + LOD | Compact displayed sites | Compact displayed selected |
| --- | ---: | ---: | ---: | ---: |
| stk versus ancestor | 4,235 | 1 + 18 = 19 | 4,234 | 18 |
| stp versus ancestor | 4,235 | 1 + 42 = 43 | 4,235 | 43 |
| stp versus stk | 4,235 | 4 + 56 = 60 | 4,235 | 60 |

IsaA S194, Pyk T537 and QoxB T552 remain labelled in all three comparisons. All nine interaction tests for these sites remain nonsignificant. Adding the interaction figures does not alter the observed-only Venn selection. The global source-model labels and the absence of parent-protein adjustment have the same interpretation described above.

- [Exact omitted estimate and reason](phospho_interaction_display_omission.tsv).
- [Source and compact-display population counts](phospho_interaction_population_and_hits.tsv).
- [Preservation checks](preservation_checks.tsv) and [label-placement checks](site_label_placement.tsv).
- [Step 7 notebook record](../../../07_phosphosite_results/notebook/2026-11-21_compact_interaction_volcanoes.md).

Independent checks verified every compact/full-range point coordinate, source-model shape and selection colour. Exactly the disclosed site was absent from the compact view, and the full-range view retained all source rows. All 27 named-site placements across the combined and individual outputs had zero plotted points underneath their label boxes. All five additional SVGs parsed and contained no embedded raster images.

In the preservation table, `original_point_layers_unchanged = FALSE` for the affected compact view and `numeric_axes_unchanged = FALSE` for compact interactions are intentional: they record the single authorised omission and the shared reduced y range. `source_estimates_unchanged = TRUE` confirms that the retained statistical estimates were not modified. The full-range view has no omission and preserves the original numeric limits.

To reproduce only these phosphosite interactions from the RStudio project:

```r
options(poster.volcano.sections = "phospho_interactions")
source("steps/10_integration_poster/scripts/poster_a0/compact_volcanoes.R",
       encoding = "UTF-8")
options(poster.volcano.sections = NULL)
```

The SVG bundle linked at the top includes these five additional files.
