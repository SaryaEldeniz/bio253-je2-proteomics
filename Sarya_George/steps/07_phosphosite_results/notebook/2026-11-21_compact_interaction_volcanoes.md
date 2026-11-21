# 2026-11-21: compact phosphosite interaction volcanoes

Request: apply the matched compact poster style to the three all-model phosphosite interaction volcanoes, allowing the single high LOD-imputed point in the first panel to be omitted from the compact display.

The existing annotated Step 7 coordinates were reused. No statistical models were refitted, no FDRs were recalculated and no new imputation was performed. The native interaction is `(A_TSB - A_SN) - (B_TSB - B_SN)` on log2 signal; positive means a more positive TSB-minus-SN response in A than B. These interaction tests differ from the within-medium mutant-versus-ancestor tests used for the Venn diagrams.

The requested all-model interaction figure retains its original source population: **4,235 sites per comparison**, comprising 170 observed and 4,065 LOD-imputed source estimates. This includes 135 non-JE2 entries per comparison, so it is deliberately different from the JE2-only population of 4,100 sites in the compact four-panel mutant figure. Source selected counts remain **19, 43 and 60** for stk versus ancestor, stp versus ancestor and stp versus stk, using FDR < 0.10 and absolute log2 effect > 1. Their observed/LOD splits are 1/18, 1/42 and 4/56. Source-model labels are global and do not establish complete replicate detection in every relevant group.

The compact view uses common x limits of -12 to +12 and y limits of 0 to 5, also in the individual exports. Exactly one point is above this display range: **Pyk S538**, `lcl_CP020619.1_prot_ARG46254.1_1820_S538`, in the stk-versus-ancestor interaction. Its effect is +10.6980053, FDR 9.9559974e-9 and -log10(FDR) 8.0019152. It has a LOD-imputed source model and localization probability 0.9227. Relevant measured detection is ancestor TSB 0/3, ancestor SN 3/3, stk TSB 3/3 and stk SN 3/3. This is a different site from the labelled Pyk T537. Omitting it makes no change to its statistical result or selection status.

The combined compact subtitle explicitly discloses this omission; the first individual panel has its own notice. Thus the first compact panel displays 4,234 points and 18 selected points, while retaining the source total of 19 selected results. The other two panels display all 4,235 points and their 43 and 60 selected results. A full-range companion uses y = 0 to 9 and includes every source estimate. The exact omitted result and the source/display counts are recorded in the [omission table](../../10_integration_poster/reports/compact_volcanoes/phospho_interaction_display_omission.tsv) and [population table](../../10_integration_poster/reports/compact_volcanoes/phospho_interaction_population_and_hits.tsv).

The three named sites, IsaA S194, Pyk T537 and QoxB T552, retain their source coordinates and labels with arrows in each panel. All nine corresponding interaction tests remain nonsignificant. Label boxes are repositioned using actual 24-pt Arial measurements, with zero plotted points underneath them. Circles identify observed source models; crosses identify LOD-imputed source models. Selected negative and positive effects remain blue and orange; unselected points remain grey. Plot point sizes are unchanged, and only legend markers are enlarged.

Combined outputs are 529.33 x 115 mm; individual outputs are 171.11 x 115 mm. Bold titles, compact legends, 26-pt panel headings and a minimum 24-pt Arial text match the other compact poster volcanoes. SVG physical dimensions are explicit millimetres. The canonical `figures/volcanoes/phospho_interactions_all_models` SVG/PNG/PDF now contains the compact view. Its original files were preserved in `90_archive/figures/step07_before_compact_volcanoes/`.

Outputs, each in SVG/PNG/PDF:

- [Compact combined figure](../figures/volcanoes/compact_poster/phospho_interactions_all_models.svg).
- [Full-range combined figure](../figures/volcanoes/compact_poster/phospho_interactions_all_models_full_range.svg).
- [stk versus ancestor panel](../figures/volcanoes/compact_poster/all_models_interaction_stk_vs_anc_at_TSB_vs_SN_compact.svg).
- [stp versus ancestor panel](../figures/volcanoes/compact_poster/all_models_interaction_stp_vs_anc_at_TSB_vs_SN_compact.svg).
- [stp versus stk panel](../figures/volcanoes/compact_poster/all_models_interaction_stp_vs_stk_at_TSB_vs_SN_compact.svg).

Independent checks compared all 12,704 compact and 12,705 full-range plotted coordinates with the source effects/FDRs. Observed/LOD shapes and selection colours matched exactly. The only missing compact source key was Pyk S538 in the first interaction; the full-range figure omitted none. All 27 named-site label placements across combined and individual files had zero points under their boxes. Five SVGs parsed successfully and contained no embedded raster images.

The preservation report marks compact interaction numeric axes as changed and the affected point layer as changed. Those flags describe the explicitly authorised display crop and single omission, rather than changed statistical estimates. `source_estimates_unchanged` remains TRUE. The [compact figure guide](../../10_integration_poster/reports/compact_volcanoes/README.md) links all reports and the SVG bundle.

Reproduce only this section from `Sarya_George.Rproj`:

```r
options(poster.volcano.sections = "phospho_interactions")
source("steps/10_integration_poster/scripts/poster_a0/compact_volcanoes.R",
       encoding = "UTF-8")
options(poster.volcano.sections = NULL)
```

AI assistance: Codex adapted the existing R rendering script, audited source estimates and the single display omission, checked label collisions and vector outputs, and wrote this record. This work changes figure presentation, with the disclosed omission retained in the full-range companion.
