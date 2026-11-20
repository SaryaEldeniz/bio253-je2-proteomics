# Phosphosite volcano reproduction

The October STYH DEA is redrawn in the compact classic-theme blue/orange style of the earlier total-proteome volcanoes. Each point is one phosphosite estimate, not one protein. Effects and FDRs are copied unchanged; no new tests or imputation were performed.

## Comparisons and selection

The three Step 7 panels test genotype-by-medium interactions. The four Step 8/9 panels compare each mutant with the ancestor separately in TSB and SN. SN is long PASN exposure, not the 15-minute pulse.

Selection is FDR < 0.10 and absolute log2 effect > 1. Source-contaminant and nonfinite entries are excluded. No new UniProt/CDS eligibility filter is imposed; it would change the population used in the original Step 7-9 figures.

All-model panels reproduce the existing Step 7/8 hit counts. Filled circles are source-observed estimates; crosses are LOD-imputed estimates. Observed-only panels exclude all LOD-flagged rows. The observed label is global to the source model, not a promise of three detections in every contrast group.

The separately named JE2_observed_only panels additionally retain only entries classified as JE2 CDS in the existing identifier map. This removes external entries that escaped the source contaminant flag. For example, P02662_S61 is bovine alpha-S1-casein: it is an observed selected stp-vs-stk interaction in the source reproduction, but is excluded from the JE2-specific plot and quality subset. Consequently that panel has three JE2 observed hits rather than four source-observed hits. The original all-model and observed-only reproductions remain available for comparison.

| Contrast | All source hits | Observed source hits | LOD-imputed hits | JE2 observed hits | Quality subset hits |
| --- | ---: | ---: | ---: | ---: | ---: |
| `interaction_stk_vs_anc_at_TSB_vs_SN` | 19 | 1 | 18 | 1 | 1 |
| `interaction_stp_vs_anc_at_TSB_vs_SN` | 43 | 1 | 42 | 1 | 1 |
| `interaction_stp_vs_stk_at_TSB_vs_SN` | 60 | 4 | 56 | 3 | 2 |
| `stk_vs_anc_at_TSB` | 60 | 5 | 55 | 5 | 5 |
| `stp_vs_anc_at_TSB` | 76 | 6 | 70 | 6 | 6 |
| `stk_vs_anc_at_SN` | 45 | 7 | 38 | 7 | 7 |
| `stp_vs_anc_at_SN` | 63 | 8 | 55 | 8 | 8 |

## Localization, replication and parent protein

The annotated selected-site table supplies the best localization probability over all samples, the normalized-abundance detection count in every relevant group, and the matched total-protein effect/FDR. The illustrative quality subset requires JE2 CDS, source-observed, best localization probability >= 0.75 and >=2/3 measured values in each relevant group (four groups for an interaction, two for a simple contrast). It is an additional annotation, not a silent filter on the reproduction plots. Best localization probability does not guarantee high localization confidence in every sample.

Glutamate dehydrogenase T269 (B7H15_04940, Q2FZU5) is a source-observed selected stp-vs-stk interaction, but only 1/3 kinase-mutant SN values were measured. Its best localization probability is 0.8957. It remains in the JE2 observed plot and is flagged outside the conservative quality subset. A large interaction estimate alone is not sufficient evidence of a well-replicated phosphorylation response.

Site and parent-protein estimates were fitted separately. A nonsignificant parent effect does not establish unchanged protein abundance. These plots do not measure phosphorylation occupancy or establish direct kinase substrates. A protein-adjusted significance test has not been performed.

Per-sample values for every selected site are saved locally in results/step07_volcano_reproduction/selected_site_sample_values.tsv.gz for inspection. No raw measurement tables are added to reports/.

## Reading the sign

In mutant-versus-ancestor panels, positive means higher site abundance in the mutant; negative means higher in the ancestor. For an interaction A vs B at TSB vs SN, the effect is (A_TSB - A_SN) - (B_TSB - B_SN). Defining R = SN - TSB, a positive effect means R_B > R_A; a negative effect means R_A > R_B. This does not by itself establish an increase or decrease within either genotype.

All panels within a figure share symmetric x limits and a common y scale. Observed-only and all-model figures each use their own data-derived limits. There is no binning, jitter, coordinate clipping or FDR recomputation. Any source FDR of zero receives a finite coordinate only for display.

## Outputs and rerun

Both the three main interaction panels and the four supplementary mutant-versus-ancestor panels label IsaA S194, Pyk T537 and QoxB T552 with small arrows terminating at the original coordinates. All nine labelled interaction tests are nonsignificant; three of the twelve labelled simple comparisons meet selection. See labelled_interaction_sites.tsv and labelled_mutant_sites.tsv. Labels are placed in nearby empty regions; conservative padded rectangles are checked against every data point, the other labels and the other arrow paths. These placement checks are saved in site_label_placement_checks.tsv. Only the coloured legend markers are enlarged; plotted points retain their original size. Both figure sets use a single explanatory subtitle and omit the technical source-filter footer. The source/JE2 selection policies remain documented above and in the file names.

The figures are in ../../figures/step07/matched_style/ as PNG, PDF and SVG. The individual panels use JE2 source-observed estimates and are available for poster assembly. See volcano_counts.tsv, selected_sites_annotated.tsv, input_checksums.tsv and session_info.txt beside this file.

From Sarya_George.Rproj run:

```r
source("steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R")
```

From the repository root prefix the path with Sarya_George/. The script reuses R/setup.R and the existing site/DEA readers; course data stay in the configured local data directory.
