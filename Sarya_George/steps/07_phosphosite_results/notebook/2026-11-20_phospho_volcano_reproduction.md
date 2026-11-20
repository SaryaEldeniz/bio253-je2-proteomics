# 2026-11-20 - Recreate phosphosite volcanoes in the total-proteome style

## Question

Recreate the teammate's Step 7-9 phosphosite comparisons using the earlier total-proteome volcano format. Keep genotype-by-medium interactions distinct from mutant-versus-ancestor comparisons within TSB and SN.

## What was run

`steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R` reads the configured October STYH DEA, FragPipe localization table, normalized phosphosite abundances and matched total-protein DEA. It reuses `R/setup.R`, the data readers, contrast plan and replicate-count helper. Source effects and FDRs are unchanged. Thresholds are FDR <0.10 and absolute log2 effect >1.

Figures use the earlier compact classic theme, grey below thresholds, blue negative and orange positive, and dashed selection lines. Full-source plots distinguish observed circles from LOD-imputed crosses. Separate observed-only and JE2 CDS observed-only versions are generated; individual JE2 observed panels are available for poster assembly. All coordinates are exact, without binning.

## Results

| Contrast | Source hits | Observed source hits | LOD-imputed hits | JE2 observed hits | JE2 observed + localization + replication subset |
| --- | ---: | ---: | ---: | ---: | ---: |
| stk/anc interaction | 19 | 1 | 18 | 1 | 1 |
| stp/anc interaction | 43 | 1 | 42 | 1 | 1 |
| stp/stk interaction | 60 | 4 | 56 | 3 | 2 |
| stk/anc in TSB | 60 | 5 | 55 | 5 | 5 |
| stp/anc in TSB | 76 | 6 | 70 | 6 | 6 |
| stk/anc in SN | 45 | 7 | 38 | 7 | 7 |
| stp/anc in SN | 63 | 8 | 55 | 8 | 8 |

Each source-filter reproduction includes 4,235 finite site estimates, including 170 source-observed estimates. The JE2 CDS observed-only figures contain 157 sites per contrast. Exact-reproduction counts agree with the existing Step 7 and Step 8 reports and with an independent R audit.

The extra quality subset requires JE2 CDS, source-observed, best localization probability >=0.75 and >=2/3 measured normalized abundances in every relevant group. It is annotated separately and does not silently filter the reproduction plots.

## Important findings and limits

- P02662_S61 is bovine alpha-S1-casein despite a false source contaminant flag. It is retained only in the clearly labelled source reproductions and excluded from the JE2 CDS versions. This reduces observed stp/stk interaction hits from four to three.
- Glutamate dehydrogenase T269 (B7H15_04940, Q2FZU5) has only 1/3 measured values in kinase-mutant SN. It remains visible in the JE2 observed plot but is outside the extra quality subset. A source-observed label alone does not establish adequate replication.
- All 32 source-observed significant rows across the seven contrasts have best localization probability >=0.75. This is the best probability across samples, not a guarantee of confident localization in every sample.
- In simple contrasts positive means higher in the mutant and negative means higher in the ancestor. For an interaction A/B at TSB/SN, positive means a more positive SN-minus-TSB response in B; negative means a more positive response in A. Neither sign alone establishes the within-genotype response direction.
- Parent-protein estimates are annotated rather than used to infer phosphorylation occupancy. No protein-adjusted statistical model was fitted.

## Verification and AI assistance

Codex wrote and ran the R reproduction script, inspected the rendered combined PNGs, checked source effects/FDRs and join cardinality, reconciled counts against the prior reports, and ran `tests/run_tests.R` successfully (46 assertions). A read-only subagent independently recomputed all seven counts and reviewed contrast directions, population filters and quality annotations. No original Step 7-9 figures or analyses were overwritten. Source measurements and complete point/sample tables remain in ignored data/results locations.

[Figures](../figures/volcanoes) - [Methods, counts and annotated hits](../reports/volcano_reproduction/README.md) - [R script](../scripts/07_recreate_phospho_volcanoes.R)

## Poster annotation update

At the user's request, each of the three interaction panels now has small labelled arrows to IsaA S194, Pyk T537 and QoxB T552. Exact site IDs are specified in the R script, avoiding nearby sites on the same proteins. All nine endpoints were independently checked against the earlier interaction-estimates table; all are source-observed JE2 sites and none passes the selected-hit thresholds. QoxB T552 in stp/ancestor has FDR 0.615981; the other eight FDRs are approximately 1. The arrows therefore point near the bottom axis and do not indicate significant interactions.

The first scope line beneath the interaction title and the source-filter footer have been removed from the poster panels. The replacement subtitle reads: "PASN response = log2(SN/TSB); interaction = response(second genotype) - response(first genotype)." The removed footer explained that the original contaminant flag did not exclude all external proteins, including bovine casein. That qualification remains in the methods report and is handled by the separate JE2-only versions.

Only the coloured legend markers were enlarged (3.8 versus 0.95 for data points). Site coordinates, point sizes, selection thresholds and counts were retained. Interaction PNG/PDF/SVG exports, including the separate observed-only panels, were regenerated. Codex implemented the formatting changes; a read-only subagent verified the nine site identities and the interaction sign convention.

## Clear labels and matching supplementary panels

The user requested labels that do not cover data, plus the same three site labels in the four mutant-versus-ancestor panels. The R script now searches nearby empty regions for the three compact labels in each panel. It checks conservative padded label rectangles against all data points and against previously placed labels, then connects each label to the original site coordinate. Data-point size and position remain unchanged. Placement details are written to `steps/07_phosphosite_results/reports/volcano_reproduction/site_label_placement_checks.tsv`.

All seven comparisons now label all three sites: nine labels in each combined interaction figure and twelve in each combined supplementary mutant figure. All 21 site-contrast combinations are observed JE2 entries. Only the three previously selected TSB comparisons are significant: IsaA S194 in stk/ancestor, and Pyk T537 and QoxB T552 in stp/ancestor. Other labels retain the appropriate nonsignificant point colour. The exact estimates are in `labelled_interaction_sites.tsv` and `labelled_mutant_sites.tsv`.

The supplementary panels now use the same enlarged coloured legend markers, the same short selection/protein-adjustment caption, and a single scientific subtitle: "Within each medium: site abundance difference = log2(mutant/ancestor)." The technical source-filter footer and model-scope subtitle are omitted consistently from both figure sets. Full-source, observed-only and JE2 observed-only selection policies remain documented in the report and file names.

The final layout also prevents leader lines from passing beneath other labels, using vectorised segment/rectangle intersection checks with horizontal/vertical control cases. The completed render produced 39 nonempty PNG/PDF/SVG exports and 84 label placements across the combined and individual figures. Every padded rectangle contains zero plotted points, and all source counts remain unchanged. Codex visually inspected the main all-model and supplementary JE2 observed-only PNGs after the final render.
