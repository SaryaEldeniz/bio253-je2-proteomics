# A0 poster: proposed figure selection

This is a consultation plan, not a completed poster. Approximately two thirds of the poster are available for results. Keep the scientific sequence visible: measured effects, functional summaries, association context, then phosphosite candidate interpretation.

## Methods strip

Use a compact version of [the parallel analysis workflow](../figures/methods_svg/03_analysis_workflow.svg), beside the experimental design. This belongs in Methods and should not consume a large Results panel. Distinguish total-protein interactions from phosphosite simple comparisons and from phosphosite interaction tests.

## Proposed Results allocation

Percentages refer only to the available Results area and are layout estimates.

| Panel | Contents | Approximate space |
| --- | --- | ---: |
| A: measured interaction effects | Three total-protein and three JE2 observed-only phosphosite interaction volcanoes, in the same contrast order | 25% |
| B: functional summaries | One compact KEGG GSEA summary for the three interactions, emphasizing the four modules; retain study FDR and NES direction | 15% |
| C: association context | The complete 75-node STRING network, with readable ring and association legends | 30% |
| D: phosphosite comparison sets | Compact TSB/SN Venn summary; move the full residue split to the supplement | 10% |
| E: three candidate sites | IsaA S194, Pyk T537, QoxB T552: quantitative evidence plus linear maps; small structure thumbnails if legible | 20% |

The [total-protein interaction volcanoes](../../04_total_proteome/figures/figure03_total_interaction_volcanoes.svg) retain observed circles and LOD-flagged crosses. The [JE2 observed-only phosphosite interaction volcanoes](../../07_phosphosite_results/figures/volcanoes/phospho_interactions_JE2_observed_only.svg) exclude source LOD-labelled estimates. Keep these populations explicit. Their corresponding thresholds and measurement support are documented in the source step reports.

For panel B, use the current [observed-only enrichment outputs](../../05_enrichment/reports/observed_only/README.md), rather than mixing them with older primary-analysis Top-N figures. The observed-only interaction KEGG GSEA results contain 7, 3 and 3 significant rows, while KEGG ORA contains only one significant interaction row (purine metabolism in stk versus ancestor, study FDR about 0.062). A single GSEA panel therefore uses space more effectively than three ORA panels. A four-module by three-contrast overview is a proposed new presentation, not a figure already generated here. If drawn, retain nonsignificant cells as such and provide complete tables behind the QR code; do not imply these four modules are the only tested sets.

Use the [75-node network](../../06_string_network/figures/poster_networks/STRING_75_poster_halos.svg) with its native legend. It contains 412 functional associations at score >=0.4, retains 10 isolated proteins and has no automatic added partners. Anchors are positioned centrally for presentation. STRING associations do not establish direct kinase substrates or physical binding; see [STRING's official network definition](https://en.string-db.org/help/scores/).

## Correct interpretation of the phosphosite branch

The exclusive sets contain **14 medium-specific comparison hits, representing 13 distinct sites on 12 proteins**. CsbD S5 occurs in both media. Exclusivity means meeting FDR <0.10 and |log2 effect| >1 in exactly one mutant-versus-ancestor comparison within the same medium. It does not mean biological absence from the other strain or establish a significant difference between the two mutants. See the [complete site list](../../09_sequence_logos/reports/step09_exclusive_phosphosite_proteins.md).

The three overlap candidates pass the relevant TSB simple comparison: IsaA S194 in stk versus ancestor, Pyk T537 and QoxB T552 in stp versus ancestor. All nine corresponding genotype-by-PASN interaction tests are nonsignificant (minimum FDR about 0.616). Thus the [three-site baseline-versus-interaction evidence figure](../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.svg) belongs near the candidate cards, or should be simplified into them. This measurement evidence has priority over large structure images. Parent-protein-adjusted phosphorylation significance has not been established.

For structure thumbnails, use AlphaFold NCTC 8325 reference-model captions and retain JE2-to-reference mapping and pLDDT. The existing PDF atlas places IsaA on 3D page 1 left, QoxB on page 2 right and Pyk on page 3 left. Linear maps are IsaA on page 1 top, QoxB on page 2 top and Pyk on page 2 middle. Predicted placement is structural context, not evidence that phosphorylation changes activity.

The pooled sequence logo contains 13 distinct windows and the four group logos contain only 3, 4, 3 and 4 windows. A small pooled logo or three aligned candidate sequence windows is reasonable if space remains; the complete group logos and alignments belong in the supplement. These are frequency displays, not a validated kinase-recognition motif.

## QR supplement

The [curated supplementary index](poster_supplement/README.md) links complete enrichment tables and plots, four simple-comparison phosphosite volcanoes, observed sample plots, site/parent context, Venn details, sequence/structure atlases, network attributes and reproducible code. It avoids linking raw course measurement tables.

After publication, the intended QR target is:

<https://github.com/SaryaEldeniz/bio253-je2-proteomics/blob/main/Sarya_George/steps/10_integration_poster/reports/poster_supplement/README.md>

This local plan and supplement index do not establish that the new files have been pushed. Open the published URL on a phone and verify every linked figure before printing a QR code. An immutable commit link can freeze the exact supplement used for the final poster.
