# Poster supplementary figures, tables and reproducible code

This curated index is prepared locally for the poster's QR code. It becomes an online supplement after these files are published to the repository; this page does not establish that they have already been uploaded. Use the published URL of this page for the QR code and verify it on a phone before printing.

The project studies *S. aureus* JE2 ancestor, kinase mutant (`stk`) and phosphatase mutant (`stp`). The poster focuses on TSB versus long PAO1-supernatant exposure (`SN`); the 15-minute pulse is a separate condition. Source effects and FDRs are retained. Unless stated otherwise, selection is FDR < 0.10 and absolute log2 effect > 1. Missing measurements remain missing.

## S1. Full volcano plots and statistical definitions

- [Latest compact poster volcanoes and separate panel index](../compact_volcanoes/README.md), with SVG/PNG/PDF exports and preservation checks. The compact all-model phosphosite interaction view omits one high LOD-imputed Pyk S538 estimate from its display; the [full-range companion SVG](../../../07_phosphosite_results/figures/volcanoes/compact_poster/phospho_interactions_all_models_full_range.svg) retains it. It remains in all statistical counts and source results.
- [Separate headings, bold G1/G2 interaction formula, and shared rectangular/square legends](../../figures/volcano_labels/README.md).
- [Total proteome: three interaction volcanoes (SVG)](../../../04_total_proteome/figures/figure03_total_interaction_volcanoes.svg) and [within-genotype TSB/SN comparisons (SVG)](../../../04_total_proteome/figures/figure02_total_PASN_volcanoes.svg). These total-proteome figures show both source-observed estimates and LOD-flagged estimates, using different symbols.
- [Phosphoproteome: three JE2 observed-only interaction volcanoes (SVG)](../../../07_phosphosite_results/figures/volcanoes/phospho_interactions_JE2_observed_only.svg).
- [Phosphoproteome: four JE2 observed-only mutant-versus-ancestor comparisons in TSB/SN (SVG)](../../../08_mutant_venn/figures/volcanoes/phospho_mutants_TSB_SN_JE2_observed_only.svg).
- [Phosphosite populations, filtering, thresholds and interpretation](../../../07_phosphosite_results/reports/volcano_reproduction/README.md), [counts](../../../07_phosphosite_results/reports/volcano_reproduction/volcano_counts.tsv) and [annotated selected-site estimates](../../../07_phosphosite_results/reports/volcano_reproduction/selected_sites_annotated.tsv).

The simple comparisons and genotype-by-medium interactions test different questions. A significant site difference in TSB does not establish a significant genotype-dependent response to SN. A source-observed label does not guarantee three measured replicates in every group.

## S2. All significant observed-only GO and KEGG enrichment results

- [Poster functional summary (SVG)](../../../05_enrichment/figures/poster_functional_summary/functional_summary_interactions.svg) and [interpretation, exact tables and R exporter](../../../05_enrichment/reports/poster_functional_summary/README.md). The matrix includes 12 distinct KEGG pathways and 13 significant pathway–interaction pairs. The GO-BP ORA panel shows three overlap-based representatives of eight significant terms for stk–ancestor; the other two interactions have no passing GO-BP ORA terms.
- [Complete ORA/GSEA index: every comparison, figure and ordered table](../../../05_enrichment/reports/observed_only/README.md).
- [All significant GSEA GO/KEGG results](../../../05_enrichment/reports/observed_only/gsea_all_significant_GO_KEGG.tsv), [all significant ORA GO/KEGG results](../../../05_enrichment/reports/observed_only/ora_all_significant_GO_KEGG.tsv) and [figure manifest](../../../05_enrichment/reports/observed_only/figure_manifest.tsv).
- [Observed-only analysis settings](../../../05_enrichment/reports/observed_only/settings.txt), [background audit](../../../05_enrichment/reports/observed_only/universe_audit.tsv) and [detailed methods](../../../05_enrichment/docs/step05_methods.md).

These analyses exclude LOD-labelled proteins before enrichment and retain all terms passing study-wide BH FDR < 0.10. They do not add a minimum-replicate filter. The observed-only analysis is distinct from the earlier primary-analysis Top-N figures. GO BP/MF/CC and KEGG are annotation collections; ORA and GSEA are the statistical methods applied to them. NES direction must be read with the native interaction formula, not as pathway activation or metabolic flux.

## S3. Individual measurements, means and parent-protein context

- [Figure 5: candidate total-protein measurements (SVG)](../../../04_total_proteome/figures/figure05_candidate_measurements.svg) and [Figure 6: interaction examples (SVG)](../../../04_total_proteome/figures/figure06_interaction_examples.svg). Points are individual observed normalized log2 values; short bars are means and `n/3` labels report measurement counts.
- [Figure guide](../../../04_total_proteome/figures/README.md) and [derived mean/count summary](../../../04_total_proteome/reports/poster_sample_means.tsv).
- [Purine-module measurements](../../../05_enrichment/figures/step05_candidate_sample_values.png) and [nitrogen-module measurements](../../../05_enrichment/figures/step05_nitrogen_sample_values.png), showing TSB/SN only. These proteins retain the original primary-based selection; these are measurement figures, not a separate enrichment test.
- [Network candidate sample-plot index](../../../06_string_network/reports/sample_plots.md).
- [IsaA S194, Pyk T537 and QoxB T552: baseline measurements and interaction estimates (SVG)](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.svg), with [figure interpretation](../../../06_string_network/reports/phosphosite_overlap/figure_control_vs_interaction.md) and [site-versus-parent results](../../../06_string_network/reports/phosphosite_overlap/README.md).

Phosphosite signal has not been formally adjusted for parent-protein abundance and does not measure phosphorylation occupancy. None of the three named sites passes the interaction FDR threshold in the three PASN interaction tests.

## S4. Exclusive-site Venn diagrams and the complete hit list

- [Residue-specific observed-site Venn figure](../../../09_sequence_logos/figures/figure09_phospho_residue_venn.png), [counts](../../../09_sequence_logos/reports/step09_phospho_residue_venn_counts.tsv) and [protein/site assignments](../../../09_sequence_logos/reports/step09_exclusive_phosphosite_proteins.md).
- [Venn and supplementary-comparison methods](../../../08_mutant_venn/README.md). The older overall Venn includes LOD-labelled hits and is a different analysis population.

The exclusive set contains **14 medium-specific comparison hits, representing 13 distinct JE2 sites on 12 proteins**. CsbD S5 occurs in both media. “Only stk/stp” means significant in exactly one mutant-versus-ancestor comparison within that medium; it does not mean detected only in that strain. A difference in significance is not itself a significant difference between the mutants.

## S5. Full sequence, linear-protein and structure atlases

- [Sequence logos and aligned JE2 windows, plus/minus 10 residues (PDF)](../../../09_sequence_logos/reports/structures_and_logos/JE2_Sequence_Logos_Plus_Minus_10_AA.pdf).
- [Linear protein models and phosphosite positions (PDF)](../../../09_sequence_logos/reports/structures_and_logos/JE2_Phosphosite_Protein_Models.pdf).
- [Reference 3D structure cartoons and highlighted sites (PDF)](../../../09_sequence_logos/reports/structures_and_logos/JE2_Phosphosites_3D_Structures.pdf).
- [Atlas guide, coordinate mapping, annotation sources and reproduction](../../../09_sequence_logos/reports/structures_and_logos/README.md).
- Poster exports: [sequence-logo SVG](../../../09_sequence_logos/reports/structures_and_logos/JE2_Sequence_Logos_Plus_Minus_10_AA.svg), [linear-map SVG](../../../09_sequence_logos/reports/structures_and_logos/JE2_Phosphosite_Protein_Models.svg), [3D-atlas SVG](../../../09_sequence_logos/reports/structures_and_logos/JE2_Phosphosites_3D_Structures.svg), [individual page SVGs](../../../09_sequence_logos/reports/structures_and_logos/svg_pages/) and [SVG ZIP bundle](../../../09_sequence_logos/reports/structures_and_logos/JE2_structures_and_logos_SVGs.zip). The original raster structure renders remain embedded; logo/linear artwork and page labels are vector.

These atlases cover the full selected set, rather than only the three network-overlap proteins. The pooled logo uses 13 distinct sites; the group logos use small sets and are exploratory sequence context, not established kinase-recognition motifs. Structures are reference-model context and do not establish a phosphosite's functional effect. JE2 and reference residue coordinates must be distinguished.

## S6. The complete 75-protein STRING network and annotation audit

- [75-node poster network (SVG)](../../../06_string_network/figures/poster_networks/STRING_75_poster_halos.svg) and [complete guide](../../../06_string_network/reports/poster_networks/complete75_README.md).
- [Current node attributes](../../../06_string_network/reports/poster_networks/nodes_75_live.csv), [associations](../../../06_string_network/reports/poster_networks/edges_75.csv), [degrees](../../../06_string_network/reports/poster_networks/degrees_75.csv) and [identifier mappings](../../../06_string_network/reports/poster_networks/identifier_mapping_75.csv).
- [Network settings](../../../06_string_network/reports/poster_networks/settings.csv), [API provenance](../../../06_string_network/reports/poster_networks/API_provenance.csv) and [course-versus-current KEGG membership audit](../../../06_string_network/reports/poster_networks/course_vs_live_KEGG_membership.csv).

The network retains 75 nodes, 412 functional associations at confidence score >=0.4 and 10 isolates, with no automatically added partners. Its pathway outlines use fresh STRING KEGG annotations; enrichment uses the frozen course memberships. Central anchor placement and functional edges do not establish direct binding or kinase/phosphatase-substrate relationships. The three halo contrasts are algebraically related, not independent replication.

## S7. Reproduction, resources and provenance

- [Project guide and RStudio entry point](../../../../README.md), [complete repository map](../../../../../REPOSITORY_GUIDE.md) and [R code index](../../../../00_shared/docs/R_code_index.tsv).
- [Identifier mapping and web-resource workflow](../../../03_identifiers_resources/README.md), [total-proteome workflow](../../../04_total_proteome/README.md), [enrichment workflow](../../../05_enrichment/README.md) and [network workflow](../../../06_string_network/README.md).
- [Phosphosite workflow](../../../07_phosphosite_results/README.md), [mutant overlaps](../../../08_mutant_venn/README.md), [sequence/structure workflow](../../../09_sequence_logos/README.md) and [integration/poster materials](../../README.md).

The repository supplies code, derived statistical summaries and rendered figures. Course raw measurements and local data caches remain outside the published supplement. Each workflow documents its required local inputs, settings and rerun commands.
