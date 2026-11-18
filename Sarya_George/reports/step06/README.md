# Step 6: choosing STRING seeds from the protein contrasts

Prepared 2026-11-17 from the October total-proteome export. 65 unique contrast-selected proteins plus Stk1 and Stp1 = 67 upload IDs. No network has yet been requested or expanded; the student will upload IDs to STRING 12.5.

## The main change from step 5

The four previously selected KEGG pathways remain the biological scope. Within that scope, protein selection now depends on protein-level interaction results, measurement support and identifier mapping. The scripts never read NES, pathway FDR or GSEA leading-edge results. Previous enrichment helped define the topic; it does not rank or filter these proteins.

In simple language: first ask which individual proteins respond differently in the mutants, then use STRING to learn what connects them. We did not choose proteins because their pathway had a large NES.

## 1. Recover every annotated module member

The inventory includes 151 pathway-protein memberships representing 148 distinct JE2 entries. Of these, 144 entries have interaction data. A protein can belong to more than one pathway. Members absent from the DEA remain visible in the inventory.

Inputs: results/step03/terms_cds.tsv.gz for the supplied KEGG assignments; results/step03/id_map.tsv for JE2 locus, name, UniProt ortholog and mapping quality; results/step04/total/annotated_all.tsv.gz for unchanged October DEA estimates. No source model was refit and no source FDR was recalculated.

Files: [complete membership inventory](01_module_inventory.csv); [all module-protein interaction rows](02_all_module_interactions.csv).

## 2. Compare each protein across all three interactions

Define R_g = log2(SN_g) - log2(TSB_g), the within-genotype response. The source interaction (A_TSB-B_TSB)-(A_SN-B_SN) is R_B-R_A:

| Saved interaction | Formula | Meaning of a positive effect |
|---|---|---|
| stk vs anc at TSB vs SN | R_anc - R_stk | SN response is more positive in anc than stk |
| stp vs anc at TSB vs SN | R_anc - R_stp | SN response is more positive in anc than stp |
| stp vs stk at TSB vs SN | R_stk - R_stp | SN response is more positive in stk than stp |

More positive can mean a stronger increase or a weaker decrease. Interaction effects do not show absolute SN abundance. All three estimates are retained for every selected protein, including nonsignificant, sparse or LOD-flagged estimates with explicit labels.

Important: I_stp/stk = I_stp/anc - I_stk/anc. These three contrasts are related views of the same data; significance in two contrasts is not independent replication.

The environmental source contrasts are TSB-SN. For the response table only, their effects are negated and confidence-interval endpoints reversed. Their original FDRs are preserved.

Files: [wide decision table with effects, FDRs, source intervals, detection counts and response direction](04_candidate_decision_table.csv); [within-genotype SN-minus-TSB responses](03_within_genotype_SN_minus_TSB.csv); [contrast direction key](contrast_direction_key.tsv).

## 3. Make an explicit contrast-based selection

A protein is shortlisted when at least one of the three interactions has: an eligible noncontaminant CDS; finite effect and supplied protein FDR < 0.1; source estimate_type = observed; at least 2 of 3 measured replicates in every relevant group; and an annotation-eligible singleton UniProt mapping. A flagged pathway assignment is held for review, rather than silently treated as a canonical pathway member.

Two tiers keep effect magnitude separate from statistical support:

- Tier A: 12 proteins also have |interaction log2 effect| > 1 in at least one eligible supported contrast.
- Tier B: 53 further proteins have a supported smaller effect; there is no minimum effect cutoff for this tier.

There is no top-N cutoff and no ranking by minimum FDR or NES. Candidate tables are sorted by pathway, tier and protein name. Protein FDR < 0.1 is the existing source per-contrast decision rule, not a newly calculated FDR for the union of three contrasts. The union is an exploratory network seed set, not a new global discovery list.

The >=2/3 rule is new prioritization for step 6. It does not change step-5 observed-only results: source-observed estimates can still have only 1/3 detections in a group. Those entries stay visible in the full decision table, but that sparse contrast cannot qualify a protein for this upload.

64 of the 65 selected proteins have 3/3 measured values in every TSB/SN group. NrdI has 2/3 in anc_SN and remains explicitly labelled; no missing replicate was filled in.

| Pathway | Annotated members | Members with interaction data | Selected | Tier A | Contrast-supported but held |
|---|---:|---:|---:|---:|---:|
| Purine metabolism | 77 | 77 | 33 | 9 | 7 |
| Oxidative phosphorylation | 32 | 29 | 8 | 0 | 1 |
| Nitrogen metabolism | 15 | 14 | 7 | 3 | 0 |
| Aminoacyl-tRNA biosynthesis | 27 | 27 | 17 | 0 | 0 |

8 otherwise contrast-supported module/protein pairs are held for annotation review. These proteins may still be biologically interesting; the uncertainty concerns their assignment to the named pathway. They have a separate upload file and sample panels. These review flags have not changed historical GSEA membership or NES values.

Files: [selected candidates](05_selected_candidates.csv); [membership review](06_annotation_review.csv); [supported proteins across the full proteome, including outside these four pathways](07_global_contrast_supported_context.csv).

## 4. Inspect the measured values and prepare the manual STRING handoff

Generated four interaction-effect heatmaps and 16 pages of individual sample plots covering all selected proteins, both anchors and the contrast-supported membership-review proteins.

Sample panels show measured normalized log2 abundances, not raw intensity or imputed abundance. They contain only TSB and long SN exposure. Dots are individual biological replicates; horizontal bars are observed group means; n/3 labels show detection counts. Missing values stay NA. No lines imply pairing of replicate numbers. Each protein has its own y-axis, so compare conditions within a protein rather than absolute heights across proteins.

The heatmaps show the same three interaction effects planned for your eventual network halos. All use a shared -2 to +2 colour scale, with printed effect values and * for supplied protein FDR < 0.1. Grey NA means the source effect is LOD/other, sparse (<2/3), or nonfinite; zero is white. Original values remain available in the decision table. These are protein effects, not NES.

[Open the full sample-plot gallery](sample_plots.md).

### Examples that explain the contrast direction

| Protein | Model SN-TSB: anc | stk | stp | Important interaction | Interpretation |
|---|---:|---:|---:|---|---|
| PurF | -2.721 | -4.035 | -3.261 | stk/anc +1.314; FDR 1.97e-5 | Decreases in all three strains; the decrease is stronger in stk |
| AtpD | +0.525 | +0.939 | +0.548 | stk/anc -0.413; FDR 1.73e-5 | Increases in all three strains; the increase is stronger in stk |
| ValS | +0.451 | +0.722 | +0.316 | stp/stk +0.406; FDR 4.49e-5 | Increases less in stp than in stk |

These examples were independently checked against the observed sample means. Do not call a pathway activated/inactivated or infer direct phosphorylation from these abundance differences.

Visual spot checks included the purine sample page containing PurF/PurC/NrdI, the nitrogen and respiratory interaction heatmaps, and the anchor sample page. NrdI has visibly wider replicate spread and one missing anc_SN value. Stp1 has a high SN replicate in the stp strain relative to the other two. These observations are retained for interpretation; no sample was deleted and no model was refit. The shortlist is exploratory and does not replace reviewing the complete sample gallery.

### What to upload to STRING 12.5

Recommended first pass: run the four pathway-specific lists separately. Every list includes Stk1 and Stp1 as biological anchors; these two are included regardless of statistical significance, and should remain visible if disconnected.

| Upload file | Contents |
|---|---|
| [purine_plus_anchors.txt](string_uploads/purine_plus_anchors.txt) | 33 selected purine proteins + 2 anchors |
| [oxidative_phosphorylation_plus_anchors.txt](string_uploads/oxidative_phosphorylation_plus_anchors.txt) | 8 selected respiratory proteins + 2 anchors |
| [nitrogen_plus_anchors.txt](string_uploads/nitrogen_plus_anchors.txt) | 7 selected nitrogen proteins + 2 anchors |
| [aminoacyl_tRNA_plus_anchors.txt](string_uploads/aminoacyl_tRNA_plus_anchors.txt) | 17 selected tRNA-related proteins + 2 anchors |
| [all_selected_plus_Stk1_Stp1.txt](string_uploads/all_selected_plus_Stk1_Stp1.txt) | Combined 67 unique IDs |
| [large_effect_subset_plus_anchors.txt](string_uploads/large_effect_subset_plus_anchors.txt) | Optional 12 Tier-A proteins + 2 anchors; omits the smaller-effect modules, so not the sole analysis |
| [membership_review_plus_anchors.txt](string_uploads/membership_review_plus_anchors.txt) | Separate review set + 2 anchors; do not automatically label these as pathway enzymes |

Use the existing reviewed reference organism namespace (taxon 93061, reference S. aureus NCTC8325) when available and verify every returned match. These are reference ortholog accessions, not necessarily JE2 protein sequences. Do not mix those STRING IDs with the course custom-JE2 STRG namespace. Preserve the original JE2 locus in every join.

1. Upload one list through the multiple-protein search in STRING 12.5.
2. Check returned organism, protein description and accession for every match; resolve missing or multiple matches instead of silently dropping them. Record decisions in [STRING_mapping_review_template.csv](string_uploads/STRING_mapping_review_template.csv). Returned STRING IDs are intentionally blank until your export.
3. Choose a functional association network, high confidence 0.7, and zero added interactors for the initial network. Keep isolated nodes. Save the organism, version, network type, threshold, evidence settings and search date.
4. Export the mapping and interaction/edge TSV files, plus the network image if helpful. Save downloads locally under data/raw/step06_string/ so raw downloads stay git-ignored. Tell me where they are saved; the next scripts can join the exported IDs to the prepared node attributes.
5. Expansion and choosing bridge proteins happen after this first network. No expansion was performed in this preparation stage.

STRING functional associations are not necessarily binding or direct Stk1/Stp1 substrate relationships. A KEGG-derived STRING edge is not independent confirmation of the KEGG pathway result. Official references: [STRING score/network definitions](https://string-db.org/help/scores/), [STRING settings](https://en.string-db.org/help/getting_started/).

### Prepared data for the eventual halos

[08_network_node_attributes.csv](08_network_node_attributes.csv) contains all 67 nodes and all three original effects, FDRs, source-model labels and detection counts. halo_effect_* columns are numeric only for source-observed effects with >=2/3 measured in relevant groups. Their order is inner=stk/anc, middle=stp/anc, outer=stp/stk. Missing/unsuitable quantitative estimates stay NA; nonsignificant eligible estimates remain numeric. Node and edge tables can be joined after manual STRING mapping review.

## Interpretation limits that matter

- These are total-protein abundances. This step does not establish phosphorylation changes or direct kinase/phosphatase targets.
- The stp strain also has strongly reduced Stk1, with a possible polar effect noted in the course guidance. Interpret stp-associated responses as that strain's phenotype, not automatically a pure phosphatase-specific effect.
- Source-observed, detection count, effect size and protein FDR describe different aspects of support. None alone proves biological causation.
- A missing STRING edge is not evidence that the proteins cannot interact. Network confidence and protein FDR are different quantities.

## Reproduce from RStudio

Open Sarya_George.Rproj and run:

```r
source("scripts/06_network/06_prepare_candidates.R")
source("scripts/06_network/06_plot_candidates.R")
source("scripts/06_network/06_write_guide.R")
```

All scripts source R/setup.R and reuse the existing local helpers. The editable module-specific annotation-review input is scripts/06_network/module_annotation_review.csv. No web API credentials are required for this preparation. Sample-level values are stored only in ignored results/step06/, and raw inputs remain outside versioned reports.

## Checks performed

Checked unique protein/contrast rows and UniProt upload IDs; preserved source effects/FDRs; reconstructed interaction effects from within-genotype responses; verified the three-contrast algebra; independently recounted detections from sample measurements; checked selected rows are observed and have >=2/3 support; compared the 65 IDs with a separate source-based audit. Original source intervals are retained without assuming an undocumented confidence level.

See [selection rules and checks](selection_rules_and_checks.txt), [input hashes including the annotation-review input](input_provenance.tsv), and [R session information](sessionInfo.txt).
