# Complete 75-protein poster network and annotation-review provenance

All eight previously held proteins are now included. The 75 unique proteins match the full sample-plot inventory: 65 previously selected candidates + eight restored candidates + Stk1/Stp1 anchors.

## Where the review label came from

The phrase 'pathway membership under review' was an analyst-authored judgment, not a warning from the lecturer, STRING, KEGG or UniProt. The course pathway membership itself was supplied in Week 2. Comparing those memberships with protein descriptions prompted a manual question; it did not establish an annotation error.

- The original course annotations come from `JE2_protein_terms.tsv` and protein descriptions from `JE2_protein_metadata.tsv` inside the archive configured by `config/annotation_zip.txt`. The reader is `R/identifiers/identifiers.R`; step 3 preserves these terms in `results/step03/terms_cds.tsv.gz`.
- The earliest hard-coded review list is `scripts/05_enrichment/05_extract_module_drivers.R`, lines 47-52.
- The expanded notes are in `scripts/06_network/module_annotation_review.csv`; exact line numbers for these eight are below.
- `scripts/06_network/06_prepare_candidates.R` originally required both quantitative support and no manual review flag. All eight passed quantitative support in at least one interaction; the manual flag alone kept them out of the 67-protein list. That additional exclusion was too restrictive for exploratory networking.

| UniProt | Plot name | Supplied JE2 description | Supplied module questioned | Manual CSV line | Current degree |
| --- | --- | --- | --- | ---: | ---: |
| Q2G010 | B7H15_04455 | thermonuclease | Purine metabolism | 5 | 0 |
| Q2FYK7 | B7H15_07520 | fatty acid-binding protein DegV | Purine metabolism | 11 | 0 |
| Q2FYG2 | B7H15_07765 | DNA-binding protein HU | Purine metabolism | 6 | 2 |
| Q2FUT0 | B7H15_14945 | DNA-directed RNA polymerase subunit delta | Purine metabolism | 7 | 0 |
| Q2FV52 | isaA | transglycosylase IsaA | Purine metabolism | 2 | 0 |
| Q2FXM9 | pyk | pyruvate kinase | Purine metabolism | 4 | 20 |
| Q2FWF8 | sceD | transglycosylase SceD | Purine metabolism | 3 | 0 |
| Q2G0W1 | dabA | hypothetical protein | Oxidative phosphorylation | 12 | 0 |

Examples: [UniProt Q2FXM9](https://www.uniprot.org/uniprotkb/Q2FXM9/entry) identifies Pyk as pyruvate kinase. [KEGG K00873](https://www.kegg.jp/entry/K00873) describes pyruvate kinase and its metabolic pathways. [RCSB annotations for Q2FV52](https://www.rcsb.org/annotations/AF_AFQ2FV52F1) include IsaA in peptidoglycan catabolism. These descriptions explain the initial question about a purine label; they do not exclude broader or indirect relationships.

## Separate, later check against STRING 12.5

The fresh 75-protein functional-annotation export does not assign any of the eight to the four selected KEGG terms (map00230, map00190, map00910, map00970). This is an annotation disagreement with the course table, not proof that the course annotation was wrong. Other KEGG terms may be present: `restored_eight_annotation_provenance.csv` lists every KEGG term returned for each restored protein. 'Context' means outside these four displayed groups; it does not mean no biological function or no pathway membership anywhere.

Evidence: `course_vs_live_KEGG_membership.csv` (set 75), `pathway_membership_live.csv`, and the raw `functional_annotations_75.tsv` under ignored `data/raw/step06_string/poster_api/`. `API_provenance.csv` records versioned request URLs and retrieval times. Exact reference STRING mapping is checked against UniProt cross-references; taxon 93061 is the selected reference, with JE2 source identities retained.

## Figure

The network has 75 nodes and 412 functional associations at score >=0.4; 10 isolated proteins are retained. No automatic interaction partners were added. Edges are functional associations and do not establish direct binding or phosphorylation.

Current four-group counts: purine 24, oxidative phosphorylation 8, nitrogen 7, aminoacyl-tRNA biosynthesis 17; 17 context proteins and two anchors. The eight restored proteins appear in context. Original course memberships and manual notes remain in the data for audit; they are not presented as official warnings in the figure.

Stk1 and Stp1 are pinned centrally for emphasis. Three halos, inside to outside, show original stk-vs-ancestor, stp-vs-ancestor and stp-vs-stk TSB-vs-SN interaction effects. Positive means the TSB-minus-SN response is more positive in the first genotype; negative means more positive in the second. The common colour scale is clipped only for display. Source log2 effects and FDR are unchanged. Grey means no eligible observed estimate, not zero. Halo dots denote source FDR<0.1; eligible halo estimates have observed source labels and >=2/3 detections in every relevant group.

- [PNG](../../../figures/step06/poster_networks/STRING_75_poster_halos.png)
- [Vector PDF](../../../figures/step06/poster_networks/STRING_75_poster_halos.pdf)
- [Vector SVG](../../../figures/step06/poster_networks/STRING_75_poster_halos.svg)

## Reproduce from the Sarya_George RStudio project

```r
source("scripts/06_network/06_prepare_complete75.R")
source("scripts/06_network/06_fetch_poster_networks.R")
source("scripts/06_network/06_compare_poster_annotations.R")
Sys.setenv(JE2_POSTER_SETS = "75")
source("scripts/06_network/06_plot_poster_networks.R")
source("scripts/06_network/06_write_complete75_guide.R")
source("scripts/06_network/06_check_poster_outputs.R")
```

The original 49/67 inputs and figure files remain available. Historical ORA/GSEA tables are unchanged. No phosphosite-overlap investigation was performed in this update.
