# User-selected STRING upload lists

These lists preserve all 47 UniProt IDs selected by the student from the sample plots. No additional NES, FDR, effect-size or automatic-shortlist filter was applied. The previous automatic lists remain unchanged.

Files ending _selected.txt contain only the chosen proteins. Files ending _selected_plus_anchors.txt additionally contain Stk1 (Q2FZ64) and Stp1 (Q2FZ65). Plain-text uploads have one accession per line and no header.

| Module | Selected | With Stk1/Stp1 | Membership flags retained |
|---|---:|---:|---:|
| Purine metabolism | 24 | 26 | 7 |
| Oxidative phosphorylation | 5 | 7 | 0 |
| Nitrogen metabolism | 3 | 5 | 0 |
| Aminoacyl-tRNA biosynthesis | 15 | 17 | 0 |

Combined: 47 chosen proteins, or 49 with both anchors. All IDs were matched uniquely to JE2 entries in the existing inventory.

Seven selected proteins retain the previous purine-membership review flag: Q2G010, Q2FYK7, Q2FYG2, Q2FUT0, Q2FV52, Q2FXM9 and Q2FWF8. They are included in the purine and combined lists exactly as requested. This flag concerns pathway assignment, not whether the protein is biologically interesting. See annotation_review_retained.csv for product descriptions and reasons.

selected_manifest.csv records names, JE2 loci, modules and selection origin. network_node_attributes.csv contains this exact 49-node set with the three original effects/FDRs and halo-ready columns. Source-LOD or <2/3-supported estimates have NA halo values; original estimates remain available with flags.

Review STRING's returned organism and mapping. Keep the reference organism namespace consistent; do not mix reference-ortholog STRING IDs with custom JE2 STRG IDs. No STRING requests or network expansion were performed by this script.

Reproduce from Sarya_George.Rproj:

```r
source("scripts/06_network/06_export_user_selection.R")
```

Input: scripts/06_network/user_selected_uniprot_ids.txt. Supporting inputs are the step-6 membership inventory and original October annotated total-proteome contrasts.
