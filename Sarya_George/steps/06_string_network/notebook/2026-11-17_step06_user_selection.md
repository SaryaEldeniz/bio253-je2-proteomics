# 2026-11-17 — User-selected STRING lists

The student selected 43 explicit UniProt IDs after inspecting the sample plots.
Preserved that exact selection, with no additional NES, FDR or effect filter.

Ran `source("steps/06_string_network/scripts/06_export_user_selection.R")` from the group
project. The input text file and R script are saved under steps/06_string_network/scripts/.

Results in steps/06_string_network/reports/string_uploads/user_selected/:

- Purine: 20 selected, 22 with Stk1/Stp1 anchors.
- Oxidative phosphorylation: 5 selected, 7 with anchors.
- Nitrogen: 3 selected, 5 with anchors.
- Aminoacyl-tRNA: 15 selected, 17 with anchors.
- Combined: 43 exact selected IDs; 45 with anchors.
- Seven previous purine-membership flags retained, with all seven proteins
  explicitly included in both module and combined uploads as the user requested.
- Matching 45-node halo-attribute table contains all three source effects/FDRs,
  source-model flags and detection counts. Unsuitable quantitative halo estimates
  remain NA. No models or source results changed.

Checks passed: all 43 requested IDs resolve to unique JE2 entries; no duplicate
upload IDs; exported combined selection matches input order exactly; 45-node
table matches combined-with-anchors file; 135 unique protein/interaction rows.
Earlier automatic lists remain intact. R-code index regenerated (73 scripts).

No STRING request/export, network expansion, commit or push performed. Next:
student uploads chosen lists, reviews mappings and provides STRING 12.5 exports.
Codex wrote and ran the export R script and prepared the outputs.

## Additional four proteins

The student added G2G2S0, Q2FYG7, Q2G086 and Q2FZI9. The local mapping uniquely
identifies Q2G2S0 as PurB; interpreted G2G2S0 as this typographical error and
recorded the correction explicitly. Other additions: Ndk, deoxyribonucleotidase
B7H15_04035 and PurF. All four belong to the supplied purine map.

Reran the export script. Current selection: 47 proteins, or 49 with anchors;
purine: 24 proteins, or 26 with anchors. Other module counts remain 5/3/15.
Regenerated manifests, mapping template and the exact 49-node attribute table.
Replaced fixed-size checks with input-derived checks for future additions.
