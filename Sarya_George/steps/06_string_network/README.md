# Step 6 · Candidate selection and STRING networks

This step brings the contrast-supported protein candidates, selected UniProt lists, STRING exports/API annotations and network figures together. Stk1/Stp1 are network anchors. The reports distinguish source/course memberships from fresh STRING KEGG assignments.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | Candidate preparation/plots, upload-list export, STRING import/expansion, API retrieval, network drawing and evidence guides. |
| [reports/](reports/README.md) | Candidate decision tables, sample-plot index, selection rules and readable protein lists. |
| [reports/string_uploads/](reports/string_uploads/) | Original module lists, user-selected lists and anchor-inclusive 49/67/75-protein inputs. |
| [reports/string_network/](reports/string_network/README.md) | Imported/exported seed networks, expansion choices, nodes/edges, mapping and protein review. |
| [reports/poster_networks/](reports/poster_networks/README.md) | Fresh STRING API provenance, current pathway memberships, coordinates and poster-layout checks. |
| [reports/phosphosite_overlap/](reports/phosphosite_overlap/README.md) | IsaA/Pyk/QoxB overlap, precise phosphosite/protein contrasts, localization and primary-literature evidence. |
| [figures/](figures/) | Candidate contrast/sample plots and the original networks. `poster_networks/` contains the 49/67/75-node PNG/PDF/SVG figures. |
| [notebook/](notebook/) | Candidate, user-list, STRING and annotation decisions. |

Required local inputs: Step 4 contrasts and sample abundances, Step 3 identifier mapping, Step 5 module-driver tables, local STRING exports and/or cached public STRING responses. The main exported network uses the functional network with minimum score 0.4. [The complete75 guide](reports/poster_networks/complete75_README.md) documents restored candidates and fresh memberships.

Important entry scripts include `06_prepare_candidates.R`, `06_export_user_selection.R`, `06_import_string_exports.R`, `06_fetch_poster_networks.R`, `06_prepare_complete75.R` and `06_plot_poster_networks.R`. Follow the relevant report's rerun order, as the import-based and API-based workflows use different source files.

The baseline-versus-interaction phosphosite figure is now with [Step 7](../07_phosphosite_results/README.md). The cross-layer overlap investigation stays here with the network evidence.

[Next: phosphosite interpretation](../07_phosphosite_results/README.md) · [Project overview](../../README.md)
