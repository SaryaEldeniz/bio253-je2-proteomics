# Step 6, parts 5-7: STRING expansion, contrast halos and protein review

The main network follows the student's explicit choice: STRING 12.5 functional associations, confidence >=0.4. All 47 selected proteins and both biological anchors are retained. Network expansion was performed with public reference identifiers; no experimental abundances or effect tables were uploaded to external services.

## What the seven supplied files contributed

| File | Use in this analysis |
|---|---|
| Short interaction TSV | Main graph edges, combined confidence and evidence channels |
| Protein annotations TSV | Exact UniProt alias mapping, STRING names/descriptions; swapped URL/description headings corrected |
| Functional annotations TSV | Node-to-term membership context; not an enrichment test and contains no FDR |
| Protein FASTA | Checked reference-ID coverage; archived without treating reference sequences as JE2 phosphosite sequences |
| Coordinates TSV | Checked exported node coverage; original coordinates preserved locally; fresh deterministic layouts used for halos |
| Degrees TSV | Independently reconciled degrees against the imported graph |
| PSI-MI XML | Archived supporting edge-evidence export; schema version 2.5 is not STRING version 12.5 |

Raw downloads were copied to the Git-ignored project data/raw/step06_string/seed_export/. Input hashes are in export_input_hashes.csv. Reports and scripts contain derived outputs rather than raw exports.

## Import checks

All 49 requested UniProt IDs match exactly one exported STRING ID by exact alias tokens, all in reference taxon 93061. No gene-name guesses were needed. The imported graph has 204 unique undirected edges and 10 isolated proteins. At >=0.7, it has 104 edges and 14 isolates; this stricter view is supplementary, not the student's main threshold.

The annotation export has its domain_summary_url and annotation content reversed in all 49 rows. Import detects the URL pattern and corrects the fields. Both STRING names and JE2 display names are retained, including asnC/asnS label differences.

The student's screenshots confirm Full STRING network, Evidence edge display, all shown evidence sources enabled, medium confidence 0.400, no first/second-shell added interactors, force-directed layout, and Leiden clustering resolution 0.7. Use clusters in network layout was unchecked; between-cluster edges were set to faded lines. The Leiden resolution is a community-detection parameter, not the edge-confidence threshold. Our >=0.7 sensitivity view is a separate filter. Node-centre colours in the R figures represent supplied pathways, not the exported Leiden communities; those community assignments were not reconstructed.

The version-specific live API exactly reproduces all 204 seed edges and scores. See [mapping](seed_identifier_mapping.csv), [import checks](import_checks.txt), and [export-vs-API differences](export_vs_live_seed_edge_differences.csv), whose differences table has zero rows.

## 5. One controlled round of expansion

Requested up to five highest-scoring partners per original seed, at confidence >=0.4, using the STRING 12.5 interaction_partners API. This returns 242 rows (2-5 per seed), including original proteins; after removing seeds and duplicates, 126 new candidates remain. These are up to five overall partners per seed, not five guaranteed new proteins per seed.

Fetched the induced network for all 175 discovered nodes so the top-five discovery cap does not masquerade as a complete degree count. This discovered-set network has 1,289 edges. It is still not the whole organism's interactome.

Mapped new STRING IDs through exact UniProt STRING cross-references, then through the existing JE2 map. Of the 126 candidates, 115 have exact reference cross-references and unique JE2 mappings, six have exact reference mappings without a JE2 match, four have ambiguous reference/JE2 mappings, and one has no returned exact cross-reference. Unresolved mapping means measurement status is unknown; it is not proof of an unmeasured protein. Original data are assigned only through unique mappings.

A compact overview has a declared display budget of 12 added nodes. First reserve a best-ranked pathway partner for each of the four selected modules, then fill by network context: anchor neighbours, links across modules, links to multiple seeds, quantitative-data availability and confidence. No NES, protein FDR or effect magnitude selects these additions. The 12-node cap is for readability, not a biological threshold or a claim that the other 114 candidates are unimportant.

The resulting compact graph has 61 nodes and 363 edges, retaining every original seed and anchor. Of these edges, 176 score >=0.7 and 187 score between 0.4 and 0.7. All 12 additions have a unique JE2 mapping and quantitative protein data.

| Added protein | UniProt | Reason retained | Links to original seeds |
|---|---|---|---:|
| pheT | Q2FZD8 | Best-ranked partner covering a selected module | 25 |
| rpsG | P48940 | Direct anchor neighbour in one-round functional network | 14 |
| B7H15_06340 | Q2FZ67 | Direct anchor neighbour in one-round functional network | 13 |
| ribBA | Q2FXG1 | Direct anchor neighbour in one-round functional network | 13 |
| ldh2 | Q2G1Y5 | Direct anchor neighbour in one-round functional network | 10 |
| ldh1 | Q2G218 | Direct anchor neighbour in one-round functional network | 10 |
| glmU | Q2G0S3 | Direct anchor neighbour in one-round functional network | 4 |
| gatA | Q2FWY9 | Direct anchor neighbour in one-round functional network | 8 |
| fmt | Q2FZ68 | Direct anchor neighbour in one-round functional network | 8 |
| guaA | Q2G0Y6 | Best-ranked partner covering a selected module | 26 |
| atpC | Q2FWF1 | Best-ranked partner covering a selected module | 6 |
| B7H15_06820 | Q2FYY6 | Best-ranked partner covering a selected module | 4 |

GuaA and PheT link to 26 and 25 original seeds, respectively, within this deliberately selected network. This makes them useful context, not demonstrated global hubs or differentially regulated proteins.

[Every candidate and selection decision](expansion_candidate_decisions.csv) · [All candidate-to-seed links](candidate_to_seed_links.csv) · [Full one-round nodes](nodes_full_one_round.csv) · [Full one-round edges](edges_full_one_round.csv) · [API query provenance](API_provenance.csv).

In simple language: we asked what else connects to your chosen proteins, checked which of those proteins exists in our data, and added a small, explained subset to the picture. Nothing was added because of its NES.

## 6. Three contrast halos

| Ring | Original interaction | Positive means |
|---|---|---|
| Inner | stk vs anc at TSB vs SN | SN-TSB response more positive in anc than stk |
| Middle | stp vs anc at TSB vs SN | SN-TSB response more positive in anc than stp |
| Outer | stp vs stk at TSB vs SN | SN-TSB response more positive in stk than stp |

Blue is negative, near-white is zero, orange-red is positive. All rings use a common -2 to +2 log2-effect scale. Colour values outside this range are clipped only for display; original numbers remain in the tables. A grey halo means no suitable quantitative estimate (LOD/other model, sparse observations or unresolved mapping), not zero. Black ring dots mark supplied protein FDR<0.1, and nonsignificant quantitative effects remain coloured.

Module colours appear in node centres; dark centres identify Stk1/Stp1; dashed centre outlines identify added proteins; daggers identify retained pathway-membership flags. Edges encode STRING confidence independently of protein effect. No ring-width variation exaggerates one comparison. Node positions are drawing aids, not biological distance or direction.

The original, sensitivity and expanded overviews share coordinates. Module and anchor details are separately laid out and explicitly scoped for readability. Isolates remain in the complete views. All plots are exported as PNG/PDF/SVG, with plotting coordinates and halo source values saved.

[Original 49-node halos](../../../figures/step06/string_network/01_selected_49_export_halos.png) · [>=0.7 sensitivity](../../../figures/step06/string_network/02_selected_49_high_confidence_halos.png) · [61-node expanded overview](../../../figures/step06/string_network/03_expanded_network_halos.png) · [Stk1/Stp1 neighbourhood](../../../figures/step06/string_network/04_Stk1_Stp1_neighbourhood_halos.png).

Module zooms: [purine](../../../figures/step06/string_network/05_purine_module_halos.png), [respiration](../../../figures/step06/string_network/06_respiration_module_halos.png), [nitrogen](../../../figures/step06/string_network/07_nitrogen_module_halos.png), [tRNA](../../../figures/step06/string_network/08_tRNA_module_halos.png). These scoped views keep the selected module seeds, both anchors and adjacent retained added partners; other module seeds are omitted for readability. Membership is listed in module_plot_node_membership.csv.

In simple language: the connections come from STRING, while the three rings come from our own experiment. A strong connection and a large protein change answer different questions.

## 7. Compare the additions with the actual experiment

All 61 nodes were checked against the original three protein interactions. The 47 original selected proteins remain: nine have a supported |effect|>1 and 38 have smaller supported effects. Their seven pathway-membership flags remain visible. The two anchors are included for biology, without requiring an interpretable abundance effect.

Among the 12 added proteins, only Ldh1, Ldh2 and B7H15_06340 pass the current protein-FDR rule in a quantitatively supported interaction. The other nine remain useful network context, not newly declared differential-abundance discoveries.

| Added protein | Interaction effect and supplied FDR | What the within-genotype responses show |
|---|---|---|
| Ldh1 / Q2G218 | stk/anc -1.158, FDR 0.000322; stp/stk +0.758, FDR 0.00904 | SN-TSB increases in anc +2.048, stk +3.207, stp +2.449: the increase is strongest in stk |
| Ldh2 / Q2G1Y5 | stp/anc -0.267, FDR 0.0532; stp/stk -0.211, FDR 0.0987 | SN-TSB increases in anc +2.325, stk +2.381, stp +2.592: modestly stronger increase in stp; FDR support is near the threshold |
| B7H15_06340 / Q2FZ67 | stp/stk +0.388, FDR 0.00828 | SN-TSB anc +0.114, stk +0.329, stp -0.058: the response is more positive in stk than stp |

Ldh1 is the clearest new large-effect candidate for following up abundance-response biology. Ldh2 offers a smaller related comparison. B7H15_06340 is annotated as a 16S rRNA methyltransferase and has a smaller supported response difference. Product labels are traceable to the recorded UniProt accessions/course annotation, not newly inferred functions.

Measured replicate plots were generated for all 12 added proteins, using TSB/SN only, with detection counts and no replacement of missing values. Inspect these before interpreting any individual effect. RibBA has only 2/3 observations in stk_TSB and appreciable replicate spread; it is retained as context, not a significant interaction candidate.

[Added-protein measurements, page 1](../../../figures/step06/string_network/added_partner_samples_page01.png) · [Page 2](../../../figures/step06/string_network/added_partner_samples_page02.png).

[Full protein review with original effects/FDRs and within-genotype responses](final_protein_review.csv) · [Readable review](protein_review_readable.md) · [Anchor links and evidence channels](anchor_links_review.csv).

### What the anchor connections do and do not establish

Stk1-Stp1 (0.998) and Stp1-Ndk (0.796) are already present in the original export. Four newly retained anchor links exceed 0.7: Fmt-Stk1 (0.863), Fmt-Stp1 (0.986), B7H15_06340-Stk1 (0.908), and B7H15_06340-Stp1 (0.978). Genomic-neighbourhood/coexpression evidence contributes strongly to these new links; this is not direct substrate validation.

Ldh1-Stp1 and Ldh2-Stp1 both score 0.442 and contain experimental-channel support, but remain moderate-confidence functional associations. GlmU-Stk1 scores 0.513 and is text-mining-only. Do not narrate all displayed anchor edges as physical binding or direct phosphorylation. Check the saved channel scores rather than judging by proximity in the drawing.

In simple language: Ldh1 looks interesting because our measurements show a reproducible mutant-dependent response. Fmt or GuaA can help explain network context even without a significant interaction effect. Being connected to Stk1/Stp1 does not prove that they phosphorylate/dephosphorylate that protein.

## Limits and validation

The three interaction effects are algebraically dependent (I_stp/stk = I_stp/anc - I_stk/anc), so multiple significant comparisons are not independent replication. These are total-protein abundances, not phosphosite changes. The stp strain also has reduced Stk1, with a possible polar effect noted in course guidance. Source FDRs are per-contrast values; no new FDR for this network-selected union was calculated.

Independent checks reproduced all 49 mappings, 204 seed edges, the 12-partner selection, and all 61 proteins' original effects/FDRs across three contrasts. Finite halos require source-observed estimates and >=2/3 relevant detections. Mapping ambiguity remains explicit; degree means degree inside the selected graph. Shared pathway/database evidence is context rather than independent confirmation of KEGG enrichment.

## Reproduce and continue

Open Sarya_George.Rproj and run:

```r
source("scripts/06_network/06_import_string_exports.R")
source("scripts/06_network/06_expand_string_network.R")
source("scripts/06_network/06_review_network_proteins.R")
source("scripts/06_network/06_plot_string_networks.R")
source("scripts/06_network/06_write_network_guide.R")
```

API responses are cached locally under ignored results/api_cache/. Reruns reuse matching requests. Downloads default to the user Downloads folder, configurable with JE2_STRING_DOWNLOAD_DIR. No extra graph packages were installed. Raw files and sequence data stay git-ignored. The CSV node/edge tables are also suitable for Cytoscape import if you want to adjust positions manually.

Next biological choice: use Ldh1 and the rRNA methyltransferase as new abundance-response hypotheses; compare them with your original seed candidates and the phosphosite data before proposing direct Stk1/Stp1 targeting. A protein missing from this compact display remains available in the full one-round tables.

Official method references: [STRING API](https://string-db.org/help/api/), [STRING score and network definitions](https://string-db.org/help/scores/), [UniProt programmatic queries](https://www.uniprot.org/help/api_queries).
