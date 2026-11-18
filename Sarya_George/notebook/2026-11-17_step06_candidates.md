# 2026-11-17 — Step 6, contrast-first network seeds

Question: which proteins in the four selected KEGG pathways should be inspected
in STRING, based on their protein contrasts rather than pathway NES?

Ran from Sarya_George.Rproj:

```r
source("scripts/06_network/06_prepare_candidates.R")
source("scripts/06_network/06_plot_candidates.R")
source("scripts/06_network/06_write_guide.R")
```

Used the October total-proteome annotated DEA, full supplied KEGG membership
tables, reviewed JE2-to-reference UniProt map and normalized sample measurements.
No NES, GSEA leading-edge or pathway-significance inputs are read by the new
selection scripts. The four pathway topics themselves were selected previously.
No source models or FDR values were changed. No network calls or expansion were
performed: the student will upload the prepared lists to STRING 12.5.

Results:

- Inventory: 151 pathway/protein memberships, 148 JE2 entries; 144 entries have
  interaction data. Four inventory entries remain visible without DEA data.
- Selection: source observed, supplied protein FDR<0.1 in at least one of the
  three interactions, >=2/3 detections in each relevant group, eligible CDS,
  noncontaminant, singleton annotation-eligible UniProt, no flagged membership.
- 65 unique selected proteins: 33 purine, 8 oxidative phosphorylation, 7 nitrogen,
  17 aminoacyl-tRNA. Twelve are large-effect Tier A (9 purine and 3 nitrogen),
  and 53 are smaller-effect Tier B. No top-N cutoff.
- Eight otherwise contrast-supported module/protein pairs held for membership
  review (7 purine, 1 oxidative phosphorylation), retained in audit and separate
  upload list. These flags are not definitive incorrect-annotation findings.
- 67 combined upload IDs including Stk1 and Stp1 biological anchors.
- All selected proteins have all three effects/FDRs exported. Halo-ready columns
  retain nonsignificant eligible effects, but use NA for sparse or source-LOD
  estimates. Complete original estimates and flags remain available.
- Four effect heatmaps and 16 sample pages, TSB/SN only; original sample rows
  stay in ignored results/step06/. Reports contain derived summaries only.
- 64/65 selected proteins have all six groups measured 3/3; NrdI has 2/3 in
  anc_SN. Sample values are normalized observations, not raw intensity.

Checks: source effects/FDRs preserved, unique joins/IDs, detection recount from
measurements, difference-of-differences reconstructed, contrast algebra verified,
upload manifest matches plain-text IDs. An independent agent recomputed the
65 selected IDs from source DEA and full memberships. A summary-only variable
shadowing bug initially counted held proteins in Tier-A module totals; corrected
before delivery, with a reconciliation assertion added. Actual candidate IDs and
upload files were unaffected.

Visual spot checks: nitrogen/respiratory effect maps, purine page 1 and anchor
sample page. NrdI shows wider replicate spread; Stp1 in stp_SN has one high
replicate. No observations were removed. The user should inspect the full gallery
before biological interpretation. Three numerical examples (PurF, AtpD, ValS)
were independently checked against sample means; model-response differences
were at most 0.0155 log2 across these nine example genotype responses.

Interpretation: an interaction is a response difference, not absolute SN
abundance. The three interactions are linearly dependent and are not independent
replication. Protein FDRs are the source per-contrast FDRs; the union is an
exploratory seed list, not a newly corrected global discovery family. stp also
has reduced Stk1 and cannot automatically isolate phosphatase-specific effects.
Total-protein results do not establish direct phosphorylation targets.

Next: student uploads one module list plus anchors at a time; reviews returned
organism/mapping; exports STRING 12.5 mapping and edge TSVs, with functional
network type, confidence 0.7 and initially zero added interactors. Save raw
downloads under local data/raw/step06_string/ (ignored), then join using the
reviewed mapping. Expansion and halo-network drawing follow the handoff.

AI assistance: Codex wrote R scripts, created tables/figures and the guide, and
used a read-only independent agent for scientific review and selection audit.
Human inspection and biological choice remain part of the project.
