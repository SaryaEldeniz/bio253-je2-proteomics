# Dated public-resource checks

These calls use explicit public protein/pathway IDs only. No sample measurements are submitted.
The frozen course STRING-v12 memberships remain the enrichment input; exact snapshot minor release and GO/KEGG dates are unknown.
Live STRING calls are pinned to version 12.5 and reference organism NCTC8325 (taxon 93061).
Only same-taxon STRING matches also verified by UniProt STRING cross-references enter the small context network.
The live functional network is not a binding assay or a direct kinase-target network. Unmeasured/untested functions are not inferred from absent edges.
The custom JE2 organism IDs in the slides/local snapshot are not sent as reference-organism IDs.
Selection: Stk1/Stp1 plus up to six proteins most frequent in the leading edges of the top two study-significant KEGG terms in each focus contrast; ties by accession.
Inspect TSV lookup status/error fields; public service failures do not change enrichment results.

Rerun with run_project_step("enrichment_resources"). Responses and HTTP provenance are cached under ignored results/api_cache/.
