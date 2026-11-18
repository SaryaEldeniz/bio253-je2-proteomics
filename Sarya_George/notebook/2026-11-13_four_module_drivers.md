# 2026-11-13: leading-edge proteins in four KEGG modules

Extracted all leading edges for purine metabolism, oxidative phosphorylation,
nitrogen metabolism and aminoacyl-tRNA biosynthesis across the three observed-only
interaction comparisons. The source fgsea result and cached leadingEdge fields
agree. No new differential models, enrichment tests or protein significance
filters were applied to the leading edges.

The five study-significant cells contain 73 module/protein/interaction entries
and 57 unique proteins. The full twelve-cell table also retains nonsignificant
module comparisons (150 entries, 84 unique proteins), explicitly flagged. Each
leading-edge protein belongs to its actual tested pathway and observed-only
universe. Interaction effects reconstruct from R_anc-R_stk, R_anc-R_stp or
R_stk-R_stp, with R=SN-TSB. All leading-edge rank signs match their NES signs.

The detailed tables include UniProt IDs, JE2 loci, supplied product labels,
source interaction effects and FDRs, ranking statistics, detection patterns and
within-genotype response estimates. Individual hit selection is labelled
separately: leading-edge membership does not require proteinFDR<0.1 and|effect|>1.

Several supplied product descriptions are not obviously canonical for their
assigned KEGG map, including IsaA/SceD/HU/RNA-polymerase delta in purine metabolism
and RbfA/peptidaseM16 in oxidative phosphorylation. These entries are retained
to reproduce the reported NES and flagged for annotation review, without claiming
the broad membership is definitively wrong or silently reassigning proteins.
The supplied asnC label is preserved alongside its asparagine-tRNA-ligase product.

R source: scripts/05_enrichment/05_extract_module_drivers.R. Outputs:
reports/step05/observed_only/four_module_drivers/. Codex extracted/formatted the
tables; a separate read-only agent independently verified all12leading edges,
counts, mapping uniqueness and rank signs. Leading edges identify statistical
contributors, not causal regulators or an additive decomposition of NES.
