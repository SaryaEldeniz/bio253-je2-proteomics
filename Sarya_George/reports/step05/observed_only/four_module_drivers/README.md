# Proteins driving the four selected KEGG GSEA patterns

These are all GSEA leading-edge proteins for the observed-only analysis. They are core statistical contributors, not proven causal regulators, direct Stk1/Stp1 substrates, or an exact additive decomposition of NES.

The significant-module tables contain all leading edges for the five module/interaction pairs passing study FDR<0.1. The all-12 table additionally includes the seven nonsignificant module/interaction pairs, explicitly flagged. Individual proteins need not pass the protein-level hit rule to contribute to GSEA.

Four modules: purine metabolism, oxidative phosphorylation, nitrogen metabolism and aminoacyl-tRNA biosynthesis. No LOD-imputed source models are included.

Interaction effect = (A_TSB-B_TSB)-(A_SN-B_SN). With R_g=SN_g-TSB_g, the effects are R_anc-R_stk, R_anc-R_stp and R_stk-R_stp. R values below are supplied model estimates, not flux measurements.

[Detailed table for significant modules](leading_edge_significant_modules.csv) / [Compact protein list](protein_list_compact.csv) / [All twelve module-comparison pairs](leading_edge_all_12_comparisons.csv) / [Module overview](module_overview.tsv)

| Module | Interaction | NES | Module FDR | Leading-edge proteins | Individually selected proteins |
| --- | --- | ---: | ---: | ---: | ---: |
| Purine metabolism | stk vs anc: PASN response difference | 1.960 | 0.00142 | 21 | 8 |
| Purine metabolism | stp vs stk: PASN response difference | -1.829 | 0.00614 | 23 | 5 |
| Oxidative phosphorylation | stk vs anc: PASN response difference | -1.891 | 0.0152 | 12 | 0 |
| Nitrogen metabolism | stp vs anc: PASN response difference | 1.865 | 0.0198 | 5 | 2 |
| Aminoacyl-tRNA biosynthesis | stp vs stk: PASN response difference | 1.870 | 0.0054 | 12 | 0 |

All listed identifiers and product descriptions come from the supplied mapping/source tables. Reference ortholog annotations do not guarantee an identical JE2 sequence or correct assignment of every broad pathway term. Flagged memberships (for example IsaA/SceD in purine metabolism or RbfA/peptidase M16 in oxidative phosphorylation) require review; they remain in the list to reproduce the reported NES.

Checks: each leading edge matches the cached fgsea output, every contributor belongs to its tested set and observed-only universe, all IDs are uniquely mapped, and each interaction equals its difference-of-differences reconstructed from the within-genotype responses.
