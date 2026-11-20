# Step 9 — Proteins with exclusive significant phosphosites

The residue-specific Venn diagram shows 14 hits across TSB and SN. These hits represent 13 distinct JE2 phosphosites on 12 proteins. CsbD-family protein site S5 is counted in both media, and the same protein has a T52 hit in SN. [Venn diagram](../figures/figure09_phospho_residue_venn.png) · [Counts by residue](step09_phospho_residue_venn_counts.tsv)

In each medium, a site qualified when its observed STYH estimate met FDR < 0.10 and |log₂ fold change| > 1 in exactly one of the two mutant-versus-ancestor contrasts. The protein names and locus tags follow the JE2 search database. UniProt IDs are the corresponding BLAST-based matches to the NCTC 8325 reference strain. [Selection code](../scripts/09_render_residue_venn.R) · [Annotation parser](../../../../course_starter/R/parse.R) · [NCTC 8325 proteome](https://www.uniprot.org/proteomes/UP000008816)

| UniProt ID (NCTC 8325) | Protein name in JE2 | JE2 locus tag | Exclusive hits by medium and mutant |
| --- | --- | --- | --- |
| Q2FV52 | Transglycosylase IsaA | `B7H15_14290` | TSB `stk` S194 |
| Q2FZY9 | CsbD family protein | `B7H15_04690` | TSB `stp` S5, SN `stp` S5 and T52 |
| Q2G227 | Phosphopentomutase | `B7H15_00775` | TSB `stk` T85 |
| Q2FZK0 | Cytochrome ubiquinol oxidase subunit I | `B7H15_05485` | TSB `stp` T552 |
| Q2FXM9 | Pyruvate kinase | `B7H15_09250` | TSB `stp` T537 |
| Q2FXL6 | Universal stress protein UspA | `B7H15_09315` | TSB `stp` T157 |
| Q2FWE5 | Serine hydroxymethyltransferase | `B7H15_11735` | TSB `stk` Y51 |
| Q2FZH6 | Phosphocarrier protein HPr | `B7H15_05615` | SN `stk` S52 |
| Q2FX13 | Hypothetical protein | `B7H15_10475` | SN `stk` S51 |
| Q2FZD2 | Thiol reductase thioredoxin | `B7H15_05950` | SN `stp` S71 |
| Q2FYI5 | Cell division protein GpsB | `B7H15_07625` | SN `stk` T93 |
| Q2FZ28 | HslU–HslV peptidase ATPase subunit | `B7H15_06530` | SN `stp` T217 |
