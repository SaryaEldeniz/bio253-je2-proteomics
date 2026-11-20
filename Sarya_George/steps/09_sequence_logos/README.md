# Step 9 · Residues, sequence logos and protein structures

This step groups the phosphoresidue counts, exclusive observed-site Venn diagrams, sequence logos and protein/structure figures. The exclusive set contains 14 medium-specific hits representing 13 distinct JE2 sites on 12 proteins; repeated sites across media are documented in the report.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | R residue-count, residue-Venn and position-0 logo entry scripts. |
| [scripts/phosphosite_figures/](scripts/phosphosite_figures/) | Retained R extraction plus Python sequence-logo, linear-protein and structure-rendering bundle, with its dependency list. |
| [reports/](reports/) | Exclusive-site protein list and residue/Venn count tables. |
| [reports/structures_and_logos/](reports/structures_and_logos/README.md) | Sequence logos using ±10-AA JE2 windows, protein models, 3D structures and rebuilding instructions. |
| [figures/](figures/) | Original residue summary, residue-specific Venn and position-0 logo PNGs. |
| [notebook/](notebook/) | Step 9 selection and interpretation notes. |

Required local inputs: `derived/dea_phospho_STYH.tsv.gz`, `derived/protein_annotation.tsv.gz` and `raw/resources/JE2_proteins.fasta`. The structural/feature bundle also uses cached or public UniProt, Pfam and AlphaFold records; follow its [existing requirements and build guide](reports/structures_and_logos/README.md).

```r
source("steps/09_sequence_logos/scripts/09_render_residue_summary.R")
source("steps/09_sequence_logos/scripts/09_render_residue_venn.R")
source("steps/09_sequence_logos/scripts/09_render_central_residue_logo.R")
```

The position-0 logo describes which residue is phosphorylated; the ±10-AA logos show surrounding sequence context. JE2 site coordinates and mapped reference-protein coordinates are distinguished in the structure/sequence documentation.

[Next: integrated report and poster material](../10_integration_poster/README.md) · [Project overview](../../README.md)
