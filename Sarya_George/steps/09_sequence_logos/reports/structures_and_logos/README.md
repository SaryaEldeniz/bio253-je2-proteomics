# Step 9 phosphosite figures

These figures cover the exclusive observed STYH hits for the two mutant-versus-ancestor
comparisons in TSB and SN. The four groups contain 14 hits across 13 distinct JE2
phosphosites on 12 proteins.

| Figure | Contents |
| --- | --- |
| [Sequence logos (±10 AA)](JE2_Sequence_Logos_Plus_Minus_10_AA.pdf) | A pooled logo, four group logos, and the aligned JE2 sequence windows. Letter height shows the observed amino-acid frequency among available residues. |
| [3D phosphosite structures](JE2_Phosphosites_3D_Structures.pdf) | AlphaFold reference cartoons with highlighted residues, arrows, group labels, and local pLDDT. |
| [Compact IsaA, Pyk and QoxB SVG](JE2_IsaA_Pyk_QoxB_3D_Compact.svg) | The three selected structure views with editable arrow labels, UniProt IDs in the boxes, and a vertical pLDDT legend within the IsaA panel. |
| [Linear protein models](JE2_Phosphosite_Protein_Models.pdf) | JE2 residue coordinates with projected UniProt features, Pfam domains, and labelled phosphosites. |

CsbD S5 belongs to both `stp` groups and contributes once to the pooled logo.
The structural and linear figures map the course's JE2 sequences onto the NCTC 8325
reference sequences. For Q2FX13, JE2 S51 aligns to reference S53.

The compact SVG shows IsaA S194 (TSB · stk), Pyk T537 (TSB · stp) and
QoxB T552 (TSB · stp). Its structure images are embedded raster layers.
Labels, arrows and the legend can be edited in an SVG editor.

## Rebuild the PDFs

Run these commands from the repository root with R and Python 3.12 available.
Install the Python dependencies in a virtual environment. The tested
PyMOL build is pinned in [requirements.txt](../../scripts/phosphosite_figures/requirements.txt).

```text
python -m pip install -r Sarya_George/steps/09_sequence_logos/scripts/phosphosite_figures/requirements.txt
python Sarya_George/steps/09_sequence_logos/scripts/phosphosite_figures/build_figures.py
```

The [build script](../../scripts/phosphosite_figures/build_figures.py)
runs [the R extraction](../../scripts/phosphosite_figures/extract_proteins.R)
from the same STYH table used by the
[residue-specific Venn plot](../../scripts/09_render_residue_venn.R).
It selects observed estimates with FDR < 0.10 and an absolute log₂ fold change > 1,
then identifies hits exclusive to one mutant contrast within each medium.

The local course data directory must contain these files.

- `derived/dea_phospho_STYH.tsv.gz`
- `derived/protein_annotation.tsv.gz`
- `raw/resources/JE2_proteins.fasta`

Set `JE2_DATA_DIR`, use the group's `config/data_dir.txt`, or pass `--data-dir`.
The default is the repository's `data/` directory. If R is installed outside the
executable search path, pass its executable with `--rscript`.

The build retrieves the public UniProt, Pfam, and AlphaFold records when needed.
Pass `--cache-dir` with a local directory outside the repository to reuse those
records and the structure renders. The intermediate table containing course
sequences is created in a temporary directory. PDFs are written to this report
folder by default. Use `--output-dir` to write them elsewhere.

## Sources

- [Step 9 protein and site assignments](https://github.com/SaryaEldeniz/bio253-je2-proteomics/blob/ff9a6546437974b2834cd182ba2c376403019bb4/Sarya_George/reports/step09_exclusive_phosphosite_proteins.md)
- [UniProt sequence annotations](https://www.uniprot.org/help/sequence_annotation)
- [Pfam in InterPro](https://www.ebi.ac.uk/interpro/entry/pfam/)
- [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)
- [WebLogo frequency-based sequence logos](https://weblogo.berkeley.edu/manual.html)
