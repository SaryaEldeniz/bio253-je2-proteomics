# Compact IsaA, Pyk and QoxB structure figure

## Question

Prepare the approved three-protein SVG for GitHub, using the structure views
from the existing Step 9 3D report and the final annotation layout.

## Work performed

The SVG contains the matching source-report views of IsaA (Q2FV52), Pyk
(Q2FXM9) and QoxB (Q2FZK0). Each editable arrow label combines the protein
name, JE2 residue, medium, mutant and UniProt accession. The pLDDT legend is
stacked in the free lower-left area of the IsaA panel.

## Checks and result

- IsaA S194 is labelled TSB · stk.
- Pyk T537 is labelled TSB · stp.
- QoxB T552 is labelled TSB · stp.
- All model images are embedded in the SVG.
- The three arrow labels and the legend are vector elements.
- The source PDF images were matched by their decoded pixels to the cached
  report renders, and the exported SVG was rendered and visually checked.

The [compact SVG](../reports/structures_and_logos/JE2_IsaA_Pyk_QoxB_3D_Compact.svg)
uses NCTC 8325 AlphaFold reference models with aligned JE2 residue numbers.
The molecular images are raster layers inside the editable SVG.

## AI assistance

Codex assembled the selected views, applied the requested layout revisions,
checked the labels and prepared the publication commit.
