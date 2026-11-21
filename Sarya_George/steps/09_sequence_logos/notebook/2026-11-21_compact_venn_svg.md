# 2026-11-21 · Compact Venn SVG for the poster

## Request and changes

Sarya requested an SVG of the Step 9 exclusive source-observed Venn diagram, followed by a compact poster layout and publication to GitHub. Codex produced and visually inspected the compact vector figure. The updated repository already uses workflow-step folders, so this export belongs in `steps/09_sequence_logos/figures/`.

The five comparisons and counts are unchanged. Only spacing, label length and the common caption were adjusted. The blank intersections continue to omit shared sites; they do not indicate zero overlap. Circle areas are schematic.

| Medium / residue | Only stk | Only stp |
| --- | ---: | ---: |
| TSB / pS | 1 | 1 |
| TSB / pT | 1 | 3 |
| TSB / pY | 1 | 0 |
| SN / pS | 2 | 2 |
| SN / pT | 1 | 2 |

## Inputs and reproduction

- Input: the existing [aggregate count table](../reports/step09_phospho_residue_venn_counts.tsv), identical to the local source table used for the original export.
- Output: [compact SVG](../figures/figure09_phospho_residue_venn_compact.svg), 6.5 × 3.45 inches, freely scalable.
- Code: [base R export script](../scripts/09_render_residue_venn_compact.R). Run from the group project directory or repository root as documented in the step README.

No raw measurements, FASTA sequences or biological reanalysis are needed for this graphics-only export. No selection thresholds or counts were changed. Base R graphics were used because the installed R library lacked the project packages; the script can still use the standard project setup when available. The SVG was checked to contain vector elements rather than an embedded raster image and rendered for visual inspection.

GitHub publication was explicitly authorised by Sarya. Unrelated local edits and raw data are excluded from this commit.
