# Standalone volcano legend and titles

Use these assets with the three total-proteome interaction panels and three phosphosite interaction panels. They use the existing Arial poster style and matching blue/orange/grey symbols. No plots or data were changed.

- [Shared boxed legend](volcano_shared_legend_box.svg): 24-pt text, enlarged symbols and the common selection rule. Approximately 278 × 28 mm; white interior and thin grey-blue border.
- [Total-proteome heading](title_total_proteome_interactions.svg): bold Arial 32 pt, approximately 260 × 14 mm; transparent background.
- [Phosphosite heading](title_phosphosite_interactions.svg): bold Arial 32 pt, approximately 252 × 14 mm; transparent background.
- [Shared interaction formula using G1/G2](interaction_formula_G.svg): bold 24-pt mathematical text, approximately 205 × 32 mm; transparent background. Separate from either heading, with no explanatory last line.
- [Square legend with selection criteria](volcano_shared_legend_square.svg): 100 × 100 mm, 24-pt text.
- [Square legend without selection criteria](volcano_shared_legend_square_no_selection.svg): 80 × 80 mm, 24-pt text. Retains the five colour/symbol labels, including the grey category.
- [All six SVGs](volcano_legend_and_titles_SVGs.zip).

Exact physical dimensions are in [asset_sizes.tsv](asset_sizes.tsv) and in the SVG roots. PNG previews are beside each SVG. Keep the default physical size when placing them to retain these font sizes. All SVG content is vector paths; no fonts need to be installed to preserve the lettering.

The sign colours refer to the interaction effect, `(A_TSB - A_SN) - (B_TSB - B_SN)`, not simply an increase/decrease under PASN. Grey identifies estimates failing one or both selection criteria. Circles and crosses identify the source-model categories. Their global model label does not guarantee complete replicate detection in every contrast group.

The shared formula uses **G1** and **G2**, the first and second genotype named in the corresponding panel (for example, stk and ancestor in `stk vs ancestor`). The ratios use unlogged abundance; their log2 difference equals the difference of the normalized log2 abundances used for the source contrast. Positive interaction values correspond to a more positive SN-versus-TSB response in G2 than G1. Mathematical subscripts use the standard plotmath size reduction. The earlier combined heading/formula/key files were replaced with heading-only versions, so no redundant combined variants remain in this asset folder or ZIP.

The separate panels already omit repeated legends and main titles, making them suitable for these assets. Retain the Pyk S538 omission note attached to the cropped left phosphosite panel; the common legend does not replace this panel-specific disclosure.

Recreate from `Sarya_George.Rproj`:

```r
source("steps/10_integration_poster/scripts/poster_a0/export_volcano_labels.R",
       encoding = "UTF-8")
```
