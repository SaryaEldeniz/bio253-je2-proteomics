# Step 4: run and explore the comparisons

Open this project's `Sarya_George.Rproj` in RStudio. Keep the Step 3 local map
and data settings. Run:

```r
source("scripts/04_analysis/04_assemble_comparisons.R")
```

This reads the existing lecturer DEA and normalized measurements, adds the
identifier map, records support counts, selects comparison rows and regenerates
the five figures. It does not refit the statistics or repeat the skipped sample-QC
module. It works entirely from local inputs; no API requests are made.

## Start with one readable table

```r
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
step04_contrast_plan() |> dplyr::select(contrast, question)

hits <- read_comparison("interaction_stk_vs_anc_at_TSB_vs_SN",
                         conservative_only = TRUE)
View(hits |> dplyr::select(display_name, locus_tag, uniprot_acc,
                           diff, FDR, direction, group_detection_pattern))
```

- Default `read_comparison(name)` returns all ordinary-CDS rows, including those
  not selected by the thresholds.
- `hits_only = TRUE` restricts to source FDR < 0.1 and |diff| > 1.
- `conservative_only = TRUE` selects conservative hits: the same thresholds,
  source estimate_type observed, and >=2/3 measurements in every relevant group.
  It implies hit selection even when hits_only is not separately supplied.
- `dataset = "phospho_STYH"` reads preliminary site-level tables. These require
  the later Step 7 localization and parent-protein evaluation before final use.

## Files and fields

`results/step04/total/annotated_all.tsv.gz` retains every source row in the ten
selected contrasts, including excluded classes for audit. `ordinary_CDS_all.tsv.gz`
is its ordinary-protein subset. Corresponding phosphosite files are under
`results/step04/phospho_STYH/` and count sites rather than proteins.

Each layer has a `contrasts/` directory with these suffixes per contrast:

| File suffix | Contents |
|---|---|
| `_all.tsv` | All ordinary-CDS rows for that comparison. |
| `_hits.tsv` | Source-FDR/effect-size selected rows. |
| `_conservative_hits.tsv` | Selected rows with the declared conservative support. |
| `_positive_hits.tsv` | Selected rows with positive native source effects. |
| `_negative_hits.tsv` | Selected rows with negative native source effects. |

Uncompressed TSV files can be opened in RStudio or imported into Excel. The
interactive reader avoids having to decompress the combined tables yourself.

| Field | Meaning |
|---|---|
| `diff`, `p.value`, `FDR`, `statistic` | Unchanged supplied statistical results. |
| `direction` | Sign interpreted in words for this particular contrast. |
| `selected_hit` | Ordinary CDS, finite diff/FDR, FDR < 0.1 and |diff| > 1. |
| `estimate_type` | Original source-model flag, preserved verbatim. |
| `group_detection_pattern` | Measured values/3 in each group actually used by this contrast. |
| `any_relevant_group_unobserved` | Whether any of those groups has no observed values. |
| `quantitative_support` | Conservative source label plus >=2 measurements per relevant group; not significance on its own. |
| `conservative_hit` | selected_hit AND quantitative_support. |
| `observed_group_contrast` | Descriptive contrast from measured group means; not a replacement for the source fitted effect. |
| `entry_class` | Richer metadata classification used to exclude contaminants and pseudo/frame entries. |
| `uniprot_acc`, `ortholog_quality` | Reference annotation and mapping confidence. |
| `analysis_status` | Protein results versus preliminary phosphosite status. |

The source LOD flag is model-wide in this export: it can reflect missingness in
groups outside the contrast being displayed. Never infer a 0/3 group from the
flag alone. Read `docs/step04_interpretation.md` for the equations and example.

## Figures and report

- `figures/step04/figure02_total_PASN_volcanoes.png`: three within-genotype PASN responses.
- `figure03_total_interaction_volcanoes.png`: three genotype-by-medium interactions.
- `figure04_total_hit_counts.png`: directional selections, split by source-model flag.
- `figure05_candidate_measurements.png`: individual measurements of the predefined candidates.
- `figure06_interaction_examples.png`: illustrative lowest-FDR conservative hit per interaction.
- `reports/step04_comparisons.md`: counts, interpretations, figures and limits.

The candidate list is the editable Step 3 file. Figure 6's examples were selected
by a declared statistical rule, not literature validation. Separate y-scales in
the sample figures allow within-protein comparison; compare points within panels.

All full data-derived tables stay in ignored `results/`. Scripts, figures,
aggregate reports and the dated notebook entry are intended GitHub contents.

Next: functional enrichment, with a declared eligible measured background. The
tables do not themselves establish pathway enrichment or direct kinase targets.
