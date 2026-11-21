# Step 8 · Mutant Venn diagrams and supplementary volcanoes

These comparisons test `stk vs ancestor` and `stp vs ancestor` separately within TSB and SN. They address within-medium genotype differences. The three genotype-by-medium interaction tests are in [Step 7](../07_phosphosite_results/README.md).

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | `08_render_venn.R`, the original overall site-overlap renderer. |
| [reports/](reports/) | `step08_phospho_venn_counts.tsv`, retaining shared same/opposite-direction counts. |
| [figures/](figures/) | Original overall STYH Venn diagram. |
| [figures/volcanoes/](figures/volcanoes/) | Four labelled supplementary comparisons, combined and individual, as PNG/PDF/SVG. |
| [notebook/](notebook/) | Selection and overlap notes. |

Required local input: `derived/dea_phospho_STYH.tsv.gz`. Supplementary volcano annotations also use the localization, replicate and parent-protein inputs listed in Step 7.

```r
source("steps/08_mutant_venn/scripts/08_render_venn.R")
# Rebuild both the main and supplementary labelled volcanoes:
source("steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R")
```

The original overall Venn includes source-observed and LOD-imputed hits meeting the stated thresholds. The residue-specific exclusive observed Venn is in [Step 9](../09_sequence_logos/README.md). “Only stk/stp” means that a site met selection in only one mutant-versus-ancestor comparison within that medium; it does not mean the site was detected in only that strain.

The `JE2_observed_only` supplementary volcanoes show the observed JE2 population used for the exclusive-site interpretation. All three named sites are labelled in all four panels, including nonsignificant comparisons. Full annotated estimates remain in the [Step 7 report](../07_phosphosite_results/reports/volcano_reproduction/README.md).

The new [compact poster volcanoes](figures/volcanoes/compact_poster/) add the JE2 LOD-imputed estimates as crosses alongside the observed circles. Combined and individual SVG/PNG/PDF exports retain the three named-site annotations. This display change does not change the observed-only exclusive-site selection or the Venn diagrams. See the [matched volcano guide](../10_integration_poster/reports/compact_volcanoes/README.md) for counts, the shared legend, all figure links and the R reproduction command.

[Next: residues and sequence logos](../09_sequence_logos/README.md) · [Project overview](../../README.md)
