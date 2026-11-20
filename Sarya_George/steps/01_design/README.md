# Step 1 · Experimental design and input inventory

Use this folder to understand the 27 samples: three genotypes × three media × three replicates. The original source names, contrast directions and design notes are retained here.

| Subfolder | What to find |
| --- | --- |
| [scripts/](scripts/) | `00_inventory_inputs.R`, `01_validate_design.R`, `01_render_design_summary.R` and `01_render_project_workflow.R`. |
| [reports/](reports/) | Input inventory, experimental-design explanation, code checks and session record. |
| [figures/](figures/) | Original design plot; `design/` contains the 27-sample/interactions table and workflow PNG/PDF/SVG files. |
| [notebook/](notebook/) | Original design-check R script, experimental-design text and dated notes. |

Required local data: the configured course data directory, particularly `derived/samples.tsv` and the source input files being inventoried. The configuration is explained in [the repository guide](../../../REPOSITORY_GUIDE.md).

From the group project, the design entry point is:

```r
source("steps/01_design/scripts/01_validate_design.R")
```

[Next: Step 2 decision](../02_validation_skipped/README.md) · [Project overview](../../README.md)
