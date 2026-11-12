# Sarya_George

R code and progress for Sarya and George's Project I: **Which functional modules
respond to PASN, and where do Stk1/Stp1 and their candidate targets sit?**

## Organization

| Folder | Contents |
|---|---|
| `R/data_access/` | File discovery, input provenance, table extraction and reading helpers. |
| `R/identifiers/` | JE2/UniProt/STRING matching, aliases and annotation lookup. |
| `R/web_resources/` | Cached UniProt, STRING and KEGG API clients. |
| `R/analysis/` | Experimental design, comparison logic, ORA and ranked GSEA. |
| `R/figures/` | Reusable R functions that draw figures. |
| `R/reference/` | Course constants, reference lookup and pHis helper functions. |
| `scripts/00_setup/` | Package installation script. |
| `scripts/01_data/` | File inventory and design-validation entry scripts. |
| `scripts/02_identifiers/` | Build the local annotation/identifier map. |
| `scripts/03_web_resources/` | Run the small live public-API demonstration. |
| `scripts/04_analysis/` | Assemble comparison tables and reconcile August/October counts. |
| `scripts/05_enrichment/` | ORA/GSEA runs, sensitivity analyses and the generated results guide. |
| `scripts/05_figures/` | Regenerate figures from saved analysis tables. |
| `scripts/90_exploratory_archive/` | Earlier exploratory R scripts, explicitly separated from final analysis. |
| `tests/` | Focused synthetic checks and the original helper tests. |
| `docs/`, `reports/`, `notebook/` | Instructions, aggregate findings and dated progress. |
| `figures/` | Generated images; the code that creates them is in the folders above. |

The existing `notebook/1_CheckDesign.R` and original notebook documents are retained.
`docs/code_provenance.tsv` maps existing code files to their organized locations.
`docs/R_code_index.tsv` lists every R file and its top-level functions.

## Open and configure

Open **`Sarya_George/Sarya_George.Rproj`** in RStudio. If packages are missing, run
`source("scripts/00_setup/install_packages.R")` from this folder.

Set the source data directory in the ignored `config/data_dir.txt` (copy the
example first), or set `JE2_DATA_DIR`. It must contain `derived/` and `raw/`.
Without an override, the code looks in the repository's root `data/` directory.
Configure the Week 2 metadata zip in `config/annotation_zip.txt`, or set
`JE2_ANNOTATION_ZIP`. No machine-specific path settings are uploaded.

## Run the workflow

From the group project in RStudio:

```r
source("run_step.R")
run_project_step("inventory")
run_project_step("design")
run_project_step("identifiers")
run_project_step("analysis")
run_project_step("figures")   # redraw without rebuilding comparison tables
run_project_step("enrichment")
run_project_step("enrichment_figures")
run_project_step("enrichment_resources") # public-ID queries; network required
run_project_step("enrichment_report")
```

Optional: `run_project_step("api")` queries public identifiers only; it needs
network access. `run_project_step("reconcile")` compares the August and October
exports when the older course workbook is available.

From the repository root, use `Rscript Sarya_George/run_step.R identifiers` or
another step name. The entry scripts support the repository root and group root.

Step 2 sample QC was skipped at the user's request because it was already
validated. The main tables use the October source export; the lecturer screenshot
matches August. The full explanation and source-model flag nuance are in the
Step 4 reports. No new source models are fitted by this pipeline.

## Interactive lookup

```r
source("R/setup.R")
lookup_je2(c("Stk1", "Stp1", "TcyA"))
lookup_je2("cystine", exact=FALSE)
local_interactions(c("Stk1", "Stp1"))$edges
uniprot_entry("Q2FZ64")$summary
hits <- read_comparison("interaction_stk_vs_anc_at_TSB_vs_SN", conservative_only=TRUE)
View(hits)
```

See `docs/step03_identifiers.md`, `docs/web_resources.md` and `docs/step04_usage.md`.

## Step 5 pathway analysis

Install the additional enrichment engine once with
`source("scripts/00_setup/install_enrichment_packages.R")`. Compatible R 4.4
Windows binaries are installed in the ignored project-local `renv/library/`.
The general installer also uses this local folder. This is a local package
library, not a complete renv lockfile; exact run versions are recorded in the report.

Open [the Step 5 results guide](reports/step05/README.md) for findings, pathway
tables and figures, and [the methods guide](docs/step05_methods.md) for every
analysis choice. `read_enrichment("GSEA", collection="KEGG")` loads primary
study-significant findings; set `significant_only=FALSE` for all tested terms.

ORA uses directional hit lists with source FDR <0.1 and absolute log2 effect >1.
GSEA uses all eligible proteins ranked by the signed source statistic. The main
background has 2,058 singleton UniProt ortholog IDs. KEGG, GO BP/MF/CC and UniProt
keywords come from the fixed supplied STRING-v12 snapshot; live reference-organism
lookups do not change those memberships. Replication, model-label and local-JE2
sensitivities remain separately labelled. Main pathway significance is study-wide
BH FDR <0.1 by method across all ten comparisons and five collections.

This step concerns protein abundance. Phosphosite localization, protein-adjusted
site changes, Venn diagrams and sequence logos remain in the later site workflow.

## Checks and provenance

Run the `tests/check_design.R`, `tests/check_identifiers.R`,
`tests/check_web_resources.R` and `tests/check_step04.R` scripts from the group
folder, or prefix their paths with `Sarya_George/` from the repository root.
These focused checks need no private data or API calls. The original full helper
suite additionally requires `testthat` and uses the course data when available.

The original course helpers are copyright 2026 Jonas Grossmann / FGCZ, reused
under the repository's MIT license. New project code and provenance are described
in the manifests and dated notebook entries. Historical exploratory scripts may
use superseded thresholds, hard-coded demonstration captions or optional older
inputs; they are not the current final inference.

## Local data and outputs

Source measurements, metadata archives, instrument files, API response caches,
full result tables, machine paths and RStudio caches stay out of Git. The code
generates detailed tables under ignored `Sarya_George/results/`; no original
measurements were copied into this project folder. Selected figures and aggregate
reports are included for review and can be regenerated.
