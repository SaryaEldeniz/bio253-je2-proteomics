# Repository reorganization - 20 November 2026

The project is organized by research step. Start with the [repository guide](../../../REPOSITORY_GUIDE.md) or the [group step index](../../README.md).

## What changed

- Moved 797 of the 839 original nonignored files into the organized layout. Every original file is retained.
- Grouped scripts, figures, result reports, methods documents and dated notes under `steps/01_design` through `steps/10_integration_poster`.
- Collected the teammate's phosphosite, Venn, sequence-logo and structural work in Steps 7, 8 and 9 alongside the related figures and tables.
- Placed the three main interaction volcanoes in Step 7 and the four supplementary mutant-versus-ancestor volcanoes in Step 8. PNG, PDF and SVG files are retained, including observed-only and all-model variants.
- Kept reusable `R/` helpers, configuration, tests, the RStudio project and private `data/` / `results/` locations at the group root.
- Moved original course templates and the original course README into `course_starter/`. Retained the original group README in `00_shared/project_history/` and exploratory code in `90_archive/`.
- Added a repository guide, group index and a guide for every step, including the explicitly skipped Step 2.
- Updated path references with the user's authorization. This includes R entry points, relative report links, output folders and the Step 9 Python directory-depth calculation. Analysis definitions, thresholds, data selection and figure appearance were not changed.
- Preserved the source measurements, private configuration, local libraries and complete intermediate tables in their existing ignored locations. No raw course data were added to Git.

## Preservation and verification

`file_moves_and_preservation.csv` lists every original path, new path, pre-move and post-update SHA-256 hash, and whether its contents are unchanged or its path references were updated. All original figures and other binary artifacts remain byte-for-byte unchanged.

The path updater initially encountered Windows text-encoding problems. Affected files were restored from verified original backups before UTF-8-safe path updates were applied. Preservation and syntax checks were repeated after restoration.

Verification included all 108 original R files: every file parses and its numeric literals are unchanged. The original helper tests and focused design, identifier, web-resource, Step 4, ORA, GSEA and saved enrichment-result checks passed. Local Markdown links were checked after relocation. No models or source FDR values were refitted for the reorganization.

`validation_artifacts/` retains a report created at a former path by an initial validation run. Only empty former directories were removed; no original files were discarded. Existing RIS/session-information whitespace was preserved.

Run the project from `Sarya_George/Sarya_George.Rproj`. The new folder layout and exact rerun commands are documented in each step's README.
