# 20 November 2026 - Repository reorganization

## Request and authorization

Publish the current phosphosite figures and organize both teammates' code, figures, tables and notes by step. The user authorized path-reference changes only so the relocated files remain usable. No original files were to be deleted and no analyses were to be changed.

## Work completed

Created a SHA-256 inventory of all 839 original nonignored files before any moves. Moved 797 files according to an explicit old-to-new map, then verified preservation before editing paths. Updated executable and document path references, including the shared named-step launcher, main/supplementary figure output directories and the Step 9 Python bundle. Added navigation guides and the detailed file-move/preservation manifest.

Original course material is in `course_starter/`; group work is in `Sarya_George/steps/`. Step 7 contains the main interaction volcanoes, Step 8 contains the supplementary mutant comparisons and Venn diagrams, and Step 9 contains sequence logos and structural context. Shared helpers and private local data/results retain their established group-root locations.

## Checks and limits

Every original is retained. Existing binary artifacts are unchanged. All 108 original R sources parse, with unchanged numeric literals. Helper and focused analysis checks passed, including saved observed-only enrichment and figure-membership checks. Local Markdown links resolve. The original analysis outputs were relocated rather than recomputed.

A UTF-8 issue in the path-update utility was caught during verification. The affected text was restored from verified backups, the utility was corrected, and preservation and syntax checks were repeated before publication preparation.

## AI assistance

Codex inventoried and reorganized the repository, updated paths and wrote the guides. Read-only agents audited path dependencies and preservation; a separate agent drafted new navigation documents. All file moves and final hashes are recorded in `00_shared/reorganization/file_moves_and_preservation.csv`.
