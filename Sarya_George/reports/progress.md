# Project I progress

Updated 11 November 2026.

## Completed: Step 1, design and comparison plan

- Created a separate code project and initialized a local Git repository on main.
- Reused the lecturer's MIT-licensed helpers with attribution. Adapted data-path
  configuration so each teammate can read their own local course-data directory.
- Checked 27 unique sample labels, nine genotype-medium groups and three replicate
  labels per group in both measurement layers.
- Confirmed nine shortened pulse sample names resolve to canonical IDs.
- Verified two phospho-specific reinjection flags: JE2stk_pulseSN_1 and _2.
- Recorded 33 available contrasts per layer and 12 planned Project I contrasts.
- Documented A-minus-B direction and explicit group-mean coefficients.
- Inventoried abundance versus DEA coverage: total 2348 versus 2324 features;
  STYH phospho 4278 versus 4235 sites. Saved omitted identities locally; no cause
  for omission is asserted at this stage.
- Generated Figure 1 directly from the layer metadata and visually checked it.
- Saved input hashes and R/package versions.
- Independent code review checked scientific direction, metadata handling and
  portability. Review led to stronger coefficient/contrast checks and allowing
  helper tests to load without private input data.

## Validation status

| Check group | Result |
|---|---|
| Actual-input design checks | 28/28 passed |
| Focused synthetic contrast/name checks | 10/10 passed |
| Same focused checks without private data | 10/10 passed |
| Full original starter testthat suite | Unavailable: testthat is not installed |

`source("scripts/00_setup/install_packages.R")` installs declared dependencies when network access
is available. Then run `Rscript tests/run_tests.R` for the full original suite.
That suite has not been reported as passing. R emitted locale and package-build
version warnings; the successful checks and figure generation still exited 0.

## Outputs and code to read

- `scripts/01_data/01_validate_design.R`: main Step 1 R script and report generation.
- `R/analysis/design.R`: comparison plan, coefficient functions and design figure.
- `R/setup.R`: load packages/helpers and locate local data.
- `tests/check_design.R`: focused checks using synthetic group means/names.
- `reports/step01_design.md`: actual findings and planned comparisons.
- `reports/step01_inputs.md`: input file hashes and sizes.
- `reports/step01_session.txt`: R/package versions.
- `reports/step01_code_checks.md`: focused test results.
- `results/step01/`: local sample, contrast, validation and feature-coverage tables.
- `figures/figure01_design.png`: experimental-design figure.
- `notebook/2026-11-11_step01.md`: dated record of decisions and AI help.

## Remaining, in order

- Step 2: skipped at the user's request; the user reports the data are already validated.
- Step 3: completed; see the section below.
- Step 4: total-proteome PASN comparisons and genotype interactions.
- Step 5: ORA and ranked functional enrichment.
- Step 6: module/network figures.
- Step 7: phosphosite reliability and parent-protein checks.
- Step 8: mutant Venn diagrams with directions.
- Step 9: residue counts and site-centred sequence logos.
- Step 10: module/site/literature integration.
- Step 11: optional pulse extension and final reproducible report.

## GitHub status

The local repository is initialized. A remote repository and upload are future
steps. `.gitignore` excludes original inputs, the local data-location setting and
detailed local result tables. Code, selected figures, aggregate reports and the
notebook are the intended version-controlled content. No commits or uploads have
been made during Step 1.

## Completed: Step 3, local identifiers and public APIs

- Built a master map of 4,189 source entries, keeping exclusions/unknowns in an
  audit and defining ordinary protein entries by metadata `entry_class == "cds"`.
- Built 55,076 searchable alias/name/description records, including explicit
  Stk1/PknB/Stp1 names, exact JE2 keys, locus tags, UniProt and local STRING IDs.
- Total DEA: 2,229 ordinary CDS entries; 2,107 good JE2-to-UniProt mappings,
  corresponding to 2,080 unique target accessions. All 2,229 have a local STRING
  alias; 2,143 have terms and 1,640 have local physical-network edges.
- Phospho DEA: 1,109 ordinary parent-protein entries; 1,068 good mappings,
  corresponding to 1,063 unique target accessions. These are parent counts, not
  phosphosite counts.
- The richer metadata identified 68 total-DEA contaminant entries whose DEA flag
  was FALSE, alongside 24 separately flagged contaminants and three frame/pseudo
  entries. The master audit keeps them; primary ordinary-protein analysis excludes them.
- Added local partner lookup and reference-edge labelling. Many-to-one matches
  stay explicit, and edge rows are not multiplied by adding labels.
- Added cached, rate-limited UniProt, STRING and KEGG R clients, with response
  provenance and unresolved/ambiguous lookup handling.
- Network access was initially restricted; after a turn-scoped permission grant,
  live requests succeeded. No experimental measurements were sent.
- Live demonstration: five exact UniProt primary accessions, five STRING IDs
  verified against UniProt cross-references, four selected functional edges, four
  selected physical edges, eight public-context partner rows, 194 annotation rows,
  four KEGG gene-pathway links and one retrieved KEGG gene entry.
- Passed five actual mapping checks, 12 synthetic identifier/edge checks and
  17 offline API checks. The full original testthat suite remains unavailable
  because testthat is not installed; it was attempted again and recorded honestly.

### Step 3 code and outputs

- `scripts/02_identifiers/03_build_identifier_map.R`: build all local annotation/alias tables.
- `scripts/03_web_resources/03_query_web_resources.R`: reproduce the live public-identifier example.
- `R/identifiers/identifiers.R`: lookup_je2, prepare_uniprot_query, resource links, local
  interactions and reference-edge labelling.
- `R/web_resources/web_resources.R`: uniprot_entry, string_resolve, string_network,
  string_partners, string_functional_annotation, kegg_links and kegg_entry.
- `config/interesting_proteins.tsv`: editable starting protein list.
- `docs/step03_identifiers.md` and `docs/web_resources.md`: interactive examples.
- `reports/step03_identifiers.md`: mapping denominators, exclusions and example IDs.
- `reports/step03_live_api.md`: verified live retrieval summary.
- `results/step03/`, `results/step03_api/`, `results/api_cache/`: local detailed
  results and caches; excluded from Git and the shareable code archive.

## Completed: Step 4, annotated comparisons and total-proteome figures

- Assembled the ten nonoptional comparisons: three within-genotype PASN effects,
  four mutant-versus-ancestor effects in TSB/SN, and all three interactions.
- Preserved original effect estimates, statistics and FDRs. Retained all entries
  for audit and used metadata CDS classification for the ordinary-protein results.
- Added explicit direction labels, two/four relevant-group observation counts,
  source-model flags and a declared conservative sensitivity subset.
- Total proteome: 2,229 ordinary protein entries per comparison. The ancestor,
  kinase and phosphatase within-genotype PASN comparisons select 443, 478 and
  446 entries with source FDR <0.1 and |effect| >1.
- The three interactions select 95, 88 and 118 protein entries; conservative
  selected subsets contain 67, 63 and 86 entries, respectively.
- The phosphosite layer has 4,100 CDS-parent site rows per comparison. Tables are
  preliminary: localization filtering and parent-protein evaluation remain Step 7.
- Important correction to the earlier shorthand about estimate_type: its LOD
  label is model-wide in these exports. In the ancestor comparison, 25 ordinary
  protein entries with this flag have 3/3 actual measurements in both groups.
  A flag does not identify an absent group within the displayed comparison.
- Generated and visually checked five figures: PASN volcanoes, interaction
  volcanoes, directional counts, predefined candidate measurements and illustrative
  interaction examples. No missing measurements were replaced by zero.
- Passed 20/20 input checks, including preserved source values and the third
  interaction's effect identity; passed 8/8 focused code checks. Tested the
  read_comparison RStudio helper on the actual output tables.
- Independent review checked signs, filtering, flag semantics and table units.
  The full original helper suite remains unavailable because testthat is missing.

### Step 4 code and outputs

- `scripts/04_analysis/04_assemble_comparisons.R`: regenerate tables, figures and report.
- `R/analysis/step04_results.R`: annotated selections, observation support and table reader.
- `R/figures/step04_plots.R`: figure functions.
- `tests/check_step04.R`: synthetic threshold/sign/missingness checks.
- `reports/step04_comparisons.md`: findings, figures and limits.
- `docs/step04_usage.md`, `docs/step04_interpretation.md`: RStudio use and interpretation.
- `results/step04/total/` and `results/step04/phospho_STYH/`: complete local tables,
  per-contrast selections, directional lists, conservative subsets and audits.
- `figures/step04/`: five figures generated by the R code.

Next: Step 5 functional enrichment, using the annotated ordinary-protein results
and a declared eligible measured background. No enriched modules have yet been inferred.
