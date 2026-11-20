# Repository guide

The group project lives in [Sarya_George/](Sarya_George/README.md). Open [Sarya_George.Rproj](Sarya_George/Sarya_George.Rproj) and run its scripts from that project directory. Paths in the examples below are relative to that directory unless stated otherwise.

## Folder map

```text
bio253-je2-proteomics/
├── README.md                    Start here
├── REPOSITORY_GUIDE.md           This guide
├── course_starter/              Original course project and teaching material
└── Sarya_George/
    ├── Sarya_George.Rproj        Open this for the group project
    ├── R/                       Shared reusable R functions
    ├── config/                  Configuration examples and local settings
    ├── data/                    Data instructions; private local inputs
    ├── results/                 Private intermediate tables and API caches
    ├── tests/                   Checks for the shared code and analysis helpers
    ├── run_step.R               Existing named-step entry point
    ├── 00_shared/               Setup, indexes and repository history
    ├── steps/
    │   ├── 01_design/
    │   ├── 02_validation_skipped/
    │   ├── 03_identifiers_resources/
    │   ├── 04_total_proteome/
    │   ├── 05_enrichment/
    │   ├── 06_string_network/
    │   ├── 07_phosphosite_results/
    │   ├── 08_mutant_venn/
    │   ├── 09_sequence_logos/
    │   └── 10_integration_poster/
    └── 90_archive/              Retained exploratory scripts
```

Within a step, `scripts/` contains executable entry scripts; `figures/` contains rendered images; `reports/` contains result tables, interpretation and run records; `notebook/` contains dated progress notes; and `docs/`, where present, contains detailed usage/methods instructions. A step only has the subfolders it uses. Cross-step functions stay in `R/` so they are shared rather than copied into several steps.

## First setup and local data

1. Open the group `.Rproj`.
2. If required, run `source("00_shared/scripts/setup/install_packages.R")`. Enrichment has the additional installer `00_shared/scripts/setup/install_enrichment_packages.R`.
3. Copy `config/data_dir.example.txt` to the ignored `config/data_dir.txt` and set the path to your local course data directory. Alternatively set `JE2_DATA_DIR`. This directory must contain `derived/` and `raw/`.
4. For the Week 2 annotation archive, copy `config/annotation_zip.example.txt` to the ignored `config/annotation_zip.txt`, or set `JE2_ANNOTATION_ZIP`.
5. Start interactive work with `source("R/setup.R")`.

The course data directory can remain outside the repository. Moving code into step folders does not require copying the measurements into those folders. See [the group data instructions](Sarya_George/data/README.md).

| Local input | Used when |
| --- | --- |
| `derived/samples.tsv` | Sample design and group labels in Step 1 and later readers. |
| `derived/protein_annotation.tsv.gz` | Protein identity/product information throughout the workflow. |
| `derived/dea_total.tsv.gz` | October total-proteome contrasts in Steps 4–7 and candidate interpretation. |
| `derived/abundances_total.tsv.gz` | Replicate detection and sample-value plots in Steps 4–6; parent-protein context. |
| `derived/dea_phospho_STYH.tsv.gz` | Phosphosite contrasts in Steps 7–9 and overlap interpretation. |
| `derived/sites_STYH.tsv.gz` | Localization probabilities, phosphosite identities and detection context in Step 7. |
| `derived/abundances_phospho_STYH.tsv.gz` | Measured phosphosite values and replicate support in Step 7. |
| Week 2 JE2 metadata/mapping archive | Identifier aliases, UniProt orthologs and supplied STRING resources in Steps 3–6. |
| Local STRING exports / cached API responses | Network import, annotation and figures in Step 6. The public report records input provenance. |
| `raw/resources/JE2_proteins.fasta` | JE2 sequence windows and protein-coordinate checks for Step 9. |

The older STY-only search and earlier DEA exports remain separate comparisons where the original scripts use them. Do not substitute them silently for the October STYH/total inputs. Raw measurements, local settings, library files and full intermediate caches remain git-ignored.

## Work through the project

| Step and location | Main code / material | Results to open |
| --- | --- | --- |
| [1 · Design](Sarya_George/steps/01_design/README.md) | Inventory and design-validation scripts; sample-design and workflow drawing scripts. | Design reports, 27-sample illustration and workflow figures. |
| [2 · Validation skipped](Sarya_George/steps/02_validation_skipped/README.md) | Records the decision to skip repeated QC; points to the course teaching notebooks. | No newly invented Step 2 analysis. |
| [3 · Identifiers/resources](Sarya_George/steps/03_identifiers_resources/README.md) | `03_build_identifier_map.R`, `03_query_web_resources.R`. | Identifier matching, aliases and public-resource/API guides. |
| [4 · Total proteome](Sarya_George/steps/04_total_proteome/README.md) | Assemble contrasts, reconcile source exports and draw total-proteome figures. | The original three interaction volcanoes, sample values and comparison reports. |
| [5 · ORA/GSEA](Sarya_George/steps/05_enrichment/README.md) | Enrichment runs, observed-only figures and module-driver extraction. | GO/KEGG tables, all significant observed-only panels and four-module protein tables. |
| [6 · STRING](Sarya_George/steps/06_string_network/README.md) | Candidate selection, UniProt upload lists, STRING import/API, network drawing. | 49/67/75-protein network figures, mappings, nodes/edges and overlap evidence. |
| [7 · Phosphosite results](Sarya_George/steps/07_phosphosite_results/README.md) | Site/parent comparison and shared labelled-volcano renderer. | Main interaction volcanoes, annotated candidate sites and localization/replicate notes. |
| [8 · Mutant Venn](Sarya_George/steps/08_mutant_venn/README.md) | Within-medium mutant Venn script. Supplementary volcanoes share the Step 7 renderer. | TSB/SN Venn diagram, counts and four labelled mutant-versus-ancestor volcanoes. |
| [9 · Sequence logos](Sarya_George/steps/09_sequence_logos/README.md) | Residue summary/Venn/logo scripts and retained structure/logo code bundle. | Exclusive-site protein list, ±10-AA sequence logos, linear protein models and 3D structures. |
| [10 · Integration/poster](Sarya_George/steps/10_integration_poster/README.md) | Integrated report, abstract preparation and methods/slide illustration scripts. | Integrated report, abstract/reference files, editable methods SVGs and presentation files. |

## Find the most requested outputs

- **Latest main interaction volcanoes:** [Step 7 `figures/volcanoes/`](Sarya_George/steps/07_phosphosite_results/figures/volcanoes/).
- **Four supplementary mutant volcanoes:** [Step 8 `figures/volcanoes/`](Sarya_George/steps/08_mutant_venn/figures/volcanoes/).
- **Volcano selection/localization/replication tables:** [Step 7 `reports/volcano_reproduction/`](Sarya_George/steps/07_phosphosite_results/reports/volcano_reproduction/README.md).
- **Every significant observed-only GO/KEGG result:** [Step 5 observed-only report](Sarya_George/steps/05_enrichment/reports/observed_only/README.md) and [figures](Sarya_George/steps/05_enrichment/figures/observed_only/).
- **Four-module protein lists:** [Step 5 module-driver report](Sarya_George/steps/05_enrichment/reports/observed_only/four_module_drivers/README.md).
- **Complete 75-protein poster network:** [Step 6 complete75 guide](Sarya_George/steps/06_string_network/reports/poster_networks/complete75_README.md).
- **IsaA/Pyk/QoxB overlap investigation:** [Step 6 phosphosite-overlap report](Sarya_George/steps/06_string_network/reports/phosphosite_overlap/README.md).
- **Sequence, protein and structure PDFs:** [Step 9 structures and logos](Sarya_George/steps/09_sequence_logos/reports/structures_and_logos/README.md).
- **Editable poster methods illustrations:** [Step 10 methods SVGs](Sarya_George/steps/10_integration_poster/figures/methods_svg/README.md).
- **Abstract draft, primary literature and RIS references:** [Step 10 abstract folder](Sarya_George/steps/10_integration_poster/reports/abstract/).

`all_models` volcano files retain the source-filtered observed and LOD-imputed estimates. `observed_only` removes LOD-imputed rows. `JE2_observed_only` also restricts to ordinary JE2 protein entries. These are different plotted populations, identified by their filenames and the accompanying report. Labelling IsaA S194, Pyk T537 and QoxB T552 does not make their interaction tests significant.

The earlier Step 7 volcano and the original Step 8/9 figures remain alongside the later poster versions. The reorganization does not replace earlier results with a new analysis.

## Shared and historical folders

| Folder | What belongs there |
| --- | --- |
| [R/data_access/](Sarya_George/R/data_access/) | Input discovery, parsing and table readers. |
| [R/identifiers/](Sarya_George/R/identifiers/) | Identifier/alias mapping and lookup. |
| [R/web_resources/](Sarya_George/R/web_resources/) | UniProt, STRING and KEGG client helpers. |
| [R/analysis/](Sarya_George/R/analysis/) | Reusable design, comparison, enrichment and parent-protein functions. |
| [R/figures/](Sarya_George/R/figures/) | Reusable figure functions. |
| [R/reference/](Sarya_George/R/reference/) | Design constants and reference helpers. |
| [config/](Sarya_George/config/) | Shareable examples and candidate configuration; personal path settings are ignored. |
| [data/](Sarya_George/data/README.md) | Local-data instructions and private inputs, not published measurements. |
| [results/](Sarya_George/results/README.md) | Full local intermediate tables, per-sample data and caches; reports contain the shareable summaries. |
| [tests/](Sarya_George/tests/) | Parsing/data tests and focused checks. |
| [00_shared/](Sarya_George/00_shared/README.md) | Installers, code/provenance indexes, original group README and repository progress. |
| [90_archive/](Sarya_George/90_archive/README.md) | Earlier exploratory scripts retained with their original context. |
| [course_starter/](course_starter/README.md) | Original teaching project: its own `R/`, `analyses/`, `scripts/`, `tests/`, data instructions and `.Rproj`. |

Use the group project for the group results. Open `course_starter/bio253-je2-proteomics.Rproj` when working specifically through the original teaching notebooks. Its `analyses/01_explore_data.qmd`, `03_mutant_validation.qmd` and `04_phospho_two_roadmaps.qmd` remain available as course guidance.
