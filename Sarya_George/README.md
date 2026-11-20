# Sarya and George · Project I

**Which functional modules respond to PASN, and where do Stk1/Stp1 and their candidate targets sit?**

Open [Sarya_George.Rproj](Sarya_George.Rproj) in RStudio. The [complete repository guide](../REPOSITORY_GUIDE.md) explains configuration, local data and every shared folder.

| Step | Open here |
| --- | --- |
| 1 | [Experimental design and input inventory](steps/01_design/README.md) |
| 2 | [Validation skipped: decision and course references](steps/02_validation_skipped/README.md) |
| 3 | [Identifier matching, annotations and web resources](steps/03_identifiers_resources/README.md) |
| 4 | [Total-proteome comparisons and volcanoes](steps/04_total_proteome/README.md) |
| 5 | [GO/KEGG ORA, GSEA and module contributors](steps/05_enrichment/README.md) |
| 6 | [Candidate proteins and STRING networks](steps/06_string_network/README.md) |
| 7 | [Reliable phosphosite interpretation and main interaction figures](steps/07_phosphosite_results/README.md) |
| 8 | [Mutant Venn diagrams and supplementary volcanoes](steps/08_mutant_venn/README.md) |
| 9 | [Residues, sequence logos and protein/structure figures](steps/09_sequence_logos/README.md) |
| 10 | [Integrated report, abstract and poster material](steps/10_integration_poster/README.md) |

Each step groups its scripts, figures, result tables and dated notes. Shared helpers stay in [R/](R/), with [configuration](config/), [private data instructions](data/README.md), [local intermediate results](results/README.md) and [tests](tests/) at the project level.

From this project directory, start with:

```r
source("R/setup.R")
lookup_je2(c("Stk1", "Stp1", "TcyA"))
```

The latest labelled interaction and supplementary volcanoes are rebuilt together:

```r
source("steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R")
```

Configure your local course data directory using `config/data_dir.txt` or `JE2_DATA_DIR`, and the Week 2 metadata archive using `config/annotation_zip.txt` or `JE2_ANNOTATION_ZIP`. Personal path files and raw measurements are not uploaded.

The [shared setup/history folder](00_shared/README.md), [retained exploratory archive](90_archive/README.md) and [original group README](00_shared/project_history/README.md) preserve the earlier project context.
