# Shared setup, indexes and history

These materials apply across the workflow. Step-specific work is in [the numbered steps](../README.md).

| Folder | Contents |
| --- | --- |
| [scripts/setup/](scripts/setup/) | General/enrichment package installers and the R code-index writer. |
| [docs/](docs/) | `R_code_index.tsv` and `code_provenance.tsv` for locating code and its origins. |
| [notebook/](notebook/) | Repository-organization/publication history and the original work-log README. |
| [reports/](reports/) | General project progress. |
| [project_history/](project_history/README.md) | The original group README, retained as historical project documentation. |

Shared executable functions remain in [../R/](../R/); local path settings in [../config/](../config/); private measurements/caches are described in [../data/](../data/README.md) and [../results/](../results/README.md); checks remain in [../tests/](../tests/).

Run setup scripts from the group `.Rproj`, for example `source("00_shared/scripts/setup/install_packages.R")`. Follow the [current repository guide](../../REPOSITORY_GUIDE.md) for current navigation.
