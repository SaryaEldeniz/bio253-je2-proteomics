# install_packages.R -- run once:  Rscript install_packages.R
# Needs R >= 4.3 and Quarto (https://quarto.org) for rendering the notebooks.

cran <- c(
  "tidyverse",   # dplyr, tidyr, readr, stringr, ggplot2, purrr, tibble
  "readxl",      # read the prolfquapp Excel files (only needed by 00_fetch_data.R)
  "ggrepel",     # non-overlapping labels in volcano plots (optional)
  "patchwork",   # combine plots (optional)
  "testthat",    # tests/
  "knitr", "rmarkdown", "quarto",  # rendering the notebooks
  "BiocManager"
)
missing <- setdiff(cran, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)

bioc <- c("limma")   # moderated t-tests
missing <- setdiff(bioc, rownames(installed.packages()))
if (length(missing) > 0) BiocManager::install(missing, update = FALSE, ask = FALSE)

# Optional: prolfqua / prolfquapp, the packages used for the lecturer's DEA.
# Not needed for the exercises. Install if you want to re-run the DEA yourself:
#   install.packages("remotes")
#   remotes::install_github("fgcz/prolfqua")
#   remotes::install_github("prolfqua/prolfquapp")

message("Done. Check: Rscript tests/run_tests.R")
