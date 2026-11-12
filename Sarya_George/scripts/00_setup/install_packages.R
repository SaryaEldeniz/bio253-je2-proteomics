# Run once in R/RStudio. Installation is separate from analysis scripts.
required <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "testthat", "readxl", "curl", "jsonlite", "digest")
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Required packages are available.")
