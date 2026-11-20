# setup.R -- source this ONE file at the top of every script / notebook:
#
#   source("R/setup.R")        # from the repository root
#   source("../R/setup.R")     # from inside analyses/
#
# It loads the packages we use everywhere, finds the repository root
# (so paths work no matter where you run from) and sources all helper
# functions in R/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# Walk up from the current folder until we find the folder that contains
# both CLAUDE.md and R/setup.R. That folder is the repository root.
find_repo_root <- function(start = getwd()) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "CLAUDE.md")) &&
        file.exists(file.path(dir, "R", "setup.R"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (parent == dir) {
      stop("Could not find the repository root. ",
           "Run this from inside the bio253-je2-proteomics folder.")
    }
    dir <- parent
  }
}

REPO_ROOT <- find_repo_root()

# Source every helper file in R/ except this one.
for (f in sort(list.files(file.path(REPO_ROOT, "R"), pattern = "\\.R$", full.names = TRUE))) {
  if (basename(f) != "setup.R") source(f)
}
rm(f)

# One consistent plot style for the whole course.
theme_set(theme_bw(base_size = 11))
