# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
files <- c(file.path(DATA_ROOT, "raw/total/DE_total.xlsx"), file.path(DATA_ROOT, "raw/phospho_STYH/DE_phospho_STYH.xlsx"))
for (f in files) {
  cat("\nSOURCE", f, "\n")
  print(readxl::excel_sheets(f))
  x <- readxl::read_excel(f, sheet = "diff_exp_analysis", n_max = 1)
  cat("DEA columns:", paste(names(x), collapse = " | "), "\n")
}
