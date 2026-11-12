# Run all tests:   Rscript tests/run_tests.R      (from the repository root)
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
library(testthat)
test_dir(file.path(REPO_ROOT, "tests", "testthat"), reporter = "summary",
         stop_on_failure = TRUE)
