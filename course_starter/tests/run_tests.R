# Run all tests:   Rscript tests/run_tests.R      (from the repository root)
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
library(testthat)
test_dir(file.path(REPO_ROOT, "tests", "testthat"), reporter = "summary",
         stop_on_failure = TRUE)
