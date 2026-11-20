#!/usr/bin/env Rscript
# check_data.R -- verify the course data files against data/checksums.md5.
# Same check as check_data.sh, for Windows or anyone without a shell.
#
# Usage (from the repository root):  Rscript data/check_data.R
# Exit status 0 = all files present and intact, 1 = something is missing or corrupt.

data_dir <- if (file.exists("data/checksums.md5")) "data" else "."
sums_file <- file.path(data_dir, "checksums.md5")
if (!file.exists(sums_file)) stop("checksums.md5 not found. Run this from the repository root.")

lines <- readLines(sums_file)
lines <- lines[nzchar(trimws(lines))]
expected <- sub("^([0-9a-f]{32}) +.*$", "\\1", lines)
files    <- sub("^[0-9a-f]{32} +", "", lines)
paths    <- file.path(data_dir, files)

present <- file.exists(paths)
actual  <- rep(NA_character_, length(paths))
actual[present] <- unname(tools::md5sum(paths[present]))

status <- ifelse(!present, "MISSING", ifelse(actual == expected, "OK", "CORRUPT"))
for (i in seq_along(files)) cat(sprintf("%-8s %s\n", status[i], files[i]))

cat(sprintf("\n%d OK, %d missing, %d corrupt\n",
            sum(status == "OK"), sum(status == "MISSING"), sum(status == "CORRUPT")))
if (any(status != "OK")) {
  cat("Data check FAILED. Extract the course data zip into data/ (so that data/raw/, data/derived/, data/reference/ exist) and run again.\n")
  quit(status = 1)
}
cat("Data check passed. Next: open analyses/01_explore_data.qmd\n")
