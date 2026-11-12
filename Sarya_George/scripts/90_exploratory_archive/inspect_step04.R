# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source('R/setup.R')
for (dataset in c('total','phospho_STYH')) {
 x <- read_dea(dataset, drop_contaminants=FALSE)
 cat('\nDATASET', dataset, '\n'); print(names(x)); print(count(x,estimate_type))
 print(x |> filter(contrast %in% c('TSB_vs_SN_at_anc','interaction_stk_vs_anc_at_TSB_vs_SN')) |> group_by(contrast,estimate_type) |> summarise(n=n(),nonfinite_diff=sum(!is.finite(diff)),nonfinite_FDR=sum(!is.finite(FDR)),.groups='drop'))
 ab <- read_abundances(dataset, drop_contaminants=FALSE)
 cat('ABUNDANCE COLUMNS\n'); print(names(ab))
}
