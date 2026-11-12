# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
contrasts <- c("interaction_stk_vs_anc_at_TSB_vs_SN",
               "interaction_stp_vs_anc_at_TSB_vs_SN")
d <- read_dea("total", contrasts) |>
  filter(locus_tag %in% c("B7H15_06050", "B7H15_02095", "B7H15_08530", "B7H15_14250")) |>
  select(locus_tag, product, contrast, diff, FDR, estimate_type)
write_tsv(d, "results/exploratory/I_marker_interactions.tsv")
print(d, n = 20, width = 170)
