#!/usr/bin/env Rscript
# 00_fetch_data.R -- copy the course data into data/raw/ and write tidy tables
# into data/derived/.
#
# Usage (from the repository root):
#   Rscript scripts/00_fetch_data.R                    # uses $BIO253_DATA_DIR or the default path
#   Rscript scripts/00_fetch_data.R /path/to/SAJE2     # explicit path
#   Rscript scripts/00_fetch_data.R /path/to/SAJE2 --link   # symlink instead of copy
#
# /path/to/SAJE2 is the data folder you got from the lecturer. It contains
# DEAs/, phospho_combined_sites/ and more_resources/.
#
# The data are UNPUBLISHED. data/raw/ and data/derived/ are git-ignored:
# never commit them, never upload them anywhere.

# Default: the course folder's SAJE2/, two levels above the repository root
# (lecturer layout: <course>/gh_starter_pack_repos/bio253-je2-proteomics/).
DEFAULT_SOURCE <- file.path("..", "..", "SAJE2")

args <- commandArgs(trailingOnly = TRUE)
use_links <- "--link" %in% args
args <- setdiff(args, "--link")
source_dir <- if (length(args) >= 1) args[1] else Sys.getenv("BIO253_DATA_DIR", DEFAULT_SOURCE)

# Load helpers (parse_* functions). This works when run from the repo root
# or from scripts/.
setup <- if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R"
source(setup)

if (!dir.exists(source_dir)) {
  stop("Data folder not found: ", source_dir, "\n",
       "Pass it as an argument or set the environment variable BIO253_DATA_DIR.")
}
message("Reading course data from: ", source_dir)

# ---- 1. Which files do we need? ---------------------------------------------
# left: name inside data/raw/   right: path inside the course data folder
files <- c(
  "total/DE_total.xlsx" =
    "DEAs/DEA_20261029_WUtotal_JE2_Oct26_vsn/Results_WU_total_JE2_Oct26/DE_WUtotal_JE2_Oct26.xlsx",
  "phospho_STYH/DE_phospho_STYH.xlsx" =
    "DEAs/DEA_20261029_WUphospho_singleSite_PTMstyh_vsn/Results_WU_phospho_singleSite_PTMstyh/DE_WUphospho_singleSite_PTMstyh.xlsx",
  "phospho_STY/DE_phospho_STY.xlsx" =
    "DEAs/DEA_20260801_WUphospho_singleSite_PTMsty_vsn/Results_WU_phospho_singleSite_PTMsty/DE_WUphospho_singleSite_PTMsty.xlsx",
  "sites/combined_site_STY.tsv"  = "phospho_combined_sites/combined_site_STY_79.9663.tsv",
  "sites/combined_site_STYH.tsv" = "phospho_combined_sites/combined_site_STYH_79.9663.tsv",
  "resources/JE2_proteins.fasta" =
    "more_resources/p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta",
  "resources/JE2_genome.gff" = "more_resources/SE_SE2_GCA_002085525.1_ASM208552v1_genomic.gff",
  "resources/string_aliases.txt" =
    "more_resources/stringDBresources/STRG0A35SHG.protein.aliases.v12.0.txt",
  "resources/string_enrichment_terms.txt" =
    "more_resources/stringDBresources/STRG0A35SHG.protein.enrichment.terms.v12.0.txt",
  "resources/string_physical_links.txt" =
    "more_resources/stringDBresources/STRG0A35SHG.protein.physical.links.v12.0.txt"
)

raw_dir <- file.path(REPO_ROOT, "data", "raw")
derived_dir <- file.path(REPO_ROOT, "data", "derived")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Copy (or link) the raw files ----------------------------------------
manifest <- list()
for (target in names(files)) {
  src <- file.path(source_dir, files[[target]])
  dst <- file.path(raw_dir, target)
  if (!file.exists(src)) stop("Expected file is missing in the course data: ", src)
  dir.create(dirname(dst), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(dst) || nzchar(Sys.readlink(dst))) unlink(dst)
  ok <- if (use_links) file.symlink(normalizePath(src), dst) else file.copy(src, dst)
  if (!ok) stop("Could not copy ", src)
  manifest[[target]] <- tibble(file = target, source = files[[target]],
                               size_MB = round(file.size(src) / 1e6, 1),
                               md5 = unname(tools::md5sum(src)))
  message("  ", if (use_links) "linked " else "copied ", target)
}
write_tsv(bind_rows(manifest), file.path(raw_dir, "MANIFEST.tsv"))

# Optional: the lecturer's reference table for the pHis exercise (notebook 05).
# Looked for in $BIO253_PHIS_REFERENCE, then next to the data folder.
phis_ref <- Sys.getenv("BIO253_PHIS_REFERENCE",
                       file.path(dirname(source_dir), "Analysis_Plan", "pHis",
                                 "pHis_sites_annotated.tsv"))
if (file.exists(phis_ref)) {
  dir.create(file.path(REPO_ROOT, "data", "reference"), showWarnings = FALSE)
  file.copy(phis_ref, file.path(REPO_ROOT, "data", "reference", "pHis_sites_annotated.tsv"),
            overwrite = TRUE)
  message("  copied pHis reference table to data/reference/")
} else {
  message("  (no pHis reference table found -- notebook 05 check will be skipped)")
}

# ---- 3. Write tidy tables ---------------------------------------------------
write_derived <- function(x, name) {
  path <- file.path(derived_dir, name)
  write_tsv(x, path)
  message(sprintf("  wrote %-34s %7d rows", name, nrow(x)))
}
raw <- function(target) file.path(raw_dir, target)

message("Writing tidy tables to data/derived/ (takes ~1 min, reading Excel is slow)")

annotation <- parse_fasta_headers(raw("resources/JE2_proteins.fasta"))
write_derived(annotation, "protein_annotation.tsv.gz")

excel <- c(total = "total/DE_total.xlsx",
           phospho_STYH = "phospho_STYH/DE_phospho_STYH.xlsx",
           phospho_STY = "phospho_STY/DE_phospho_STY.xlsx")
for (ds in names(excel)) {
  write_derived(parse_dea_excel(raw(excel[[ds]])), paste0("dea_", ds, ".tsv.gz"))
  abund <- parse_abundance_excel(raw(excel[[ds]]))
  if (anyNA(abund$sample)) stop("Some sample names could not be parsed in ", excel[[ds]])
  write_derived(abund, paste0("abundances_", ds, ".tsv.gz"))
}

for (s in c("STY", "STYH")) {
  sites <- parse_combined_site(raw(paste0("sites/combined_site_", s, ".tsv")), search = s)
  if (anyNA(sites$sample)) stop("Some sample names could not be parsed in the ", s, " site table")
  write_derived(sites, paste0("sites_", s, ".tsv.gz"))
}

# Sample sheet (27 samples) from the total proteome annotation sheet.
annot_sheet <- readxl::read_excel(raw("total/DE_total.xlsx"), sheet = "annotation")
samples <- annot_sheet |>
  transmute(raw_file_total = raw.file, sample_name_in_excel = sampleName) |>
  bind_cols(select(parse_sample_name(annot_sheet$raw.file), -input, -reinjected)) |>
  arrange(group, replicate) |>
  relocate(sample, genotype, medium, replicate, group)
write_derived(samples, "samples.tsv")

message("Done. Next: open analyses/01_explore_data.qmd")
