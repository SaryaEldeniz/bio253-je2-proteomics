# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

week2 <- Sys.getenv("JE2_WEEK2_DIR", unset=file.path(dirname(dirname(DATA_ROOT)),"Week 2"))
for (basename in c("SAstrainSpecificIDs_to_Uniprot.xlsx",
                   "p3404_JE2_prolfqua_2023_db5wPseudo_allResultsWithMetaNWideMat_nocontNnorevs.xlsx")) {
  file <- file.path(week2, basename)
  sheets <- readxl::excel_sheets(file)
  cat("\nWORKBOOK", basename, "\n")
  print(sheets)
  for (sheet in sheets) {
    rows <- readxl::read_excel(file, sheet = sheet, n_max = 2, .name_repair = "minimal")
    cat("\nSHEET", sheet, "\nHEADERS", paste(names(rows), collapse = " | "), "\n")
    print(as.data.frame(rows[, seq_len(min(ncol(rows), 8))]), row.names = FALSE)
  }
}

contrasts <- c("TSB_vs_SN_at_anc", "interaction_stk_vs_anc_at_TSB_vs_SN",
               "interaction_stp_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_stk_at_TSB_vs_SN")
current <- read_dea("total", contrasts) |>
  filter(is.finite(diff), is.finite(FDR)) |>
  mutate(has_accession = !is.na(staa8_uniprot) & nzchar(staa8_uniprot),
         hit = FDR < 0.1 & abs(diff) > 1)
cat("\nCURRENT TOTAL PROTEOME: NCTC8325 ANNOTATION COVERAGE\n")
print(current |> group_by(contrast) |>
        summarise(features = n(), mapped_features = sum(has_accession),
                  distinct_accessions = n_distinct(staa8_uniprot[has_accession]),
                  hit_features = sum(hit), mapped_hits = sum(hit & has_accession),
                  distinct_hit_accessions = n_distinct(staa8_uniprot[hit & has_accession]),
                  .groups = "drop"), width = Inf)
print(current |> filter(contrast == "TSB_vs_SN_at_anc", locus_tag %in% c("B7H15_06350", "B7H15_06355", "B7H15_05965")) |>
        select(locus_tag, protein_id, staa8_uniprot, product), width = Inf)
