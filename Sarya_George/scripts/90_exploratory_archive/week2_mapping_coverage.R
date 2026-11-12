# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

archive <- annotation_zip_path()
meta <- read_tsv(unz(archive, "JE2_protein_metadata.tsv"),
                 col_select = c(protein_key, entry_class, locus_tag, uniprot_acc,
                                ortholog_quality, up_gene_primary, up_length_ratio, in_string),
                 show_col_types = FALSE)
stopifnot(!anyDuplicated(meta$protein_key))

contrasts <- c("TSB_vs_SN_at_anc", "interaction_stk_vs_anc_at_TSB_vs_SN",
               "interaction_stp_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_stk_at_TSB_vs_SN")
joined <- read_dea("total", contrasts) |>
  filter(is.finite(diff), is.finite(FDR)) |>
  select(protein_key, contrast, diff, FDR, estimate_type, staa8_uniprot) |>
  left_join(meta, by = "protein_key", relationship = "many-to-one") |>
  mutate(cds = !is.na(entry_class) & entry_class == "cds",
         good = cds & coalesce(ortholog_quality == "good", FALSE) & !is.na(uniprot_acc),
         hit = FDR < 0.1 & abs(diff) > 1)
print(joined |> count(entry_class, contrast) |> filter(contrast == contrasts[1]), width = Inf)
print(joined |> group_by(contrast) |>
        summarise(tested = n(), joined = sum(!is.na(entry_class)), eligible_cds = sum(cds),
                  good_je2_features = sum(good), good_unique_accessions = n_distinct(uniprot_acc[good]),
                  hits_cds = sum(hit & cds), hits_with_good_mapping = sum(hit & good),
                  good_hit_accessions = n_distinct(uniprot_acc[hit & good]), .groups = "drop"), width = Inf)
print(joined |> filter(contrast == contrasts[1], locus_tag %in% c("B7H15_06350", "B7H15_06355")) |>
        select(locus_tag, uniprot_acc, ortholog_quality, up_gene_primary, up_length_ratio, in_string), width = Inf)

terms <- read_tsv(unz(archive, "JE2_protein_terms.tsv"),
                  col_select = c(protein_key, category), show_col_types = FALSE)
cat("Functional membership rows:", nrow(terms), "\n")
print(terms |> count(category, sort = TRUE), n = 20, width = Inf)
