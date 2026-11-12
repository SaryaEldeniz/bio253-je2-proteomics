# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

# A small, read-only example: one contrast, KEGG pathway memberships from STRING.
out_dir <- "results/exploratory"
x <- read_dea("total", "TSB_vs_SN_at_anc") |>
  filter(is.finite(diff), is.finite(FDR))
aliases <- read_tsv(file.path(DATA_ROOT, "raw/resources/string_aliases.txt"), show_col_types = FALSE)
names(aliases)[1] <- "string_id"
mapping <- aliases |>
  filter(alias %in% x$locus_tag) |>
  distinct(locus_tag = alias, string_id)
stopifnot(!anyDuplicated(mapping$locus_tag))
mapped <- x |> inner_join(mapping, by = "locus_tag")
terms <- read_tsv(file.path(DATA_ROOT, "raw/resources/string_enrichment_terms.txt"), show_col_types = FALSE)
names(terms)[1] <- "string_id"
terms <- terms |>
  filter(category == "KEGG (Kyoto Encyclopedia of Genes and Genomes)") |>
  distinct(string_id, category, term, description)

run_ora <- function(tab, label) {
  universe <- unique(tab$string_id)
  pathways <- terms |>
    filter(string_id %in% universe) |>
    group_by(category, term, description) |>
    summarise(members = list(unique(string_id)), .groups = "drop") |>
    mutate(background_in_term = lengths(members)) |>
    filter(background_in_term >= 5, background_in_term <= 500)
  directions <- list(
    "Higher in SN" = tab$string_id[tab$FDR < 0.1 & tab$diff < -1],
    "Higher in TSB" = tab$string_id[tab$FDR < 0.1 & tab$diff > 1]
  )
  result <- bind_rows(lapply(names(directions), function(direction) {
    foreground <- unique(directions[[direction]])
    n <- length(foreground)
    N <- length(universe)
    bind_rows(lapply(seq_len(nrow(pathways)), function(i) {
      members <- pathways$members[[i]]
      hit_ids <- intersect(foreground, members)
      k <- length(hit_ids)
      K <- length(members)
      loci <- tab$locus_tag[match(hit_ids, tab$string_id)]
      tibble(analysis = label, direction = direction,
             category = pathways$category[i], term = pathways$term[i],
             description = pathways$description[i],
             foreground_in_term = k, foreground_total = n,
             background_in_term = K, background_total = N,
             expected_hits = n * K / N,
             enrichment_ratio = if (n > 0) (k / n) / (K / N) else NA_real_,
             p_value = if (n > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else 1,
             hit_locus_tags = paste(loci, collapse = ";"))
    }))
  }))
  # Adjust jointly over both directions and every eligible KEGG pathway.
  result |> mutate(ORA_FDR = p.adjust(p_value, method = "BH")) |>
    arrange(ORA_FDR, p_value)
}

primary <- run_ora(mapped, "all estimates")
sensitivity <- run_ora(filter(mapped, estimate_type == "observed"), "observed only")
write_tsv(bind_rows(primary, sensitivity), file.path(out_dir, "I_preview.tsv"))

coverage <- x |>
  mutate(mapped_to_string = locus_tag %in% mapping$locus_tag,
         higher_SN = FDR < 0.1 & diff < -1,
         higher_TSB = FDR < 0.1 & diff > 1) |>
  group_by(estimate_type) |>
  summarise(tested_proteins = n(), mapped_to_string = sum(mapped_to_string),
            higher_in_SN = sum(higher_SN), higher_in_TSB = sum(higher_TSB),
            .groups = "drop")
write_tsv(coverage, file.path(out_dir, "I_mapping_summary.tsv"))
print(coverage)
cat("Eligible KEGG pathways", nrow(primary)/2, "\n")
print(primary |> group_by(direction) |> slice_head(n = 5) |>
        select(direction, term, description, foreground_in_term, foreground_total,
               background_in_term, background_total, enrichment_ratio, ORA_FDR),
      n = 10, width = Inf)
cat("OBSERVED ONLY\n")
print(sensitivity |> group_by(direction) |> slice_head(n = 5) |>
        select(direction, term, description, foreground_in_term, foreground_total,
               background_in_term, background_total, enrichment_ratio, ORA_FDR),
      n = 10, width = Inf)
