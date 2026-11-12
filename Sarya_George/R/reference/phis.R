# phis.R -- small building blocks for judging phospho-histidine (pHis) sites
# (notebook 05). They only compute EVIDENCE; deciding what counts as
# "convincing" is your job.

#' A FragPipe "Peptide" entry can list several sequences separated by ";"
#' (e.g. a missed-cleavage form and a short form). FragPipe writes possibly
#' phosphorylated residues in lower case, so we upper-case before comparing.
peptide_variants <- function(peptide) {
  str_split(toupper(peptide), ";")
}

#' Number of alternative phospho-acceptors (S, T or Y) in a peptide entry.
#' Conservative: the FEWEST over all listed variants (if one variant has no
#' S/T/Y at all, the His is the only possible acceptor in that variant).
n_sty_acceptors <- function(peptide) {
  vapply(peptide_variants(peptide), \(v) min(str_count(v, "[STY]")), numeric(1))
}

#' Is any S/T/Y on the peptide ALSO marked (lower case) as phosphorylated?
co_phospho_sty <- function(peptide) {
  str_detect(peptide, "[sty]")
}

#' Is the peptide "explained" by the STY search? TRUE only if EVERY variant of
#' the entry also appears as a phosphopeptide in the STY-only search.
all_variants_in <- function(peptide, peptide_set) {
  vapply(peptide_variants(peptide), \(v) all(v %in% peptide_set), logical(1))
}

#' TRUE if at least one variant appears in the STY-only search.
any_variant_in <- function(peptide, peptide_set) {
  vapply(peptide_variants(peptide), \(v) any(v %in% peptide_set), logical(1))
}

#' All (upper-cased) phosphopeptide sequences of a read_sites() table.
phosphopeptide_set <- function(sites) {
  unique(unlist(peptide_variants(unique(sites$peptide))))
}

#' Localisation-probability classes used in the course.
lp_class <- function(lp) {
  cut(lp, c(-Inf, 0.5, 0.75, 0.9, 0.99, Inf),
      labels = c("<0.5", "0.5-0.75", "0.75-0.9", "0.9-0.99", ">=0.99"), right = FALSE)
}

#' Compare YOUR tier assignment with the lecturer's reference table
#' (data/reference/pHis_sites_annotated.tsv, if it was provided).
#'
#' @param my_tiers a table with columns site_id and tier, where tier starts
#'   with the tier number 1-6 (e.g. "1 strongest ...").
#' @return a cross-table: rows = your tier, columns = reference tier. Only the
#'   counts are shown, so it does not give away which site is which.
check_his_tiers <- function(my_tiers) {
  ref_path <- file.path(REPO_ROOT, "data", "reference", "pHis_sites_annotated.tsv")
  if (!file.exists(ref_path)) {
    message("No reference table at ", ref_path, " -- ask the lecturer.")
    return(invisible(NULL))
  }
  ref <- read_tsv(ref_path, show_col_types = FALSE) |>
    transmute(site_id = Index, reference_tier = str_extract(tier, "^[0-9]"))
  mine <- transmute(my_tiers, site_id, my_tier = str_extract(tier, "^[0-9]"))
  joined <- full_join(mine, ref, by = "site_id")
  agree <- sum(joined$my_tier == joined$reference_tier, na.rm = TRUE)
  message(agree, " of ", nrow(ref), " His sites agree with the reference.")
  table(my_tier = joined$my_tier, reference_tier = joined$reference_tier, useNA = "ifany")
}
