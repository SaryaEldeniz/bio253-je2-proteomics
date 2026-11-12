# find.R -- look proteins up and explain contrasts.

#' Find proteins by locus tag, protein_id, product text or STAA8 ortholog.
#'
#' @param query text or regular expression, case-insensitive. Several queries
#'   can be given as a vector, e.g. c("B7H15_13375", "B7H15_13370").
#' @return annotation rows of all matching proteins.
#' @examples
#'   find_protein("B7H15_13375")          # locus tag
#'   find_protein("ARG46976")             # GenBank protein id
#'   find_protein("histidine kinase")     # product text
#'   find_protein("Q2FVL2")               # UniProt id of the NCTC8325 ortholog
find_protein <- function(query, annotation = read_annotation()) {
  pattern <- regex(paste(query, collapse = "|"), ignore_case = TRUE)
  searchable <- with(annotation, paste(protein_key, locus_tag, protein_id, product,
                                       staa8_uniprot, staa8_entry, sep = " | "))
  hits <- annotation[str_detect(searchable, pattern), ]
  if (nrow(hits) == 0) message("No protein matches '", paste(query, collapse = "|"), "'.")
  hits
}

#' A short, human readable label for plots: "B7H15_13375 amino acid ABC tra..."
protein_label <- function(locus_tag, product, width = 40) {
  paste0(coalesce(locus_tag, "?"), "\n", str_trunc(coalesce(product, ""), width))
}

#' Say in words what a contrast name means.
#'
#' prolfquapp contrast names follow "A_vs_B" (optionally "_at_<level>"):
#'   diff = mean log2 abundance of A minus mean log2 abundance of B.
explain_contrast <- function(contrast) {
  if (str_starts(contrast, "interaction_")) {
    m <- str_match(contrast, "^interaction_(.+)_vs_(.+)_at_(.+)_vs_(.+)$")
    return(sprintf(paste0("%s: difference of differences = (%s - %s in %s) - (%s - %s in %s). ",
                          "Positive: the %s-vs-%s difference is larger in %s than in %s."),
                   contrast, m[2], m[3], m[4], m[2], m[3], m[5], m[2], m[3], m[4], m[5]))
  }
  m <- str_match(contrast, "^(.+)_vs_(.+?)(?:_at_(.+))?$")
  where <- if (is.na(m[4])) "averaged over the other factor" else paste("within", m[4])
  sprintf(paste0("%s: diff = log2(%s) - log2(%s), %s. ",
                 "diff > 0 means HIGHER in %s; diff < 0 means HIGHER in %s."),
          contrast, m[2], m[3], where, m[2], m[3])
}
