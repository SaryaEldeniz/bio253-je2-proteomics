# New Project I helpers: planned contrasts and the design figure.

project_i_contrasts <- function() {
  tibble::tribble(
    ~contrast, ~priority, ~question,
    "TSB_vs_SN_at_anc", "core", "Long PASN response in the ancestor",
    "TSB_vs_SN_at_stk", "core", "Long PASN response in the kinase mutant",
    "TSB_vs_SN_at_stp", "core", "Long PASN response in the phosphatase mutant",
    "stk_vs_anc_at_TSB", "core", "Kinase mutant effect in control medium; first Venn",
    "stp_vs_anc_at_TSB", "core", "Phosphatase mutant effect in control medium; first Venn",
    "stk_vs_anc_at_SN", "core", "Kinase mutant effect during long PASN exposure; second Venn",
    "stp_vs_anc_at_SN", "core", "Phosphatase mutant effect during long PASN exposure; second Venn",
    "interaction_stk_vs_anc_at_TSB_vs_SN", "core", "Does the kinase mutant differ in its PASN response?",
    "interaction_stp_vs_anc_at_TSB_vs_SN", "core", "Does the phosphatase mutant differ in its PASN response?",
    "interaction_stp_vs_stk_at_TSB_vs_SN", "follow_up", "Do the two mutants differ in their PASN response?",
    "TSB_vs_SNpulse_at_anc", "optional", "Early 15-minute PASN response in the ancestor",
    "SNpulse_vs_SN_at_anc", "optional", "Early versus long PASN exposure in the ancestor"
  )
}

# Coefficients refer to the nine genotype-medium cell means, not paired samples.
# Pooled contrasts are intentionally excluded: their model weighting needs review.
contrast_weights <- function(contrast) {
  weights <- setNames(rep(0, length(GROUPS)), GROUPS)
  interaction <- str_match(contrast, "^interaction_(anc|stk|stp)_vs_(anc|stk|stp)_at_(TSB|SNpulse|SN)_vs_(TSB|SNpulse|SN)$")
  environmental <- str_match(contrast, "^(TSB|SNpulse|SN)_vs_(TSB|SNpulse|SN)_at_(anc|stk|stp)$")
  mutant <- str_match(contrast, "^(anc|stk|stp)_vs_(anc|stk|stp)_at_(TSB|SNpulse|SN)$")
  if (!is.na(interaction[1, 1])) {
    a <- interaction[1, 2]; b <- interaction[1, 3]
    x <- interaction[1, 4]; y <- interaction[1, 5]
    weights[paste(a, x, sep = "_")] <- 1
    weights[paste(b, x, sep = "_")] <- -1
    weights[paste(a, y, sep = "_")] <- -1
    weights[paste(b, y, sep = "_")] <- 1
  } else if (!is.na(environmental[1, 1])) {
    weights[paste(environmental[1, 4], environmental[1, 2], sep = "_")] <- 1
    weights[paste(environmental[1, 4], environmental[1, 3], sep = "_")] <- -1
  } else if (!is.na(mutant[1, 1])) {
    weights[paste(mutant[1, 2], mutant[1, 4], sep = "_")] <- 1
    weights[paste(mutant[1, 3], mutant[1, 4], sep = "_")] <- -1
  } else stop("Unsupported condition-specific contrast: ", contrast)
  weights
}

contrast_formula <- function(contrast) {
  weights <- contrast_weights(contrast)
  positive <- names(weights)[weights == 1]
  negative <- names(weights)[weights == -1]
  paste0("(", paste(positive, collapse = " + "), ") - (",
         paste(negative, collapse = " + "), ")")
}
