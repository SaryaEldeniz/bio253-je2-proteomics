summarize_step07_parent_comparison <- function(phospho, total, contrasts,
                                               fdr_cutoff = 0.10, log2fc_cutoff = 1) {
  summarize_contrast <- function(contrast_name) {
    sites <- subset(phospho, contrast == contrast_name & !is_contaminant)
    proteins <- subset(total, contrast == contrast_name & !is_contaminant)
    stopifnot(nrow(sites) > 0, nrow(proteins) > 0, !anyDuplicated(proteins$protein_key))

    # Match each phosphosite to the corresponding protein in this contrast.
    parent_row <- match(sites$protein_key, proteins$protein_key)
    parent_diff <- proteins$diff[parent_row]
    parent_fdr <- proteins$FDR[parent_row]
    parent_type <- proteins$estimate_type[parent_row]

    site_available <- is.finite(sites$diff) & is.finite(sites$FDR)
    parent_available <- is.finite(parent_diff) & is.finite(parent_fdr)
    site_significant <- site_available & sites$FDR < fdr_cutoff & abs(sites$diff) > log2fc_cutoff
    parent_significant <- parent_available & parent_fdr < fdr_cutoff & abs(parent_diff) > log2fc_cutoff
    same_sign <- sign(sites$diff) == sign(parent_diff)

    data.frame(
      contrast = contrast_name,
      sites_tested = sum(site_available),
      sites_with_parent = sum(site_available & parent_available),
      significant_sites = sum(site_significant),
      phospho_lod = sum(site_significant & sites$estimate_type == "lod_imputed"),
      parent_lod = sum(site_significant & parent_available & parent_type == "lod_imputed"),
      same_direction = sum(site_significant & parent_significant & same_sign),
      opposite_direction = sum(site_significant & parent_significant & !same_sign),
      parent_below_thresholds = sum(site_significant & parent_available & !parent_significant),
      parent_unavailable = sum(site_significant & !parent_available)
    )
  }

  do.call(rbind, lapply(contrasts, summarize_contrast))
}
