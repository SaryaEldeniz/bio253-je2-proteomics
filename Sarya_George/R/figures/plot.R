# plot.R -- the plots you will make again and again.

#' Per-sample abundances of one or a few proteins in all 9 groups.
#'
#' Always look at this plot BEFORE you trust a fold change: are the three
#' replicates consistent? Are values missing? Is one sample an outlier?
#'
#' @param query passed to find_protein() (locus tag, product text, ...)
#' @param abundances table from read_abundances(); default: total proteome.
#' @param max_proteins safety limit so a vague query does not draw 300 panels.
#' @return a ggplot. Numbers on top = detected replicates out of 3.
plot_protein <- function(query, abundances = read_abundances("total"), max_proteins = 12) {
  keys <- find_protein(query)$protein_key
  x <- filter(abundances, protein_key %in% keys)
  if (nrow(x) == 0) stop("None of the matching proteins is in this abundance table.")
  if (n_distinct(x$protein_key) > max_proteins) {
    warning("More than ", max_proteins, " proteins match; showing the first ", max_proteins,
            ". Make your query more specific.")
    x <- filter(x, protein_key %in% head(unique(x$protein_key), max_proteins))
  }
  x <- mutate(x, label = protein_label(locus_tag, product)) |> add_x_position()
  # how many of the 3 replicates have a value?
  n_detected <- x |>
    group_by(label, genotype, medium, x_pos) |>
    summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop") |>
    mutate(text = paste0(n, "/3"))
  ggplot(x, aes(x_pos, normalized_abundance, colour = genotype)) +
    geom_boxplot(aes(group = interaction(genotype, medium)), width = 0.22,
                 outlier.shape = NA, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.04, height = 0, seed = 1), na.rm = TRUE) +
    geom_text(data = n_detected, aes(label = text, y = Inf), vjust = 1.3, size = 2.7,
              show.legend = FALSE) +
    facet_wrap(~label, scales = "free_y") +
    medium_x_axis() +
    scale_colour_manual(values = GENOTYPE_COLOURS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = "medium", y = "normalised abundance (log2)",
         caption = "numbers = replicates with a value (of 3); missing = not detected")
}

# Fixed x positions: one slot per medium, genotypes side by side. (With
# ggplot's automatic dodging, boxes jump around when a group is completely
# missing -- exactly the case we care about for knock-out mutants.)
add_x_position <- function(x) {
  offset <- c(anc = -0.27, stk = 0, stp = 0.27)
  mutate(x, x_pos = as.numeric(factor(medium, levels = MEDIA)) +
           unname(offset[as.character(genotype)]))
}
medium_x_axis <- function() {
  scale_x_continuous(breaks = seq_along(MEDIA), labels = MEDIA,
                     limits = c(0.6, length(MEDIA) + 0.4))
}

#' Per-sample log2 intensities of phospho sites (from read_sites()).
plot_sites <- function(site_ids, sites) {
  x <- filter(sites, site_id %in% site_ids) |>
    mutate(label = paste0(coalesce(locus_tag, protein_key), " ", residue, position, "\n",
                          str_trunc(coalesce(product, ""), 35))) |>
    add_x_position()
  if (nrow(x) == 0) stop("None of these site_ids is in the table.")
  ggplot(x, aes(x_pos, log2_intensity, colour = genotype)) +
    geom_point(position = position_jitter(width = 0.04, height = 0, seed = 1), na.rm = TRUE) +
    facet_wrap(~label, scales = "free_y") +
    medium_x_axis() +
    scale_colour_manual(values = GENOTYPE_COLOURS) +
    labs(x = "medium", y = "site intensity (log2, not normalised)",
         caption = "missing points = site not detected in that sample")
}

#' Volcano plot for one contrast of a read_dea() table.
#'
#' @param highlight optional protein_keys or locus tags to label in the plot.
volcano <- function(dea, contrast_name, fdr_cutoff = FDR_CUTOFF, lfc_cutoff = LFC_CUTOFF,
                    highlight = NULL) {
  x <- filter(dea, contrast == contrast_name) |>
    mutate(significant = FDR < fdr_cutoff & abs(diff) > lfc_cutoff,
           class = case_when(significant & diff > 0 ~ "up (A higher)",
                             significant & diff < 0 ~ "down (B higher)",
                             TRUE ~ "not significant"))
  if (nrow(x) == 0) stop("Contrast '", contrast_name, "' not found.")
  p <- ggplot(x, aes(diff, -log10(FDR), colour = class, shape = estimate_type)) +
    geom_point(alpha = 0.6, size = 1.4) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("up (A higher)" = "#D55E00", "down (B higher)" = "#0072B2",
                                   "not significant" = "grey70")) +
    scale_shape_manual(values = c(observed = 16, lod_imputed = 4)) +
    labs(title = contrast_name, subtitle = explain_contrast(contrast_name),
         x = "diff = log2 fold change (A - B)", y = "-log10(FDR)") +
    theme(plot.subtitle = element_text(size = 7))
  if (!is.null(highlight)) {
    lab <- filter(x, protein_key %in% highlight | locus_tag %in% highlight)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = locus_tag),
                                        colour = "black", size = 3, show.legend = FALSE)
    } else {
      p <- p + geom_text(data = lab, aes(label = locus_tag), colour = "black", size = 3,
                         vjust = -0.6, show.legend = FALSE)
    }
  }
  p
}

#' Count significant proteins/sites per contrast.
count_hits <- function(dea, fdr_cutoff = FDR_CUTOFF, lfc_cutoff = LFC_CUTOFF) {
  dea |>
    group_by(contrast) |>
    summarise(n_tested = n(),
              n_up = sum(FDR < fdr_cutoff & diff > lfc_cutoff, na.rm = TRUE),
              n_down = sum(FDR < fdr_cutoff & diff < -lfc_cutoff, na.rm = TRUE),
              n_hits = n_up + n_down,
              .groups = "drop")
}
