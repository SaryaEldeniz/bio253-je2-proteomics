# parse.R -- turn the raw, messy file formats into tidy tables.
#
# You normally do NOT call these yourself: scripts/00_fetch_data.R uses them
# once to write the tidy tables in data/derived/. They are separate functions
# so that we can test them (tests/testthat/).

# ---- Sample names -----------------------------------------------------------

#' Parse genotype / medium / replicate out of any string that contains a
#' sample name, e.g. a raw file name
#'   "20260612_003_C41881_S1143451_JE2anc_TSB_1_input"
#' or a (truncated) sample name "2anc_pulseSN_1".
#'
#' Why: the sample names in the total-proteome tables of the pulse samples
#' were cut off ("2anc_pulseSN_1" instead of "JE2anc_pulseSN_1"), and the raw
#' file names say "pulseSN" while the contrasts say "SNpulse". This function
#' gives every sample ONE canonical name.
parse_sample_name <- function(x) {
  m <- str_match(x, "(anc|stk|stp)_(TSB|SN|pulseSN|SNpulse)_([0-9]+)")
  genotype  <- m[, 2]
  medium    <- if_else(m[, 3] %in% c("pulseSN", "SNpulse"), "SNpulse", m[, 3])
  replicate <- as.integer(m[, 4])
  ok <- !is.na(genotype)
  tibble(
    input     = x,
    # canonical sample name, written like in the lab: JE2anc_pulseSN_1
    sample    = if_else(ok, paste0("JE2", genotype, "_",
                                   if_else(medium == "SNpulse", "pulseSN", medium),
                                   "_", replicate), NA_character_),
    genotype  = factor(genotype, levels = GENOTYPES),
    medium    = factor(medium, levels = MEDIA),
    replicate = replicate,
    group     = factor(if_else(ok, paste(genotype, medium, sep = "_"), NA_character_),
                       levels = GROUPS),
    # stk_pulseSN_1 and _2 of the phospho data were re-injected on another day
    reinjected = str_detect(x, "_rep_")
  )
}

# ---- Protein identifiers and descriptions -----------------------------------

#' The FASTA / DEA protein ids look like "tr|KEY|KEY". We use KEY as the join
#' key everywhere. Ids without "tr|" (e.g. contaminants) are returned as is.
protein_key <- function(protein_Id) {
  if_else(str_detect(protein_Id, "^tr\\|"),
          str_split_fixed(protein_Id, "\\|", 3)[, 2],
          protein_Id)
}

#' Parse the FGCZ FASTA description into useful columns.
#'
#' A typical description is
#'   FGCZ [locus_tag=B7H15_13375] [protein=amino acid ABC transporter
#'   substrate-binding protein] [protein_id=ARG46976.1]
#'   [location=complement(2535968..2536747)] [gbkey=CDS]
#'   BLASTORTHO tr_Q2FVL2_Q2FVL2_STAA8 evalue=0
#'
#' BLASTORTHO is the best BLAST hit in the reference strain NCTC 8325 (STAA8)
#' in UniProt. Use it to look a protein up in UniProt (e.g. Q2FVL2). For
#' reviewed (sp_) entries the entry name also gives a gene name (e.g. SASG).
parse_description <- function(description) {
  d <- coalesce(description, "")
  bracket <- function(field) str_match(d, paste0("\\[", field, "=([^\\]]*)\\]"))[, 2]
  ortho <- str_match(d, "BLASTORTHO (sp|tr)_([A-Z0-9]+)_(\\S+)")
  tibble(
    description   = description,
    # fall back to any B7H15_ number for the 6-frame translated pseudo-genes
    locus_tag     = coalesce(bracket("locus_tag"), str_extract(d, "B7H15_[0-9]+")),
    product       = coalesce(bracket("protein"),
                             if_else(str_detect(d, "frameshifted|pseudo"),
                                     "frameshifted / pseudogene translation", NA_character_)),
    protein_id    = bracket("protein_id"),
    pseudo        = coalesce(bracket("pseudo") == "true", FALSE) | str_detect(d, "frameshifted"),
    location      = bracket("location"),
    strand        = if_else(str_detect(coalesce(bracket("location"), ""), "complement"), "-", "+"),
    staa8_uniprot = ortho[, 3],
    staa8_entry   = ortho[, 4],
    # Swiss-Prot entry names are GENE_STAA8; TrEMBL names are just ACC_STAA8
    staa8_gene    = if_else(ortho[, 2] == "sp", str_remove(ortho[, 4], "_STAA8$"), NA_character_)
  )
}

#' Read all headers of a FASTA file into a protein annotation table.
parse_fasta_headers <- function(fasta_path) {
  headers <- grep("^>", readLines(fasta_path), value = TRUE)
  headers <- str_remove(headers, "^>")
  id   <- str_extract(headers, "^\\S+")
  desc <- str_trim(str_remove(headers, "^\\S+"))
  bind_cols(tibble(protein_Id = id, protein_key = protein_key(id)),
            select(parse_description(desc), -description)) |>
    mutate(is_contaminant = str_detect(protein_Id, "^zz|^CON|contam|Cont_")) |>
    # the first entry of this FASTA is a comment line, not a protein
    filter(!str_detect(protein_key, "^aa_p41881"))
}

# ---- prolfquapp Excel output ------------------------------------------------

#' Read the "diff_exp_analysis" sheet of a prolfquapp DE_*.xlsx into a slim
#' long table: one row per protein (or site) x contrast.
parse_dea_excel <- function(xlsx) {
  x <- suppressWarnings(readxl::read_excel(xlsx, sheet = "diff_exp_analysis"))
  is_site <- "site" %in% names(x)
  out <- x |>
    mutate(protein_key = protein_key(protein_Id)) |>
    select(protein_key, any_of(c("site", "Index", "Peptide", "modAA", "posInProtein",
                                 "SequenceWindow")),
           is_contaminant = CON, contrast, diff, conf.low, conf.high, avgAbd,
           std.error, statistic, df, p.value, FDR, any_of("estimate_type"), modelName)
  # Older prolfquapp runs (the STY phospho DEA from August) used a moderated
  # linear model WITHOUT imputation and have no estimate_type column.
  if (!"estimate_type" %in% names(out)) out$estimate_type <- "observed"
  if (is_site) out <- rename(out, site_id = Index, residue = modAA, position = posInProtein)
  out
}

#' Read the "normalized_abundances" sheet (long format, one row per protein
#' (or site) x sample) and attach the parsed sample information.
parse_abundance_excel <- function(xlsx) {
  x <- suppressWarnings(readxl::read_excel(xlsx, sheet = "normalized_abundances"))
  is_site <- "site" %in% names(x)
  raw_file_col <- intersect(c("raw.file", "file"), names(x))[1]
  raw_col <- intersect(c("exp_medpolish", "Intensity"), names(x))[1]
  samples <- distinct(parse_sample_name(unique(x[[raw_file_col]])))
  out <- x |>
    mutate(protein_key = protein_key(protein_Id),
           input = .data[[raw_file_col]],
           raw_intensity = .data[[raw_col]]) |>
    select(protein_key, any_of(c("site", "Index", "Peptide", "modAA", "posInProtein")),
           is_contaminant = CON, input, raw_intensity,
           any_of("Localization_Probability"), normalized_abundance) |>
    left_join(samples, by = "input") |>
    select(-input)
  if (is_site) out <- rename(out, site_id = Index, residue = modAA, position = posInProtein)
  out
}

# ---- FragPipe combined_site table -------------------------------------------

#' Turn a FragPipe combined_site_*.tsv (one row per site, three columns per
#' sample) into a long table: one row per site x sample.
#'
#' Why: in the wide table, "not detected" is written as intensity 0 and
#' localisation probability 0. A 0 is NOT a measurement, so we set it to NA.
parse_combined_site <- function(tsv, search = NA_character_) {
  x <- read_tsv(tsv, show_col_types = FALSE, guess_max = 10000)
  sample_cols <- str_subset(names(x), " (Localization Probability|Intensity)$")
  # keep only per-sample columns (not "Best Localization Probability", not MaxLFQ)
  sample_cols <- str_subset(sample_cols, "MaxLFQ|^Best ", negate = TRUE)
  long <- x |>
    mutate(is_contaminant = str_detect(Protein, "contam|^zz|CON_")) |>
    select(site_id = Index, protein_key = `Protein ID`, peptide = Peptide,
           best_loc_prob = `Best Localization Probability`, is_contaminant,
           all_of(sample_cols)) |>
    pivot_longer(all_of(sample_cols), names_to = c("input", ".value"),
                 names_pattern = "^(.*) (Localization Probability|Intensity)$") |>
    rename(loc_prob_fragpipe = `Localization Probability`, intensity = Intensity) |>
    # FragPipe sometimes reports a localisation probability for a sample in
    # which the site was NOT quantified (intensity 0). So: decide "detected"
    # from the intensity only, and keep the per-sample LP only where detected.
    mutate(intensity = na_if(intensity, 0),
           loc_prob = if_else(is.na(intensity), NA_real_, loc_prob_fragpipe),
           log2_intensity = log2(intensity))
  site_parts <- str_match(long$site_id, "_([A-Z])([0-9]+)$")
  samples <- distinct(parse_sample_name(unique(long$input)))
  long |>
    mutate(residue = site_parts[, 2],
           position = as.integer(site_parts[, 3]),
           search = search) |>
    left_join(samples, by = "input") |>
    select(search, site_id, protein_key, residue, position, peptide, best_loc_prob,
           sample, genotype, medium, replicate, group, reinjected,
           intensity, log2_intensity, loc_prob, loc_prob_fragpipe, is_contaminant)
}
