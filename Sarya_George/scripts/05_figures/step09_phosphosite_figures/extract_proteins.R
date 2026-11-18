# Select the hits shown in the Step 9 residue-specific Venn diagrams.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript extract_proteins.R DATA_DIRECTORY OUTPUT_TSV", call. = FALSE)
}
data_dir <- args[1]
media <- c("TSB", "SN")
residues <- c("S", "T", "Y", "H")
mutants <- c("stk", "stp")
contrasts <- as.vector(outer(mutants, media, function(g, m) paste0(g, "_vs_anc_at_", m)))

# Follow the same source-estimate criteria as 09_render_residue_venn.R.
phospho <- read.delim(file.path(data_dir, "derived", "dea_phospho_STYH.tsv.gz"))
tested <- subset(phospho, contrast %in% contrasts & !is_contaminant &
                   is.finite(FDR) & is.finite(diff))
stopifnot(setequal(tested$contrast, contrasts), !anyNA(tested$site_id),
          !anyDuplicated(tested[c("contrast", "site_id")]), all(tested$residue %in% residues))
hits <- subset(tested, estimate_type == "observed" & FDR < 0.10 & abs(diff) > 1)

exclusive_hits <- function(medium, site_residue, mutant) {
  contrast_name <- paste0(mutant, "_vs_anc_at_", medium)
  other_contrast <- paste0(setdiff(mutants, mutant), "_vs_anc_at_", medium)
  group_hits <- subset(hits, contrast == contrast_name & residue == site_residue)
  other_ids <- hits$site_id[hits$contrast == other_contrast & hits$residue == site_residue]
  group_hits <- group_hits[!group_hits$site_id %in% other_ids, ]
  group_hits$medium <- rep(medium, nrow(group_hits))
  group_hits$mutant <- rep(mutant, nrow(group_hits))
  group_hits
}

selected <- do.call(rbind, lapply(media, function(medium) {
  do.call(rbind, lapply(residues, function(residue) {
    do.call(rbind, lapply(mutants, function(mutant) exclusive_hits(medium, residue, mutant)))
  }))
}))
stopifnot(nrow(selected) > 0, !anyNA(selected$position))

annotation <- read.delim(file.path(data_dir, "derived", "protein_annotation.tsv.gz"))
stopifnot(!anyDuplicated(annotation$protein_key))
selected_annotation <- annotation[match(selected$protein_key, annotation$protein_key), ]
stopifnot(!anyNA(selected_annotation[c("staa8_uniprot", "product", "locus_tag")]))

fasta <- readLines(file.path(data_dir, "raw", "resources", "JE2_proteins.fasta"), warn = FALSE)
headers <- grep("^>", fasta)
sequence_ends <- c(headers[-1] - 1L, length(fasta))
fasta_keys <- sub("^>tr\\|([^|]+)\\|.*", "\\1", fasta[headers])
protein_keys <- unique(selected$protein_key)
fasta_records <- match(protein_keys, fasta_keys)
stopifnot(!anyNA(fasta_records), all(vapply(
  protein_keys, function(key) sum(fasta_keys == key) == 1L, logical(1)
)))

sequences <- vapply(fasta_records, function(i) {
  paste0(fasta[seq.int(headers[i] + 1L, sequence_ends[i])], collapse = "")
}, character(1))
sequence_by_hit <- sequences[match(selected$protein_key, protein_keys)]
stopifnot(all(substr(sequence_by_hit, selected$position, selected$position) == selected$residue))

# This table is an intermediate passed to the drawing code through a temporary directory.
result <- data.frame(
  accession = selected_annotation$staa8_uniprot,
  name = selected_annotation$product,
  locus_tag = selected_annotation$locus_tag,
  sequence = sequence_by_hit,
  position = selected$position,
  residue = selected$residue,
  medium = selected$medium,
  mutant = selected$mutant
)
write.table(result, args[2], sep = "\t", quote = TRUE, qmethod = "double",
            row.names = FALSE, fileEncoding = "UTF-8")
