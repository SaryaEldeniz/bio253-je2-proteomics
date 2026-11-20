samples <- read.delim("data/derived/samples.tsv")

samples$genotype <- factor(samples$genotype, levels = c("anc", "stk", "stp"))
samples$medium <- factor(samples$medium, levels = c("TSB", "SNpulse", "SN"))

design <- with(samples, table(genotype, medium))
print(design)
cat("Total samples:", nrow(samples), "\n")

# Check that each group contains one sample per replicate.
stopifnot(
  nrow(samples) == 27L,
  !anyNA(samples$sample),
  all(nzchar(samples$sample)),
  anyDuplicated(samples$sample) == 0L,
  setequal(samples$replicate, 1:3),
  all(design == 3L),
  all(with(samples, table(genotype, medium, replicate)) == 1L)
)

cat("Design check passed.\n")
