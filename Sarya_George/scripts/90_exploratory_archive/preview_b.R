# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

out_dir <- "figures/exploratory"
work_dir <- "results/exploratory"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
targets <- tibble::tribble(
  ~locus_tag, ~short_name,
  "B7H15_13375", "TcyA",
  "B7H15_13370", "TcyB",
  "B7H15_13365", "TcyC",
  "B7H15_02105", "TcyP-like (alternative importer)"
)
ab <- read_abundances("total") |>
  inner_join(targets, by = "locus_tag")
cat("Abundance columns:", paste(names(ab), collapse = " | "), "\n")
cat("Records:", nrow(ab), "Detected:", sum(!is.na(ab$normalized_abundance)), "\n")

# Each facet shows measured sample values; missing observations are never zero-filled.
x <- add_x_position(ab) |>
  mutate(short_name = factor(short_name, levels = targets$short_name))
n_detect <- x |>
  group_by(short_name, genotype, medium, x_pos) |>
  summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop")
means <- x |>
  group_by(short_name, genotype, medium, x_pos) |>
  summarise(avg = if (all(is.na(normalized_abundance))) NA_real_ else
              mean(normalized_abundance, na.rm = TRUE), .groups = "drop")
p <- ggplot(x, aes(x_pos, normalized_abundance, colour = genotype)) +
  geom_point(position = position_jitter(width = 0.045, height = 0, seed = 22),
             size = 2.4, alpha = 0.9, na.rm = TRUE) +
  geom_point(data = means, aes(y = avg), shape = 95, size = 8, na.rm = TRUE) +
  geom_text(data = n_detect, aes(y = Inf, label = paste0(n, "/3")),
            vjust = 1.6, size = 3, show.legend = FALSE) +
  facet_wrap(~short_name, ncol = 2, scales = "free_y") +
  medium_x_axis() +
  scale_colour_manual(values = GENOTYPE_COLOURS,
                      labels = c(anc = "Ancestor", stk = "Kinase mutant", stp = "Phosphatase mutant")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(title = "B | Cystine importers: actual sample measurements",
       subtitle = "TSB = control; SNpulse = 15 min; SN = about 15 h of supernatant exposure",
       x = NULL, y = "Normalized protein abundance (log2 scale)", colour = NULL,
       caption = "Dots = measured samples; bars = group means; labels = detected replicates out of 3.\nEach protein has its own y-axis scale. Missing measurements are omitted, not replaced by zero.") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 17),
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 9),
        plot.margin = margin(14, 18, 12, 14))
ggsave(file.path(out_dir, "B_transporter_data.png"), p, width = 11, height = 7.6, dpi = 160)

dea <- read_dea("total", c("TSB_vs_SN_at_anc", "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp")) |>
  inner_join(targets, by = "locus_tag") |>
  select(short_name, locus_tag, contrast, diff, FDR, estimate_type,
         any_of(c("conf.low", "conf.high"))) |>
  arrange(short_name, contrast)
write_tsv(dea, file.path(work_dir, "B_preview_dea.tsv"))
print(dea, n = 20, width = 160)
write_tsv(n_detect, file.path(work_dir, "B_detection_counts.tsv"))
cat("Total proteins tested per ancestor SN contrast:",
    nrow(read_dea("total", "TSB_vs_SN_at_anc")), "\n")
