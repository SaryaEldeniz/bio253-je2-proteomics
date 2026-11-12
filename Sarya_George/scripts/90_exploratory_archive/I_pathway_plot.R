# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

scratch <- "results/exploratory"
output <- "figures/exploratory"
dir.create(output, showWarnings = FALSE)
x <- read_tsv(file.path(scratch, "I_preview.tsv"), show_col_types = FALSE) |>
  filter(analysis == "all estimates", term %in% c("map00020", "map01200", "map00650"),
         (term == "map00650" & direction == "Higher in SN") |
           (term != "map00650" & direction == "Higher in TSB")) |>
  mutate(label = case_when(
    term == "map00020" ~ "Citrate cycle (TCA cycle)\nHigher in TSB / lower in SN",
    term == "map01200" ~ "Carbon metabolism\nHigher in TSB / lower in SN",
    TRUE ~ "Butanoate metabolism\nHigher in SN"),
    label = factor(label, levels = c("Butanoate metabolism\nHigher in SN",
                                    "Carbon metabolism\nHigher in TSB / lower in SN",
                                    "Citrate cycle (TCA cycle)\nHigher in TSB / lower in SN")),
    annotation = sprintf("%.2fx   |   ORA FDR %.4f", enrichment_ratio, ORA_FDR))

p <- ggplot(x, aes(enrichment_ratio, label, fill = direction)) +
  geom_col(width = 0.52, show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#74808c", linewidth = 0.45) +
  geom_text(aes(label = annotation), hjust = 0, nudge_x = 0.10,
            size = 4.0, colour = "#172936") +
  scale_fill_manual(values = c("Higher in TSB" = "#2e6f9e", "Higher in SN" = "#ba653e")) +
  scale_x_continuous(limits = c(0, 7), breaks = 0:6, expand = expansion(mult = c(0, 0))) +
  labs(title = "A first pathway enrichment preview",
       subtitle = "JE2 ancestor | TSB versus ~15 h supernatant (SN) | total proteome",
       x = "Overrepresentation relative to the measured background (fold enrichment)", y = NULL,
       caption = paste(
         "Protein selection: FDR < 0.10 and |log2(TSB / SN)| > 1; directions tested separately.",
         "ORA: 2,231 unique measured STRING proteins; 91 KEGG terms x 2 directions; joint BH adjustment over 182 tests.",
         "Shown: three terms passing exploratory ORA FDR < 0.10. Butanoate does not pass 0.05 in this primary analysis.",
         "Primary analysis includes imputed estimates. TCA and carbon metabolism persist with observed-only estimates.",
         "Pathway annotation enrichment does not establish metabolic flux. Dashed line: no enrichment (1x).",
         sep = "\n")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(size = 22, face = "bold", colour = "#172936", margin = margin(b = 8)),
        plot.subtitle = element_text(size = 13, colour = "#46596a", margin = margin(b = 20)),
        axis.text.y = element_text(size = 12, colour = "#172936", margin = margin(r = 12)),
        axis.text.x = element_text(colour = "#46596a"),
        axis.title.x = element_text(size = 12, margin = margin(t = 13)),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e2e7eb"),
        plot.caption = element_text(hjust = 0, size = 9.3, colour = "#46596a", lineheight = 1.3,
                                    margin = margin(t = 22)),
        plot.margin = margin(25, 28, 20, 25))
ggsave(file.path(output, "I_pathway_preview.png"), p, width = 13.2, height = 7.7,
       dpi = 160, bg = "white")
