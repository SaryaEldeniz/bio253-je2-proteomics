# 2026-11-21: separate titles and shared legend

Request: one boxed SVG legend for the six interaction volcanoes, plus a separate SVG title for each triplet.

Added `scripts/poster_a0/export_volcano_labels.R`, using R/grid and Cairo vector output. Exported the matching existing total-proteome and phosphosite titles as tight transparent SVGs at bold Arial 32 pt. Exported a compact white boxed common legend with 24-pt Arial text, enlarged matching symbols and the shared FDR/effect selection rule. Physical sizes are measured using the actual font and declared in mm. PNG previews and a three-SVG ZIP accompany the assets in `figures/volcano_labels/`.

Visual inspection confirmed all symbols and text fit, with minimal margins. SVG export asserts absence of embedded raster images. Source plot files, estimates, FDRs and the phosphosite display-omission policy were not changed. When assembling the poster, retain the existing left phosphosite panel's Pyk S538 omission notice. Codex prepared the reproducible R exporter and checked the three previews.

Follow-up: added two square alternatives with the same symbols and 24-pt type: 100 × 100 mm with a selection-criteria block, and 80 × 80 mm with only the five colour/model labels. Both PNG previews were visually checked for fit. The original rectangular legend and titles remain available; the ZIP and size manifest now include all five SVGs.

Title follow-up: updated both transparent title SVGs to include the genotype-by-medium heading and the delta-delta log2 abundance formula. A/B are written as Genotype 1/Genotype 2, with a key identifying them as the first/second genotype named in each panel. Formula orientation remains TSB/PASN for genotype 1 minus TSB/PASN for genotype 2, equivalent to the native source contrasts. The taller 60-mm exports separate the fraction parentheses from the title and key. Heading size is 32 pt; formula/key are 24 pt (standard mathematical subscript scaling). Updated PNG previews were checked for clipping and spacing; the five-asset ZIP was regenerated. No volcano coordinates or legends changed.

Final requested revision: replaced the combined title/formula files with heading-only SVGs/PNGs at the same filenames. Exported one separate shared `interaction_formula_G.svg`, using G1/G2 and medium subscripts, without the explanatory line. Updated the exporter, size manifest, README and six-asset ZIP; inspected all three previews. The superseded combined layouts are overwritten, avoiding duplicate title variants. No plot or legend changes.

Formula typography follow-up: made the standalone formula bold, including abundance, G labels and medium subscripts, at the same nominal 24-pt size. Re-measured its bounds, regenerated the SVG/PNG and ZIP, and checked the preview for fit. Headings and legends retain their previous appearance.
