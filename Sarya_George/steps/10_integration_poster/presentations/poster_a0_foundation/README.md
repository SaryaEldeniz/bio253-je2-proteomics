# Editable A0 poster foundation

Open **JE2_A0_poster_foundation.pptx**. Its custom slide size is exactly **841 × 1189 mm, portrait**. The PDF and PNG are viewing copies of the same draft. This is a foundation for figure placement, not a completed results poster.

## Filled content and intentional placeholders

The header uses the supplied ETH–UZH–FGCZ image, all six author names and four affiliations, and `* Correspondence: george.kahvedjian@uzh.ch`. The proposed neutral title adds **phosphosite** to cover both measurement layers; species names are italicized.

Background and hypothesis are filled. Methods contains only the Project analysis SVG, plus an editable text slot beneath it. The sample-design matrix and complementary-proteomics-layers figure are excluded as requested. The smaller workflow labels print at approximately 25.7 pt at their placed size. The original workflow SVG is preserved; the compact poster variant is `../../figures/methods_svg/03_analysis_workflow_poster.svg`.

Results contain named, editable figure frames. Discussion/conclusion, follow-up experiments, acknowledgements and the two QR areas remain visibly reserved. No numerical results, conclusions, funding acknowledgements or scannable QR destinations have been fabricated.

## Typography and formatting

| Element | PowerPoint size |
| --- | ---: |
| Title | 72 pt, bold; 65-mm title block |
| Authors | 36 pt, bold |
| Affiliations and contact | 24 pt |
| Main section headings | 42 pt |
| Panel headings and hypothesis heading | 36 pt |
| Background and hypothesis | 28 pt |
| Figure-slot labels and future caption minimum | 24 pt |

The lecturer specifies A0 portrait, Arial or Times New Roman, title 72–120 pt, headings at least 36 pt, body at least 24 pt, left-aligned text and legibility from about 1.5 m. This draft uses Arial. Reference: `Lecture_Academic_posters.pdf`, pages 61–63, 68 and 70. Margins of 20 mm, 14-mm column gutters and rules/borders of 1.5–2 pt are consistent design choices; the lecture does not mandate numeric line weights. The scientific workflow uses blue, teal and navy with direct text labels.

## Figure placement for the next stage

`figure_slots.json` records exact positions and dimensions in millimetres. The main available frames are:

| Location | Intended content | Frame size |
| --- | --- | --- |
| A, upper row | Three total-protein interaction volcanoes | Each 171.1 × 66 mm |
| A, lower row | Three phosphosite interaction volcanoes | Each 171.1 × 66 mm |
| B | Compact KEGG GSEA summary, optional ORA support | 257.7 × 161 mm |
| C | Complete 75-node STRING network | 529.3 × 244 mm |
| D | TSB and PASN mutant Venn summaries | Each 124.8 × 88 mm |
| E | IsaA S194, Pyk T537 and QoxB T552 candidate cards | Each 81.9 × 114 mm |

These are starting positions, and every frame can be moved or resized. Remove the placeholder labels before placing figures. Use poster-specific plot layouts and labels at least 24 pt at their final physical size; simply shrinking the earlier figures would make many labels too small. The complete 75-node network may need fewer printed node labels or a revised layout, while retaining all nodes and linking the full labelled version through the supplement. If quantitative panels need more height, adjust the frame allocation instead of squeezing their text.

Keep simple phosphosite contrasts and PASN interactions distinct. The Venn set has 14 medium-specific hits representing 13 distinct sites on 12 proteins. None of the nine corresponding interaction tests for the three named candidate sites passes the significance threshold. A supplementary-index draft is in `../../reports/poster_supplement/README.md`; real QR codes should be inserted only after the intended online destinations are published and checked.

## Scientific sources for the background

[1] L. Schwyter et al., *mSystems*, vol. 10, e00262-25, 2025, doi: [10.1128/msystems.00262-25](https://doi.org/10.1128/msystems.00262-25).

[2] N. Prust et al., *Molecular & Cellular Proteomics*, vol. 20, 100034, 2021, doi: [10.1074/mcp.RA120.002232](https://doi.org/10.1074/mcp.RA120.002232).

[3] M. Huemer et al., *Science Signaling*, vol. 16, eabj8194, 2023, doi: [10.1126/scisignal.abj8194](https://doi.org/10.1126/scisignal.abj8194).

These are original research papers. The lecture, compilation and example posters were used for writing/layout guidance rather than copied experimental findings. The detailed sources and analysis distinctions also remain in the PowerPoint's speaker notes.
