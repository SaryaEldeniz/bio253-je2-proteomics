# Editable SVG illustrations for poster methods

The student requested the full project design as BioRender-like SVG illustrations
in a blue-green palette. Created original vector drawings in navy, blue, teal
and mint, with live text, labels and shapes supplementing colour. No BioRender
assets, external images or embedded bitmaps are used.

Verified the design against the original course Proteomics-Data README and the
configured sample sheet: 27 unique biological samples, three genotypes by three
conditions by three replicates. Marked transposon insertion in stk1 or stp1;
included TSB,15min pulse and approximately15h SN. Long exposure uses70%TSB+
30%cell-free sterile-filtered PAO1 spent medium. The pulse mixture was not
independently specified, so only its verified timing is labelled.

Panel01 illustrates PAO1 culture/spent medium, sterile filtration and separate
JE2 exposure, then the9-group/27-replicate matrix. PAO1 rods remain upstream of
filtration; the JE2 cells are cocci. Each matrix sample has a unique SVG ID.

Panel02 illustrates generic protein processing and the two complementary
measurement layers, total peptides versus Fe-IMAC-enriched phosphopeptides,
each followed by LC-MS/MS. Both layers concern the same27samples. STYH is labelled
as a search setting. No digestion enzyme, instrument model, extraction chemistry
or processing duration was inferred from another experiment.

Panel03 separates total-protein interactions from within-medium phosphosite
comparisons. The STRING branch comes directly from protein contrasts, bypassing
enrichment, to avoid implying selection by NES. The completed integration is
parent-protein/site comparison and prioritisation, not occupancy modelling or
functional validation. Workflow symbols resembling plots/logos/networks are
schematic, not newly computed results.

Files: steps/10_integration_poster/figures/methods_svg/01_experimental_design.svg,
02_proteomics_layers.svg,03_analysis_workflow.svg,00_complete_methods.svg.
Individual panels are transparent; combined version has a removable white
background. PNG previews and a ZIP with all4SVGs and source/editing notes are
included. R generators are in steps/10_integration_poster/scripts/methods_illustrations/; optional PNG rendering
uses Node/sharp without altering vector originals.

Checks passed: sample-sheet counts,27replicate symbols, valid XML, unique IDs,
resolved arrow markers and live text; no image/script/foreignObject elements.
Visually inspected all three panels and combined preview. Fixed initial genotype
label crowding and clarified the STYH search label. The SVGs remain editable
and scale independently of PNG resolution. Existing analyses/figures unchanged;
no GitHub publication attempted in this session.

AI assistance: Codex authored design/analysis vectors and assembly/checks. One
agent independently verified design sources and samples, another checked method
branch logic, and a third authored and inspected the proteomics acquisition
illustration. All artwork was generated specifically for this project.

## 2026-11-21 · Approved icon-led project-analysis panel

Sarya requested a more compact Project analysis graphic with larger illustrations and less text. Following her review, Codex added the three protein interaction comparisons (stk vs anc, stp vs anc, stp vs stk) and replaced the JE2 file symbol with an A/T/S/P methods symbol matching her reference image. Sarya then explicitly authorised uploading this revised version to GitHub.

The approved output is `figures/methods_svg/03_analysis_workflow_visual.svg`, with its PNG preview and `scripts/methods_illustrations/05_render_visual_analysis_svg.R`. The methods guide links this poster version prominently; the original detailed artwork remains unchanged. The earlier text-heavy compact draft is not part of this publication.

The illustration retains separate protein and phosphosite branches, direct contrast-based STRING selection bypassing enrichment, within-medium phosphosite comparisons and comparison of parent-protein/site signals. Detailed parameters and the separate phosphosite interaction tests remain documented in the original panel. The icons and letter heights are schematic method symbols, not measured plots, sequences or motifs. No biological data or selection results were changed.

Codex generated the SVG using the existing R palette/icon helpers without importing raw sample data, rendered it to PNG and inspected the layout. Checked: editable SVG text; three comparison labels; unique IDs and resolved arrow markers; no embedded image, script or foreignObject elements. SVG, preview, reproduction code and documentation were prepared together for the authorised GitHub update.
