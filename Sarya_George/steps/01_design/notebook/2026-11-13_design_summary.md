# 2026-11-13: compact design and interaction figure

Created a minimal schematic following the user's dark-blue example. It shows
JE2 ancestor, JE2 stk mutant and JE2 stp mutant across TSB, SNpulse (15min PASN)
and SN (~15h PASN). Each of the nine groups contains three sample dots, giving
27 samples. Labels/counts were checked against read_samples(); no source data
or statistical results were changed.

Also added the requested workflow overview, omitting steps 2 and 11 while
retaining the original numbering. The diagram shows stages, without marking
all stages complete. R source: steps/01_design/scripts/01_render_project_workflow.R;
PNG/SVG: steps/01_design/figures/design/JE2_project_workflow.*. The output was visually checked.

The neighbouring panel lists the stk/ancestor, stp/ancestor and stp/stk genotype
interactions, all comparing TSB with SN. SNpulse is part of the full design but
outside these three focal interaction contrasts. PASN is labelled as cell-free
PAO1 supernatant.

Reproduce with scripts/05_figures/01_render_design_summary.R. PNG, SVG and PDF
are saved under steps/01_design/figures/design/. Codex used R/grid for exact text and dot counts.
