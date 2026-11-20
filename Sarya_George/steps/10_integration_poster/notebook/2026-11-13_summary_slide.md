# 2026-11-13: one-slide summary through Step 5

Created one PowerPoint slide with three native editable total-proteome volcano
charts, a native KEGG GSEA heatmap table and a native GO BP ORA bar chart.
All three volcanoes use the same 1,858-protein observed-only mapped universe
as the displayed enrichment. Hit counts are stk/anc +51/-27, stp/anc +29/-32,
and stp/stk +25/-60. Source effects/FDRs are preserved; no source models were fit.

The KEGG examples cover purine metabolism, oxidative phosphorylation, nitrogen
metabolism and aminoacyl-tRNA biosynthesis. The ORA examples are cytolysis,
branched-chain amino-acid biosynthesis and de novo IMP biosynthesis from the
stk/ancestor interaction. These are selected examples for a one-slide overview;
the uncapped complete results remain unchanged. Notes document exact selection,
directions, source files and the study-wide correction family.

The slide states the main purine/nitrogen findings and the caveat that the stp
strain also has reduced Stk1. Native chart coordinates use Excel-compatible
15-significant-digit display precision; original saved statistics retain full
precision. No imputed sample measurements were plotted.

Evidence preparation is saved in scripts/05_figures/05_prepare_slide_summary.R.
Native slide authoring is in steps/10_integration_poster/scripts/methods_illustrations/05_build_summary_slide.mjs,
used in the connected Artifact Session artifact_5so-FneVUotbNIJTDT_DxRHm6D80cBy3.
The Windows Node finalizer cannot canonicalize this workspace (EPERM), so the
export script runs the same bundled package/layout/font/native-chart/title/table
validators directly. PowerPoint mode validates complete native literal chart
data; native PowerPoint itself was not opened. The connected slide PNG was
visually inspected, including the final table-margin and bar-order corrections.

The deliverable is presentations/JE2_PASN_progress_summary.pptx. Preview PNG
is alongside it. Detailed chart-coordinate JSON and validation records stay in
ignored results/slide_summary/. Codex authored the slide, and a separate agent
prepared and independently checked the R evidence extraction.
