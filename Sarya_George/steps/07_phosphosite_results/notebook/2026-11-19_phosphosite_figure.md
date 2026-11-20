# Visualise baseline phosphosite differences and supernatant interaction tests

The student asked how to show the abstract statement about IsaA S194, Pyk T537
and QoxB T552. Created a six-panel figure with observed phosphosite measurements
above and source interaction estimates below. PNG, PDF and SVG are saved under
`steps/06_string_network/figures/phosphosite_overlap/three_sites_control_vs_interaction.*`.

Panel A displays all 54 sample slots (three sites, six TSB/SN genotype groups,
three replicates); lines connect group means, not individual paired samples.
Missing values remain missing. IsaA has 2/3 ancestor detections in both media;
QoxB has 2/3 stp_SN detections; other displayed groups have 3/3. TSB labels show
the three selected significant mutant-versus-ancestor source FDRs.

Panel B displays all nine source interactions and their source confidence limits.
Verified each interval against estimate +/- t(0.975,df)*standard_error before
labelling it an unadjusted 95% interval. No source model, effect or FDR changed.
Every interaction FDR is >=0.1. QoxB stp/ancestor has a nominal interval excluding
zero but FDR0.616; both are displayed and their different correction scope is
explicit. Significance is judged by FDR, not by the unadjusted interval alone.

The figure does not claim no supernatant response, genotype equivalence,
phosphorylation occupancy or direct Stk1/Stp1 substrates. The measurements are
not adjusted for parent-protein abundance. Localization limitations remain in
the accompanying figure notes and original overlap report.

Reproduce from the group project:
`source("steps/07_phosphosite_results/scripts/06_plot_overlap_phosphosite_tests.R")`.
The script exports the nine plotted estimates and a detailed caption alongside
the existing overlap evidence. Root visually inspected the PNG: all six panels,
sample counts, effect intervals and FDR labels are visible without clipping.
Source/interval/detection assertions passed; only existing locale/package-build
warnings were emitted. No Git publication was attempted in this session.

AI assistance: Codex wrote and ran the R script, checked the figure and supplied
interpretation guidance. No sub-agents were used.
