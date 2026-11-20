# Baseline phosphosite differences and interaction tests

Panel A displays observed source-normalized phosphosite sample values for all three genotypes in TSB and SN. Each line connects group means, not paired individual cells/samples. FDR labels apply only to the specified mutant-versus-ancestor contrast in TSB. Missing values remain NA, with detection counts shown.

Panel B displays all three source genotype-by-medium interactions for each site. Points are estimates; bars are source unadjusted95% confidence intervals. FDR labels, not whether a bar crosses zero, determine significance at FDR<0.1. QoxB stp/anc illustrates this distinction: its nominal confidence interval excludes zero but its source FDR is0.616 after multiple-testing correction.

Interaction=(first_TSB-second_TSB)-(first_SN-second_SN), using log2 phosphosite signal. Positive indicates a larger first-minus-second genotype gap in TSB; negative indicates a larger gap in SN. No new model, confidence interval or FDR was fitted. Supplied intervals were verified against t(0.975,source_df)*source_standard_error.

The findings concern three specified sites only. They do not establish no response to supernatant, equivalence between genotypes, direct kinase/phosphatase targeting or changes independent of parent-protein abundance. IsaA and QoxB have variable source localization probabilities; investigate site assignments before functional interpretation.

Suggested legend: Observed phosphosite signals for IsaA S194, Pyk T537 and QoxB T552 show significant mutant-ancestor differences in control broth (A). Corresponding genotype-by-medium interaction estimates provide insufficient evidence for genotype-dependent responses to supernatant at source FDR<0.1 (B). Lines inA connect group means; error bars inB are unadjusted95% confidence intervals. Phosphosite signals are not adjusted for parent-protein abundance.

[PNG](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.png) | [PDF](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.pdf) | [SVG](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.svg)

Reproduce: source("steps/07_phosphosite_results/scripts/06_plot_overlap_phosphosite_tests.R").
