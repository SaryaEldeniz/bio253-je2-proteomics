# Show the baseline differences and the separate genotype-by-medium tests.
# Read source estimates; do not fit models, adjust FDRs again or impute values.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
report <- file.path(REPO_ROOT,"steps/06_string_network/reports/phosphosite_overlap")
fig <- file.path(REPO_ROOT,"steps/07_phosphosite_results/figures/baseline_vs_interaction")
hits <- read_csv(file.path(report,"three_overlapping_sites_and_parent_TSB.csv"),show_col_types=FALSE)
interactions <- read_csv(file.path(report,"three_sites_PASN_interaction_estimates.csv"),show_col_types=FALSE)
stopifnot(nrow(hits)==3L,nrow(interactions)==9L,
  all(hits$site_estimate_type=="observed"),all(hits$site_FDR<.1),
  all(interactions$estimate_type=="observed"),all(interactions$FDR>=.1))
# Check the meaning of the supplied confidence limits before labelling them.
halfwidth <- qt(.975,interactions$df)*interactions$std.error
stopifnot(all(abs(interactions$conf.low-(interactions$diff-halfwidth))<1e-7),
  all(abs(interactions$conf.high-(interactions$diff+halfwidth))<1e-7))
labels <- c("IsaA S194","Pyk T537","QoxB T552")
make_label <- function(protein,residue,position) factor(paste(protein,paste0(residue,position)),levels=labels)
hits <- hits |> mutate(site_label=make_label(protein,residue,position),
  baseline_label=paste0("TSB: ",mutant," > anc; FDR = ",formatC(site_FDR,digits=3,format="fg")))
abundance <- read_tsv(file.path(DATA_ROOT,"derived/abundances_phospho_STYH.tsv.gz"),show_col_types=FALSE) |>
  semi_join(hits,by="site_id") |> filter(medium %in% c("TSB","SN")) |>
  left_join(hits |> select(site_id,protein,site_label),by="site_id",relationship="many-to-one") |>
  mutate(genotype=factor(genotype,levels=c("anc","stk","stp")),medium=factor(medium,levels=c("TSB","SN")))
stopifnot(nrow(abundance)==54L,!anyDuplicated(abundance[c("site_id","sample")]))
counts <- abundance |> group_by(site_label,genotype,medium) |>
  summarise(n_detected=sum(is.finite(normalized_abundance)),.groups="drop")
stopifnot(all(counts$n_detected %in% 2:3))
palette <- c(anc="#333333",stk="#0072B2",stp="#D55E00")
theme <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill="#F0F3F5",colour="#BCC5CB"),
  strip.text=element_text(face="bold"),plot.title=element_text(face="bold",size=13),
  plot.subtitle=element_text(size=10),legend.position="top",legend.title=element_blank(),
  plot.margin=margin(5,10,5,10))
top <- ggplot(abundance,aes(medium,normalized_abundance,colour=genotype,shape=genotype,group=genotype)) +
  stat_summary(fun=mean,geom="line",linewidth=.6,position=position_dodge(.45),na.rm=TRUE,show.legend=FALSE) +
  geom_point(position=position_jitterdodge(jitter.width=.07,dodge.width=.45,seed=19),size=2.3,alpha=.9,na.rm=TRUE) +
  stat_summary(fun=mean,geom="point",shape=95,size=6,position=position_dodge(.45),na.rm=TRUE,show.legend=FALSE) +
  geom_text(data=counts,aes(y=-Inf,label=paste0(n_detected,"/3")),
    position=position_dodge(.45),vjust=-.4,size=3,show.legend=FALSE) +
  geom_text(data=hits,aes(x=1.5,y=Inf,label=baseline_label),inherit.aes=FALSE,vjust=1.4,size=3.2,colour="#202C35") +
  facet_wrap(~site_label,nrow=1,scales="free_y") +
  scale_colour_manual(values=palette) + scale_shape_manual(values=c(anc=16,stk=17,stp=15)) +
  scale_y_continuous(expand=expansion(mult=c(.14,.21))) +
  labs(title="A  |  Mutant-ancestor differences in control broth",
    subtitle="Significant TSB comparisons are labelled; SN samples show the exposure context.",
    x=NULL,y="Normalized phosphosite signal (log2)") + theme
contrast_codes <- c(interaction_stk_vs_anc_at_TSB_vs_SN="stk - anc",
  interaction_stp_vs_anc_at_TSB_vs_SN="stp - anc",interaction_stp_vs_stk_at_TSB_vs_SN="stp - stk")
interactions <- interactions |> mutate(site_label=make_label(protein,residue,position),
  comparison=factor(unname(contrast_codes[contrast]),levels=c("stp - stk","stp - anc","stk - anc")),
  fdr_label=paste0("FDR = ",formatC(FDR,digits=3,format="f")))
stopifnot(!anyNA(interactions$comparison))
bottom <- ggplot(interactions,aes(diff,comparison)) +
  geom_vline(xintercept=0,linetype="dashed",colour="#78858D",linewidth=.5) +
  geom_segment(aes(x=conf.low,xend=conf.high,yend=comparison),linewidth=.85,colour="#384D60") +
  geom_point(shape=21,size=3.1,fill="#FFFFFF",colour="#243B50",stroke=1) +
  geom_text(aes(x=4.05,label=fdr_label),hjust=0,size=3.2,colour="#202C35") +
  facet_wrap(~site_label,nrow=1) + scale_x_continuous(limits=c(-3.25,5.65),breaks=c(-2,0,2,4)) +
  scale_y_discrete(expand=expansion(add=.65)) +
  labs(title="B  |  Does the supernatant response differ between genotypes?",
    subtitle="None of these nine phosphosite interaction tests passes FDR < 0.1.",
    x="Source interaction effect (log2 units)",y=NULL) + theme + theme(legend.position="none")
draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.text("Baseline site differences do not establish different responses to supernatant",
    x=.02,y=.983,just=c("left","top"),gp=grid::gpar(fontsize=14,fontface="bold",col="#172E40"))
  print(top,vp=grid::viewport(x=.5,y=.685,width=1,height=.54))
  print(bottom,vp=grid::viewport(x=.5,y=.245,width=1,height=.345))
  grid::grid.text(paste(
    "TSB: control broth. SN: Pseudomonas supernatant. Points: observed samples; lines: group means; counts: detected/3.",
    "Interaction = (first - second) in TSB minus (first - second) in SN. Positive: larger genotype gap in TSB; negative: larger gap in SN.",
    "Intervals: unadjusted 95% confidence intervals; significance uses source FDR < 0.1. Signals are not adjusted for parent-protein abundance.",
    sep="\n"),x=.022,y=.014,just=c("left","bottom"),gp=grid::gpar(fontsize=8.2,col="#425564",lineheight=1.16))
}
stem <- file.path(fig,"three_sites_control_vs_interaction")
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
png(paste0(stem,".png"),width=12,height=8.8,units="in",res=300,bg="white")
draw_figure(); invisible(dev.off())
pdf(paste0(stem,".pdf"),width=12,height=8.8,useDingbats=FALSE)
draw_figure(); invisible(dev.off())
svg(paste0(stem,".svg"),width=12,height=8.8)
draw_figure(); invisible(dev.off())
write_csv(interactions |> select(protein,residue,position,contrast,diff,conf.low,conf.high,FDR,estimate_type),
  file.path(report,"figure_control_vs_interaction_estimates.csv"))
writeLines(c("# Baseline phosphosite differences and interaction tests","",
  "Panel A displays observed source-normalized phosphosite sample values for all three genotypes in TSB and SN. Each line connects group means, not paired individual cells/samples. FDR labels apply only to the specified mutant-versus-ancestor contrast in TSB. Missing values remain NA, with detection counts shown.","",
  "Panel B displays all three source genotype-by-medium interactions for each site. Points are estimates; bars are source unadjusted95% confidence intervals. FDR labels, not whether a bar crosses zero, determine significance at FDR<0.1. QoxB stp/anc illustrates this distinction: its nominal confidence interval excludes zero but its source FDR is0.616 after multiple-testing correction.","",
  "Interaction=(first_TSB-second_TSB)-(first_SN-second_SN), using log2 phosphosite signal. Positive indicates a larger first-minus-second genotype gap in TSB; negative indicates a larger gap in SN. No new model, confidence interval or FDR was fitted. Supplied intervals were verified against t(0.975,source_df)*source_standard_error.","",
  "The findings concern three specified sites only. They do not establish no response to supernatant, equivalence between genotypes, direct kinase/phosphatase targeting or changes independent of parent-protein abundance. IsaA and QoxB have variable source localization probabilities; investigate site assignments before functional interpretation.","",
  "Suggested legend: Observed phosphosite signals for IsaA S194, Pyk T537 and QoxB T552 show significant mutant-ancestor differences in control broth (A). Corresponding genotype-by-medium interaction estimates provide insufficient evidence for genotype-dependent responses to supernatant at source FDR<0.1 (B). Lines inA connect group means; error bars inB are unadjusted95% confidence intervals. Phosphosite signals are not adjusted for parent-protein abundance.","",
  "[PNG](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.png) | [PDF](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.pdf) | [SVG](../../../07_phosphosite_results/figures/baseline_vs_interaction/three_sites_control_vs_interaction.svg)","",
  "Reproduce: source(\"steps/07_phosphosite_results/scripts/06_plot_overlap_phosphosite_tests.R\")."),
  file.path(report,"figure_control_vs_interaction.md"))
message("Created3-format figure:54 sample slots,three selected significant TSB comparisons,nine source interaction estimates; none interaction-significant.")
