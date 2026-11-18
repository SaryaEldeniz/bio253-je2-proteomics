# Step 6 part 7: review network candidates against experimental contrasts.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"reports/step06/string_network")
fig <- file.path(REPO_ROOT,"figures/step06/string_network")
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
nodes <- read_csv(file.path(out,"nodes.csv"),show_col_types=FALSE)
edges <- read_csv(file.path(out,"edges_expanded.csv"),show_col_types=FALSE)
source <- read_tsv(file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),show_col_types=FALSE)
codes <- c("stk_anc","stp_anc","stp_stk")
env <- source |> filter(contrast %in% paste0("TSB_vs_SN_at_",c("anc","stk","stp"))) |>
  transmute(protein_key,genotype=sub("TSB_vs_SN_at_","",contrast),SN_minus_TSB=-diff) |>
  pivot_wider(names_from=genotype,values_from=SN_minus_TSB,names_prefix="R_")
nodes <- nodes |> left_join(env,by="protein_key",relationship="many-to-one")
anchors <- filter(nodes,node_role=="Biological anchor")$string_id
anchor_edges <- edges |> filter(node_a %in% anchors | node_b %in% anchors) |>
  mutate(name_a=nodes$display_name[match(node_a,nodes$string_id)],
    name_b=nodes$display_name[match(node_b,nodes$string_id)],
    high_confidence=combined_score>=.7,
    evidence_note=case_when(experimental>0~"Nonzero experimental channel; not direct-substrate proof",
      database>0~"Curated database channel present; inspect underlying source",
      textmining>0~"Text-mining contributes; functional context only",
      TRUE~"Other computational evidence; functional context only"))
write_csv(anchor_edges,file.path(out,"anchor_links_review.csv"))
for(code in codes) nodes[[paste0("supported_",code)]] <-
  is.finite(nodes[[paste0("halo_effect_",code)]]) & is.finite(nodes[[paste0("halo_FDR_",code)]]) &
  nodes[[paste0("halo_FDR_",code)]]<.1
nodes$n_supported_contrasts <- rowSums(as.data.frame(nodes[paste0("supported_",codes)]))
nodes$large_effect_supported <- Reduce(`|`,lapply(codes,function(code)
  nodes[[paste0("supported_",code)]] & abs(nodes[[paste0("halo_effect_",code)]])>1))
nodes$experimental_review <- case_when(nodes$node_role=="Biological anchor"~"Biological anchor; do not require abundance significance",
  nodes$large_effect_supported~"Prioritize contrast interpretation: supported large effect",
  nodes$n_supported_contrasts>0~"Inspect supported smaller effect",
  TRUE~"Network context only under current protein FDR rule")
nodes$membership_review <- if_else(nodes$annotation_review,
  "Retain user selection; resolve supplied pathway membership before pathway-specific claims",
  "No specific membership flag; not an independent annotation validation")
write_csv(nodes,file.path(out,"final_protein_review.csv"))

# Plots are needed before interpreting newly added protein abundance effects.
added <- filter(nodes,node_role=="Added partner")
values <- read_abundances("total") |> filter(medium %in% c("TSB","SN")) |>
  semi_join(added,by="protein_key") |> left_join(select(added,protein_key,display_name,uniprot_acc),by="protein_key",relationship="many-to-one") |>
  mutate(protein_label=paste0(display_name," | ",uniprot_acc),
    plot_group=factor(paste(genotype,medium,sep="_"),levels=c("anc_TSB","anc_SN","stk_TSB","stk_SN","stp_TSB","stp_SN")))
stats <- values |> group_by(protein_key,protein_label,plot_group,genotype) |>
  summarise(n=sum(is.finite(normalized_abundance)),
    mean=if(any(is.finite(normalized_abundance))) mean(normalized_abundance,na.rm=TRUE) else NA_real_,.groups="drop")
write_csv(stats,file.path(out,"added_partner_measurement_summary.csv"))
stopifnot(nrow(values)==18L*nrow(added),all(added$protein_key %in% stats$protein_key))
pages <- split(added$protein_key,ceiling(seq_len(nrow(added))/6))
for(i in seq_along(pages)) {
  x <- filter(values,protein_key %in% pages[[i]])
  s <- filter(stats,protein_key %in% pages[[i]])
  p <- ggplot(x,aes(plot_group,normalized_abundance,colour=genotype)) +
    geom_point(position=position_jitter(width=.12,height=0,seed=617),size=2,na.rm=TRUE) +
    geom_point(data=s,aes(y=mean),shape=95,size=7,na.rm=TRUE) +
    geom_text(data=s,aes(y=Inf,label=paste0(n,"/3")),vjust=1.2,size=2.5,show.legend=FALSE) +
    facet_wrap(~protein_label,ncol=2,scales="free_y") +
    scale_colour_manual(values=c(anc="#52677F",stk="#0072B2",stp="#D55E00")) +
    scale_x_discrete(labels=c("anc\nTSB","anc\nSN","stk\nTSB","stk\nSN","stp\nTSB","stp\nSN"),drop=FALSE) +
    scale_y_continuous(expand=expansion(mult=c(.12,.28))) + theme_bw(base_size=10) +
    labs(title=paste0("Added STRING partners | page ",i,"/",length(pages)),
      subtitle="Selected from network context; measured TSB/SN values shown independently",
      x=NULL,y="Measured normalized abundance (log2)",
      caption="Dots: measured samples; bars: observed means; n/3: detections. Missing values remain NA.\nFree y-axis per protein. STRING association does not prove a direct kinase/phosphatase substrate.") +
    theme(legend.position="bottom",panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0,size=8))
  stem <- file.path(fig,paste0("added_partner_samples_page",sprintf("%02d",i)))
  ggsave(paste0(stem,".png"),p,width=11,height=9,dpi=200,bg="white")
  ggsave(paste0(stem,".pdf"),p,width=11,height=9,bg="white")
}

# A compact, readable review keeps every original selection, not just hits.
md <- c("# Protein review after one-round STRING expansion","",
  "All 49 original seeds/anchors are retained. Twelve partners were added for module coverage and network context. Protein-contrast evidence is evaluated separately; multiple significant interactions are not independent replication.","",
  "| Protein | UniProt | Role | Supported contrasts | Review |",
  "|---|---|---|---|---|")
for(i in seq_len(nrow(nodes))) {
  r <- nodes[i,]
  supported <- codes[as.logical(unlist(r[paste0("supported_",codes)]))]
  md <- c(md,paste0("| ",r$display_name,if_else(r$annotation_review," †","")," | ",r$uniprot_acc,
    " | ",r$node_role," | ",if(length(supported)) paste(supported,collapse=", ") else "None under this rule",
    " | ",r$experimental_review," |"))
}
md <- c(md,"","† Supplied pathway membership flagged for review; protein remains included as selected.","",
  "Supported means source-observed quantitative halo, >=2/3 relevant detections, supplied protein FDR<0.1. Large additionally means |interaction effect|>1. Missing quantitative halos are not zero effects. Original effects, FDRs, detection counts and within-genotype responses are in final_protein_review.csv.")
writeLines(md,file.path(out,"protein_review_readable.md"))
message("Reviewed ",nrow(nodes)," proteins; added-partner sample panels: ",length(pages),".")
print(nodes |> count(node_role,experimental_review),width=140)
