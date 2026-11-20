# Slide evidence from saved October observed-only results. No models are refit.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out<-file.path(REPO_ROOT,"steps/10_integration_poster/figures/slide_summary")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
data_out<-file.path(REPO_ROOT,"results","slide_summary")
dir.create(data_out,recursive=TRUE,showWarnings=FALSE)
contrasts<-c("interaction_stk_vs_anc_at_TSB_vs_SN","interaction_stp_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_stk_at_TSB_vs_SN")
titles<-c("stk vs ancestor","stp vs ancestor","stp vs stk")
cache_paths<-file.path(REPO_ROOT,"results","step05","observed_only","cache",paste0(contrasts,".rds"))
caches<-lapply(cache_paths,readRDS)
points<-bind_rows(lapply(seq_along(contrasts),function(i) {
  d<-caches[[i]]$prepared
  stopifnot(nrow(d)==1858L,all(d$estimate_type=="observed"),
    !anyDuplicated(d$target_id),all(is.finite(d$diff)),all(is.finite(d$FDR)))
  d |> transmute(contrast=contrasts[i],protein_key,target_id,diff,FDR,
    group=case_when(FDR<.1 & diff>1~"positive",FDR<.1 & diff< -1~"negative",TRUE~"other"))
}))
minimum_positive<-min(points$FDR[points$FDR>0])
display_floor<-max(.Machine$double.xmin,minimum_positive/10)
points<-points |> mutate(y=-log10(pmax(FDR,display_floor)))
x_limit<-ceiling(max(abs(points$diff))); y_limit<-ceiling(max(points$y))
series_order<-c("other","negative","positive")
volcanoes<-lapply(seq_along(contrasts),function(i) {
  d<-filter(points,contrast==contrasts[i])
  metadata<-caches[[i]]$metadata
  list(contrast=contrasts[i],title=titles[i],n=nrow(d),
    positive_direction=metadata$positive_direction[1],negative_direction=metadata$negative_direction[1],
    counts=as.list(setNames(vapply(series_order,function(g)sum(d$group==g),integer(1)),series_order)),
    series=lapply(series_order,function(g) {
      z<-filter(d,group==g)
      list(name=g,xvalues=unname(z$diff),yvalues=unname(z$y))
    }))
})
gsea<-read_tsv(file.path(REPO_ROOT,"results","step05","observed_only","gsea_all.tsv.gz"),show_col_types=FALSE)
ora<-read_tsv(file.path(REPO_ROOT,"results","step05","observed_only","ora_all.tsv.gz"),show_col_types=FALSE)
selected_kegg<-c("map00230","map00190","map00910","map00970")
kegg_labels<-c("Purine metabolism","Oxidative phosphorylation","Nitrogen metabolism","Aminoacyl-tRNA biosynthesis")
kegg<-lapply(seq_along(selected_kegg),function(j) {
  rows<-gsea |> filter(collection=="KEGG",term_id==selected_kegg[j],contrast %in% contrasts)
  rows<-rows[match(contrasts,rows$contrast),]
  stopifnot(nrow(rows)==3,!anyNA(rows$contrast),all(is.finite(rows$NES)))
  list(term_id=selected_kegg[j],term=kegg_labels[j],
    cells=lapply(seq_len(nrow(rows)),function(i)list(contrast=rows$contrast[i],NES=rows$NES[i],
      FDR=rows$FDR_study[i],significant=rows$FDR_study[i]<.1,direction=rows$direction[i])))
})
ora_ids<-c("GO:0009082","GO:0006189","GO:0019835")
ora_selected<-ora |> filter(collection=="GO_BP",contrast==contrasts[1],term_id %in% ora_ids,FDR_study<.1) |>
  arrange(desc(fold_enrichment))
stopifnot(nrow(ora_selected)==3,setequal(ora_selected$term_id,ora_ids))
ora_rows<-lapply(seq_len(nrow(ora_selected)),function(i) {
  row<-ora_selected[i,]
  list(term_id=row$term_id,term=row$term_description,fold_enrichment=row$fold_enrichment,
    overlap=row$overlap,effect_sign=row$effect_sign,FDR=row$FDR_study,direction=row$direction)
})
evidence<-list(source="October total proteome, saved observed_only analysis",
  policy="observed_only",n_proteins=1858L,protein_FDR_cutoff=.1,abs_log2_effect_cutoff=1,
  enrichment_FDR_cutoff=.1,x_limit=x_limit,y_limit=y_limit,display_floor=display_floor,
  zero_FDR_count=sum(points$FDR==0),volcanoes=volcanoes,ora=list(contrast=contrasts[1],rows=ora_rows),
  gsea=list(contrasts=contrasts,rows=kegg),
  selection_note=paste("The four KEGG examples represent nucleotide metabolism, respiration, nitrogen response and translation.",
    "Only cells with study FDR<0.1 should be highlighted. The three ORA examples are nonredundant observed-only GO BP discoveries.",
    "Study-level enrichment FDR includes ten comparisons and five original annotation collections separately by method."))
jsonlite::write_json(evidence,file.path(data_out,"evidence.json"),pretty=TRUE,auto_unbox=TRUE,digits=16,na="null")

# Compact white-background fallbacks. All 1,858 source-observed proteins remain.
colours<-c(other="#B4BEC9",negative="#0072B2",positive="#D55E00")
for(i in seq_along(contrasts)) {
  d<-filter(points,contrast==contrasts[i])
  p<-ggplot(d,aes(diff,y,colour=factor(group,levels=series_order))) +
    geom_vline(xintercept=c(-1,1),linetype="dashed",colour="#8A97A6",linewidth=.35) +
    geom_hline(yintercept=1,linetype="dashed",colour="#8A97A6",linewidth=.35) +
    geom_point(size=.85,alpha=.8) + scale_colour_manual(values=colours,guide="none") +
    scale_x_continuous(limits=c(-x_limit,x_limit),breaks=pretty(c(-x_limit,x_limit),5)) +
    scale_y_continuous(limits=c(0,y_limit),expand=expansion(mult=c(0,.02))) +
    labs(title=titles[i],subtitle=paste0("Positive ",volcanoes[[i]]$counts$positive,
      "; negative ",volcanoes[[i]]$counts$negative),
      x="Interaction effect (log2)",y="-log10(protein FDR)") +
    theme_classic(base_size=12) + theme(plot.title=element_text(face="bold",size=15),
      plot.subtitle=element_text(size=10),plot.margin=margin(7,9,7,7))
  ggsave(file.path(out,paste0("volcano_",i,".png")),p,width=4.7,height=3.4,dpi=300,bg="white")
}
writeLines(c("# Slide summary evidence", "",
  "All volcano coordinates are unchanged October source DEA effects/FDRs for the saved observed-only UniProt universe.",
  "Ignored results/slide_summary/evidence.json contains coordinates for native PowerPoint chart authoring; no raw sample measurements are exported.",
  evidence$selection_note,"",
  paste0("Common axes: x ±",x_limit,"; y 0 to ",y_limit,". Zero-FDR rows: ",sum(points$FDR==0),".")),
  file.path(out,"README.md"))
message("Slide evidence written: ",file.path(data_out,"evidence.json"))
print(bind_rows(lapply(volcanoes,function(v)tibble(contrast=v$contrast,other=v$counts$other,
  negative=v$counts$negative,positive=v$counts$positive))))
print(ora_selected |> select(term_description,fold_enrichment,overlap,FDR_study,effect_sign))
