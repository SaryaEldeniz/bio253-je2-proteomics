# Sample-level checks for selected module contributors; no enrichment is rerun.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
gsea<-read_enrichment("GSEA",collection="KEGG")
focus<-c("TSB_vs_SN_at_anc",step04_contrast_plan() |> filter(contrast_type=="interaction") |> pull(contrast))
selected<-gsea |> filter(contrast %in% focus) |> group_by(contrast) |>
  arrange(FDR_study,pathway,.by_group=TRUE) |> slice_head(n=2) |> ungroup()
leading<-selected |> select(leading_edge_ids) |> separate_rows(leading_edge_ids,sep=";") |>
  filter(nzchar(leading_edge_ids)) |> count(uniprot_acc=leading_edge_ids,sort=TRUE) |>
  arrange(desc(n),uniprot_acc) |> slice_head(n=6)
map<-read_identifier_map() |> filter(entry_class=="cds",eligible_uniprot_annotation,in_total_DEA,n_JE2_good_total_DEA==1)
candidate_keys<-unique(c(lookup_je2(c("Stk1","Stp1"))$protein_key,
  filter(map,uniprot_acc %in% leading$uniprot_acc)$protein_key))
nitrogen<-filter(gsea,contrast=="interaction_stp_vs_anc_at_TSB_vs_SN",term_id=="map00910")
nitrogen_ids<-if(nrow(nitrogen)) strsplit(nitrogen$leading_edge_ids[1],";",fixed=TRUE)[[1]] else character()
nitrogen_keys<-filter(map,uniprot_acc %in% nitrogen_ids)$protein_key
# Display the two media used by the long-exposure comparisons; keep the same proteins.
abundances<-read_abundances("total") |> filter(as.character(medium) %in% c("TSB","SN"))
figure_dir<-file.path(REPO_ROOT,"figures","step05"); dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
plot_two_media_evidence<-function(keys) {
  offsets<-c(anc=-.27,stk=0,stp=.27)
  x<-abundances |> filter(protein_key %in% keys) |>
    mutate(label=protein_label(locus_tag,product),
      x_pos=as.numeric(factor(as.character(medium),levels=c("TSB","SN")))+unname(offsets[as.character(genotype)]))
  if(!nrow(x)) stop("No measured-abundance rows for the selected proteins.")
  counts<-x |> group_by(label,genotype,medium,x_pos) |>
    summarise(n=sum(is.finite(normalized_abundance)),.groups="drop") |>
    mutate(text=paste0(n,"/3"))
  p<-ggplot(x,aes(x_pos,normalized_abundance,colour=genotype)) +
    geom_boxplot(aes(group=interaction(genotype,medium)),width=.22,outlier.shape=NA,na.rm=TRUE) +
    geom_point(position=position_jitter(width=.04,height=0,seed=1),na.rm=TRUE) +
    geom_text(data=counts,aes(label=text,y=Inf),vjust=1.3,size=2.7,show.legend=FALSE) +
    facet_wrap(~label,scales="free_y") +
    scale_x_continuous(breaks=c(1,2),labels=c("TSB","SN"),limits=c(.6,2.4)) +
    scale_colour_manual(values=GENOTYPE_COLOURS) +
    scale_y_continuous(expand=expansion(mult=c(.05,.2))) +
    labs(x="medium",y="normalised abundance (log2)")
  # Verify that removal is from plot data and axis, not just hidden tick labels.
  stopifnot(setequal(as.character(unique(x$medium)),c("TSB","SN")),
    all(x$x_pos>=.6 & x$x_pos<=2.4),n_distinct(x$protein_key)==length(unique(keys)))
  p
}
save_evidence<-function(keys,stem,title,subtitle,height) {
  if(!length(keys)) return(invisible(NULL))
  p<-plot_two_media_evidence(keys) +
    labs(title=title,subtitle=subtitle,
      caption="Points are individual samples; labels count measured values out of three. Missing values stay missing. Protein abundance does not establish flux or direct kinase targets.")
  ggsave(file.path(figure_dir,paste0(stem,".png")),p,width=15,height=height,dpi=300,bg="white")
  ggsave(file.path(figure_dir,paste0(stem,".pdf")),p,width=15,height=height,useDingbats=FALSE)
}
save_evidence(candidate_keys,"step05_candidate_sample_values",
  "Measured abundances of pathway contributors and Stk1/Stp1",
  "Six contributors chosen by leading-edge frequency in the top two significant KEGG terms per focus contrast; ties by accession",11)
save_evidence(nitrogen_keys,"step05_nitrogen_sample_values",
  "Measured values behind the stp-dependent nitrogen-metabolism result",
  "All leading-edge proteins of the study-significant stp-versus-ancestor nitrogen-metabolism interaction",8)
source_data<-read_tsv(file.path(REPO_ROOT,"results","step04","total","annotated_all.tsv.gz"),show_col_types=FALSE)
audit<-source_data |> filter(protein_key %in% unique(c(candidate_keys,nitrogen_keys)),contrast %in% focus) |>
  select(contrast,display_name,locus_tag,uniprot_acc,diff,FDR,estimate_type,group_detection_pattern)
write_tsv(audit,file.path(REPO_ROOT,"reports","step05","selected_candidate_measurement_support.tsv"))
manifest_path<-file.path(REPO_ROOT,"reports","step05_figure_manifest.tsv")
if(file.exists(manifest_path)) {
  stems<-c("step05_candidate_sample_values","step05_nitrogen_sample_values")
  manifest<-read_tsv(manifest_path,show_col_types=FALSE) |> filter(!figure %in% stems)
  extra<-tibble(figure=stems,method="Measured sample values",collection="KEGG contributors and Stk1/Stp1",
    scope=c("Selected focus contributors plus Stk1/Stp1; TSB and SN only","Nitrogen leading edge: stp-vs-ancestor interaction; TSB and SN only"),
    policy="primary-based selection",namespace="JE2",selection=c("Six highest leading-edge frequencies plus Stk1/Stp1","All nitrogen leading-edge proteins"),
    significant_candidates=c(length(candidate_keys),length(nitrogen_keys)),
    png=paste0("figures/step05/",stems,".png"),pdf=paste0("figures/step05/",stems,".pdf"))
  write_tsv(bind_rows(manifest,extra),manifest_path)
}
capture.output(sessionInfo(),file=file.path(REPO_ROOT,"reports","step05","figure_sessionInfo.txt"))
message("Module sample-evidence plots complete.")
