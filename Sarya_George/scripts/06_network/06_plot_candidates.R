# Step 6 evidence plots: protein interaction estimates and measured sample values.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
prepared <- file.path(REPO_ROOT,"results/step06/step06_prepared.rds")
if (!file.exists(prepared)) stop("Run scripts/06_network/06_prepare_candidates.R first.")
d <- readRDS(prepared)
out <- file.path(REPO_ROOT,"figures/step06")
report <- file.path(REPO_ROOT,"reports/step06")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
codes <- c("stk_anc","stp_anc","stp_stk")
contrast_labels <- c(stk_anc="stk / anc",stp_anc="stp / anc",stp_stk="stp / stk")
direction_caption <- "Positive: SN-TSB response more positive in anc (stk/anc, stp/anc) or stk (stp/stk).\nNegative: converse. Effects are differences of responses, not absolute SN abundance."
palette <- c(anc="#52677F",stk="#0072B2",stp="#D55E00")
save_plot <- function(p,stem,width,height) {
  ggsave(file.path(out,paste0(stem,".png")),p,width=width,height=height,dpi=200,bg="white",limitsize=FALSE)
  ggsave(file.path(out,paste0(stem,".pdf")),p,width=width,height=height,bg="white",limitsize=FALSE)
}
labels <- bind_rows(d$candidates |> distinct(protein_key,display_name,locus_tag,uniprot_acc),
  select(d$anchors,protein_key,display_name,locus_tag,uniprot_acc)) |> distinct(protein_key,.keep_all=TRUE) |>
  mutate(protein_label=if_else(display_name==locus_tag,paste0(locus_tag," | ",uniprot_acc),
    paste0(display_name," | ",uniprot_acc)))

# Keep original values in source tables. Display NA where precise quantitative
# interpretation is unsuitable; zero is a valid measured effect, not missing.
selected_rows <- d$interactions |> semi_join(filter(d$candidates,selected),by=c("protein_key","term_id")) |>
  left_join(select(labels,protein_key,protein_label),by="protein_key",relationship="many-to-one") |>
  mutate(display_effect=if_else(observed_model & coalesce(all_relevant_groups_2of3,FALSE) & is.finite(diff),diff,NA_real_),
    cell_label=if_else(is.finite(display_effect),paste0(sprintf("%+.2f",display_effect),if_else(FDR<.1," *","")),"NA"))
stopifnot(all(abs(selected_rows$display_effect)<=2 | is.na(selected_rows$display_effect)))
for (i in seq_len(nrow(d$modules))) {
  mod <- d$modules[i,]
  tab <- selected_rows |> filter(term_id==mod$term_id)
  order <- filter(d$candidates,selected,term_id==mod$term_id)$protein_key
  axis_labels <- setNames(labels$protein_label,labels$protein_key)
  tab <- tab |> mutate(protein_key=factor(protein_key,levels=rev(order)),contrast_code=factor(contrast_code,levels=codes))
  p <- ggplot(tab,aes(contrast_code,protein_key,fill=display_effect)) +
    geom_tile(colour="white",linewidth=.7) + geom_text(aes(label=cell_label),size=3.3) +
    scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,
      limits=c(-2,2),na.value="#C8CDD2",name="Interaction\neffect (log2)") +
    scale_x_discrete(labels=contrast_labels,position="top") + scale_y_discrete(labels=axis_labels) +
    labs(title=mod$module,subtitle="Selected by protein contrasts; no NES or leading-edge filter",
      x=NULL,y=NULL,caption=paste0("* Supplied protein FDR < 0.1. All three estimates shown, including nonsignificant ones.\n",direction_caption)) +
    theme_minimal(base_size=11) + theme(panel.grid=element_blank(),legend.position="bottom",
      axis.text.y=element_text(size=9),plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=8))
  save_plot(p,paste0("interactions_",mod$file_stem),9.6,max(5,2.7+nrow(tab)/3*.24))
}

# Stable halo-ready node table for eventual STRING network joins by UniProt.
node_source <- d$base |> semi_join(bind_rows(filter(d$candidates,selected) |> select(protein_key),
  select(d$anchors,protein_key)),by="protein_key") |>
  mutate(halo_effect=if_else(observed_model & coalesce(all_relevant_groups_2of3,FALSE) & is.finite(diff),diff,NA_real_),
    halo_FDR=if_else(is.finite(halo_effect),FDR,NA_real_)) |>
  select(protein_key,contrast_code,diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR) |>
  pivot_wider(names_from=contrast_code,values_from=c(diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR))
nodes <- read_csv(file.path(report,"string_uploads/upload_manifest.csv"),show_col_types=FALSE) |>
  left_join(node_source,by="protein_key",relationship="one-to-one")
write_csv(nodes,file.path(report,"08_network_node_attributes.csv"))

# Individual dots are measured normalized log2 abundances. Means are not paired
# biological samples; there are no lines connecting replicate IDs across groups.
group_levels <- c("anc_TSB","anc_SN","stk_TSB","stk_SN","stp_TSB","stp_SN")
group_labels <- c("anc\nTSB","anc\nSN","stk\nTSB","stk\nSN","stp\nTSB","stp\nSN")
sample_values <- d$values |> left_join(select(labels,protein_key,protein_label),by="protein_key",relationship="many-to-one") |>
  mutate(plot_group=factor(paste(genotype,medium,sep="_"),levels=group_levels))
stats <- sample_values |> group_by(protein_key,protein_label,genotype,medium,plot_group) |>
  summarise(n=sum(is.finite(normalized_abundance)),
    mean=if(any(is.finite(normalized_abundance))) mean(normalized_abundance,na.rm=TRUE) else NA_real_,
    SD=if(sum(is.finite(normalized_abundance))>1) sd(normalized_abundance,na.rm=TRUE) else NA_real_,
    range=if(any(is.finite(normalized_abundance))) diff(range(normalized_abundance,na.rm=TRUE)) else NA_real_,.groups="drop")
write_csv(stats |> select(-plot_group),file.path(report,"09_replicate_measurement_summary.csv"))

plot_sample_page <- function(keys,title) {
  lv <- labels$protein_label[match(keys,labels$protein_key)]
  x <- sample_values |> filter(protein_key %in% keys) |> mutate(protein_label=factor(protein_label,levels=lv))
  s <- stats |> filter(protein_key %in% keys) |> mutate(protein_label=factor(protein_label,levels=lv))
  ggplot(x,aes(plot_group,normalized_abundance,colour=genotype)) +
    geom_point(position=position_jitter(width=.12,height=0,seed=611),size=2,alpha=.85,na.rm=TRUE) +
    geom_point(data=s,aes(y=mean),shape=95,size=7,na.rm=TRUE) +
    geom_text(data=s,aes(y=Inf,label=paste0(n,"/3")),vjust=1.2,size=2.5,show.legend=FALSE) +
    facet_wrap(~protein_label,ncol=2,scales="free_y",drop=FALSE) +
    scale_x_discrete(labels=setNames(group_labels,group_levels),drop=FALSE) + scale_colour_manual(values=palette) +
    scale_y_continuous(expand=expansion(mult=c(.12,.28))) +
    labs(title=title,subtitle="TSB and long SN exposure; SNpulse excluded",
      x=NULL,y="Measured normalized abundance (log2)",
      caption="Dots: individual measured samples; bars: observed group means; labels: detected / 3 replicates.\nMissing measurements remain NA. Each protein has its own y-axis; these are not phosphosite values.") +
    theme_bw(base_size=10) + theme(legend.position="bottom",panel.grid.minor=element_blank(),
      strip.text=element_text(face="bold",size=9),axis.text.x=element_text(size=8),
      plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=8))
}
plot_index <- list()
make_pages <- function(keys,title,stem,kind) {
  pages <- split(keys,ceiling(seq_along(keys)/6L))
  for (i in seq_along(pages)) {
    p <- plot_sample_page(pages[[i]],paste0(title," | ",i,"/",length(pages)))
    file <- paste0(stem,"_page",sprintf("%02d",i))
    save_plot(p,file,11,9)
    plot_index[[length(plot_index)+1L]] <<- tibble(kind=kind,title=title,page=i,
      file=paste0("../../figures/step06/",file,".png"),protein_keys=paste(pages[[i]],collapse=";"))
  }
}
for (i in seq_len(nrow(d$modules))) {
  mod <- d$modules[i,]
  keys <- filter(d$candidates,selected,term_id==mod$term_id)$protein_key
  make_pages(keys,mod$module,paste0("sample_values_",mod$file_stem),"Selected module candidates")
}
make_pages(d$anchors$protein_key,"Stk1 and Stp1: biological anchors","sample_values_anchors","Biological anchors")
held <- d$candidates |> filter(annotation_review,n_supported_interactions>0) |> distinct(protein_key) |> pull(protein_key)
make_pages(held,"Contrast-supported: module annotation review","sample_values_annotation_review","Held for membership review")
index <- bind_rows(plot_index)
write_csv(index,file.path(report,"sample_plot_index.csv"))
md <- c("# Step 6 measured sample plots","",
  "Every selected protein and both biological anchors have a measured-value panel. Separate pages contain contrast-supported proteins held for module-membership review.","",
  "Only TSB and SN are shown. These are normalized individual measurements, not model-imputed values. Counts show missing detections. Bars are group means, not paired measurements. Inspect within-protein patterns and missingness; y-axis ranges vary between proteins.","")
for (i in seq_len(nrow(index))) md <- c(md,paste0("## ",index$title[i]," — page ",index$page[i]),"",
  paste0("![",index$title[i],"](",index$file[i],")"),"")
writeLines(md,file.path(report,"sample_plots.md"))
message(nrow(index)," sample-value pages plus four interaction heatmaps. No STRING network requested.")
