# Display every significant observed-only term, without an arbitrary Top-N cap.
# These helpers never change the enrichment statistics or correction family.
step05_all_plot_rows <- function(x, method = c("GSEA", "ORA"), fdr_cutoff = 0.1) {
  method <- match.arg(method)
  needed <- c("policy", "namespace", "contrast", "contrast_label", "collection",
              "pathway", "term_id", "term_description", "FDR_study", "direction")
  needed <- c(needed, if(method == "GSEA") c("NES", "size") else c("effect_sign", "fold_enrichment", "overlap"))
  if(!all(needed %in% names(x))) stop("All-results plot input is missing: ",paste(setdiff(needed,names(x)),collapse=", "))
  if(any(x$policy != "observed_only") || any(x$namespace != "UniProt")) {
    stop("These figures require the separately calculated observed_only UniProt analysis.")
  }
  if(n_distinct(x$contrast)>1 || n_distinct(x$collection)>1) stop("Choose one comparison and collection per plot.")
  x <- x |> filter(is.finite(FDR_study),FDR_study>=0,FDR_study<fdr_cutoff)
  if(method == "GSEA") {
    x <- x |> filter(is.finite(NES),size>0) |>
      arrange(desc(NES),FDR_study,pathway) |> mutate(draw_order=row_number())
  } else {
    x <- x |> filter(is.finite(fold_enrichment),fold_enrichment>0,overlap>0) |>
      mutate(effect_sign=factor(effect_sign,levels=c("positive","negative"))) |>
      arrange(effect_sign,desc(fold_enrichment),FDR_study,pathway) |>
      group_by(effect_sign) |> mutate(draw_order=row_number()) |> ungroup()
  }
  x
}

plot_step05_all <- function(x, method = c("GSEA", "ORA"), fdr_cutoff = 0.1) {
  method <- match.arg(method)
  original <- x
  x <- step05_all_plot_rows(x,method,fdr_cutoff)
  label <- if(nrow(original)) paste(unique(original$contrast_label),collapse="; ") else "Selected comparison"
  collection <- if(nrow(original)) paste(unique(original$collection),collapse="; ") else "Selected collection"
  subtitle <- paste0(collection,"; ",label,"\nAll study-significant terms (FDR < ",fdr_cutoff,
                     "); observed source models only; no Top-N limit")
  if(!nrow(x)) return(ggplot() +
    annotate("text",x=0,y=0,label="No terms pass study FDR < 0.1 in this observed-only analysis",size=4) +
    labs(title=paste(method,"functional enrichment"),subtitle=subtitle,
         caption="Source lod_imputed models excluded before enrichment. An empty result does not establish biological absence.") +
    theme_void() + theme(plot.caption=element_text(hjust=0)))
  x <- x |> mutate(.row_key=paste(if(method=="ORA") as.character(effect_sign) else "GSEA",pathway,sep="___"),
    .term_label=str_wrap(paste0(term_description," [",term_id,"]"),55),
    .logq=-log10(pmax(FDR_study,.Machine$double.xmin)))
  # Discrete levels run bottom-to-top: reversing the requested descending rows
  # puts the largest signed NES/fold enrichment at the top of every panel.
  x$.row_key <- factor(x$.row_key,levels=rev(unique(x$.row_key)))
  row_labels <- setNames(x$.term_label,as.character(x$.row_key))
  if(method=="GSEA") {
    p <- ggplot(x,aes(NES,.row_key,colour=.logq,size=size)) +
      geom_vline(xintercept=0,colour="grey65",linewidth=.35) + geom_point(alpha=.9) +
      scale_size_continuous(range=c(2,6),name="Eligible pathway\nmembers") +
      labs(x="Normalized enrichment score (signed NES)",
           caption=paste0("Rows: largest positive NES to most negative NES. ",.step05_direction_key(original,105),
                          "\n",str_wrap("Source lod_imputed models excluded before enrichment; no minimum-replicate filter added. Enrichment is not a direct measure of pathway activity.",120)))
  } else {
    x$.panel <- factor(x$direction,levels=unique(x$direction))
    p <- ggplot(x,aes(fold_enrichment,.row_key,colour=.logq,size=overlap)) +
      geom_vline(xintercept=1,colour="grey65",linetype="dashed",linewidth=.35) + geom_point(alpha=.9) +
      facet_wrap(~.panel,ncol=min(2,n_distinct(x$.panel)),scales="free_y",
                 labeller=labeller(.panel=function(v) str_wrap(v,44))) +
      scale_size_continuous(range=c(2,6),name="Overlap\nproteins") +
      labs(x="Fold enrichment relative to the eligible tested background",
           caption=str_wrap("Rows within each direction: largest to smallest fold enrichment. Source lod_imputed models excluded before enrichment; no minimum-replicate filter added. Categories can share proteins; enrichment does not establish direct kinase targets.",120))
  }
  p + scale_y_discrete(labels=row_labels) +
    scale_colour_viridis_c(option="C",name="-log10(study FDR)") +
    labs(title=paste(method,"functional enrichment"),subtitle=subtitle,y=NULL) +
    theme_bw(base_size=11) +
    theme(legend.position="bottom",axis.text.y=element_text(size=9),
      plot.subtitle=element_text(size=11),strip.text=element_text(size=10),
      panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0,size=9),
      plot.margin=margin(10,15,10,10))
}

step05_all_plot_height <- function(rows,method=c("GSEA","ORA")) {
  method<-match.arg(method)
  if(!nrow(rows)) return(4.5)
  line_count<-str_count(str_wrap(paste0(rows$term_description," [",rows$term_id,"]"),55),"\n")+1L
  row_inches<-.13*line_count+.14
  if(method=="ORA") {
    panel_inches<-tapply(row_inches,rows$effect_sign,sum)
    panel_inches<-panel_inches[is.finite(panel_inches)]
    height<-max(panel_inches)+3.2
  } else height<-sum(row_inches)+3.8
  max(5,height)
}
