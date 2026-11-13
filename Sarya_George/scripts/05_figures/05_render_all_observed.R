# Every significant GO/KEGG result in the observed-only analysis, sorted by effect.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

render_all_observed <- function() {
  result_dir<-file.path(REPO_ROOT,"results","step05","observed_only")
  report_dir<-file.path(REPO_ROOT,"reports","step05","observed_only")
  figure_dir<-file.path(REPO_ROOT,"figures","step05","observed_only")
  dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
  paths<-file.path(result_dir,c("gsea_all.tsv.gz","ora_all.tsv.gz"))
  if(any(!file.exists(paths))) stop("Run the observed-only enrichment analysis before rendering all its results.")
  data<-list(GSEA=read_tsv(paths[1],show_col_types=FALSE),ORA=read_tsv(paths[2],show_col_types=FALSE))
  plan<-step04_contrast_plan()
  collections<-c("KEGG","GO_BP","GO_MF","GO_CC")
  manifest<-list(); row_tables<-list(); index<-0L
  for(method in names(data)) for(collection_name in collections) for(contrast_name in plan$contrast) {
    index<-index+1L
    x<-filter(data[[method]],collection==collection_name,contrast==contrast_name)
    rows<-step05_all_plot_rows(x,method)
    plot<-plot_step05_all(x,method)
    height<-step05_all_plot_height(rows,method)
    width<-if(method=="ORA" && n_distinct(rows$direction)>1) 20 else 14.5
    stem<-paste0("step05_",collection_name,"__",contrast_name,"__",method,"__all_observed")
    png<-file.path(figure_dir,paste0(stem,".png")); pdf<-file.path(figure_dir,paste0(stem,".pdf"))
    # Tall complete plots are kept intact. Lower export DPI only when needed to
    # keep the raster dimension reasonable; the PDF remains resolution independent.
    dpi<-if(height>65) 180L else if(height>45) 220L else 300L
    message("[",index,"/80] ",stem," : ",nrow(rows)," significant rows")
    ggsave(png,plot,width=width,height=height,units="in",dpi=dpi,bg="white",limitsize=FALSE)
    ggsave(pdf,plot,width=width,height=height,units="in",device=grDevices::pdf,useDingbats=FALSE,bg="white",limitsize=FALSE)
    # Store exactly the plotted ordering; ORA ranks restart within direction.
    ordered<-rows |> mutate(method=method,.before=1)
    write_tsv(ordered,file.path(report_dir,paste0(stem,".tsv")))
    row_tables[[index]]<-ordered
    manifest[[index]]<-tibble(method,collection=collection_name,contrast=contrast_name,
      significant_rows=nrow(rows),positive_rows=sum(rows$direction==if(nrow(x)) x$positive_direction[1] else "",na.rm=TRUE),
      negative_rows=sum(rows$direction==if(nrow(x)) x$negative_direction[1] else "",na.rm=TRUE),
      policy="observed_only",threshold="FDR_study < 0.1",top_limit="none",
      row_order=if(method=="GSEA") "descending signed NES" else "descending fold enrichment within direction",
      width_inches=width,height_inches=height,png_dpi=dpi,
      png=paste0("figures/step05/observed_only/",stem,".png"),
      pdf=paste0("figures/step05/observed_only/",stem,".pdf"),
      table=paste0("reports/step05/observed_only/",stem,".tsv"))
  }
  manifest<-bind_rows(manifest)
  write_tsv(manifest,file.path(report_dir,"figure_manifest.tsv"))
  write_tsv(bind_rows(row_tables),file.path(report_dir,"all_plotted_significant.tsv"))
  heading<-c("# All significant observed-only GO/KEGG results","",
    "This folder is the complete observed-only presentation requested after the original Step 5 figures. Earlier Top-N plots in figures/step05/ remain an explicitly different primary-analysis view; use this folder for the no-LOD, all-significant results.","",
    "These figures exclude source-model lod_imputed proteins before enrichment. They use the observed-only analysis, without adding a minimum-replicate filter. Missing sample measurements remain missing.","",
    "All terms passing study-wide BH FDR <0.1 are displayed. There is no Top-N cap. The enrichment correction family is defined by the observed-only analysis, before selecting the four displayed GO/KEGG collections.","",
    "GSEA rows run from the largest positive NES at the top to the most negative NES at the bottom. ORA rows run from highest to lowest fold enrichment within each directional panel. Empty plots explicitly report that no terms pass the threshold.","",
    "GO_BP means Biological Process; GO_MF means Molecular Function; GO_CC means Cellular Component. KEGG supplies pathway maps. GO/KEGG categories are the membership definitions; ORA and GSEA are the statistical methods applied to those memberships.","",
    "Definitions: [Gene Ontology documentation](https://geneontology.org/docs/ontology-documentation/) and [KEGG PATHWAY](https://www.kegg.jp/kegg/pathway.html). These official resource definitions do not refresh the frozen course membership snapshot.","",
    "Each comparison/collection has a separate complete PNG and PDF, with an ordered TSV containing exactly its displayed rows. Long plots have a proportionally tall canvas; use the PDF to zoom. PNG resolution is recorded in the manifest.","",
    "[All ordered significant rows](all_plotted_significant.tsv) | [Figure manifest](figure_manifest.tsv)","",
    "| Method | Collection | Comparison | Significant rows | PNG | PDF | Ordered table |",
    "| --- | --- | --- | ---: | --- | --- | --- |")
  entries<-vapply(seq_len(nrow(manifest)),function(i) {
    row<-manifest[i,]
    paste0("| ",row$method," | ",row$collection," | ",row$contrast," | ",row$significant_rows,
      " | [PNG](../../../",row$png,") | [PDF](../../../",row$pdf,") | [TSV](",basename(row$table),") |")
  },character(1))
  writeLines(c(heading,entries),file.path(report_dir,"README.md"))
  writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(report_dir,"figure_sessionInfo.txt"))
  message("All-results plots complete: ",nrow(manifest)," PNG/PDF pairs; no significant rows omitted.")
  invisible(manifest)
}

render_all_observed()
