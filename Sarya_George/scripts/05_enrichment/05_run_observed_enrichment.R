# Explicit source-observed-only enrichment; no extra replication-count restriction.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(!requireNamespace("fgsea",quietly=TRUE)) stop("Install fgsea with scripts/00_setup/install_enrichment_packages.R.")
out<-file.path(REPO_ROOT,"results","step05","observed_only")
report<-file.path(REPO_ROOT,"reports","step05","observed_only")
for(path in c(out,report,file.path(out,"cache"))) dir.create(path,recursive=TRUE,showWarnings=FALSE)
source_file<-file.path(REPO_ROOT,"results","step04","total","annotated_all.tsv.gz")
terms_file<-file.path(REPO_ROOT,"results","step03","terms_cds.tsv.gz")
x<-read_tsv(source_file,show_col_types=FALSE,progress=FALSE)
terms<-read_tsv(terms_file,show_col_types=FALSE,progress=FALSE)
plan<-step04_contrast_plan()
settings<-list(policy="observed_only",min_size=10L,max_size=500L,protein_FDR=.1,protein_abs_log2FC=1,
  pathway_FDR=.1,seed_base=20261113L,eps=1e-10,sample_size=101L,n_perm_simple=1000L,
  family="BH separately by method across ten contrasts and all five original collections; ORA includes both directions")
inputs<-c(source_file,terms_file,file.path(REPO_ROOT,"R","analysis",c("step05_enrichment.R","step05_gsea.R")))
hashes<-tibble(file=substring(inputs,nchar(REPO_ROOT)+2),md5=unname(tools::md5sum(inputs)))
write_tsv(hashes,file.path(report,"input_hashes.tsv"))
ora_rows<-gsea_rows<-audit_rows<-logs<-list()
for(i in seq_len(nrow(plan))) {
  contrast<-plan$contrast[i]; message("[",i,"/10] observed_only : ",contrast)
  prepared<-step05_prepare(filter(x,.data$contrast==.env$contrast),"observed_only")
  stopifnot(all(prepared$estimate_type=="observed"),!anyDuplicated(prepared$target_id))
  dictionary<-step05_pathways(terms,prepared)
  ranks<-setNames(prepared$statistic,prepared$target_id)
  ranks<-ranks[order(-ranks,names(ranks),method="radix")]
  metadata<-plan[i,] |> mutate(policy="observed_only",namespace="UniProt")
  seed<-settings$seed_base+i
  fingerprint<-digest::digest(list(hashes,settings,seed,ranks,dictionary$membership,
    fgsea=as.character(packageVersion("fgsea"))),algo="sha256")
  cache_file<-file.path(out,"cache",paste0(contrast,".rds"))
  cache<-if(file.exists(cache_file)) readRDS(cache_file) else NULL
  if(is.null(cache)||!identical(cache$fingerprint,fingerprint)) {
    g<-run_step05_gsea(ranks,dictionary$sets,seed=seed,min_size=settings$min_size,max_size=settings$max_size,
      eps=settings$eps,sample_size=settings$sample_size,n_perm_simple=settings$n_perm_simple)
    cache<-list(fingerprint=fingerprint,gsea=g,ranks=ranks,pathways=dictionary$sets,
      metadata=metadata,seed=seed,prepared=prepared)
    saveRDS(cache,cache_file)
  } else message("  Reusing matching deterministic cache")
  o<-bind_rows(lapply(c("positive","negative"),function(sign) {
    fg<-prepared |> filter(FDR<settings$protein_FDR,
      if(sign=="positive") diff>settings$protein_abs_log2FC else diff< -settings$protein_abs_log2FC) |> pull(target_id)
    run_step05_ora(fg,prepared$target_id,dictionary$sets,min_size=settings$min_size,max_size=settings$max_size) |>
      mutate(effect_sign=sign)
  })) |> left_join(dictionary$metadata,by="pathway",relationship="many-to-one")
  for(field in names(metadata)) o[[field]]<-metadata[[field]][1]
  o<-o |> mutate(direction=if_else(effect_sign=="positive",positive_direction,negative_direction),
    overlap_names=vapply(overlap_ids,step05_member_labels,character(1),mapping=prepared))
  g<-cache$gsea |> mutate(leading_edge_ids=vapply(leadingEdge,function(v) paste(v,collapse=";"),character(1))) |>
    select(-leadingEdge) |> left_join(dictionary$metadata,by="pathway",relationship="many-to-one")
  for(field in names(metadata)) g[[field]]<-metadata[[field]][1]
  g<-g |> mutate(effect_sign=case_when(NES>0~"positive",NES<0~"negative",TRUE~"unavailable"),
    direction=case_when(NES>0~positive_direction,NES<0~negative_direction,TRUE~"Unavailable"),
    leading_edge_names=vapply(leading_edge_ids,step05_member_labels,character(1),mapping=prepared),
    background_size=nrow(prepared),seed=seed)
  ora_rows[[i]]<-o; gsea_rows[[i]]<-g
  audit_rows[[i]]<-tibble(contrast,universe=nrow(prepared),source_observed=sum(prepared$estimate_type=="observed"),
    source_LOD=sum(prepared$estimate_type=="lod_imputed"),sparse_observed_retained=sum(!prepared$all_relevant_groups_2of3),
    positive_hits=sum(prepared$FDR<.1 & prepared$diff>1),negative_hits=sum(prepared$FDR<.1 & prepared$diff< -1),
    tested_terms=nrow(g),failed_GSEA=sum(!is.finite(g$pval)))
  warnings<-attr(cache$gsea,"gsea_warnings")
  logs[[i]]<-tibble(contrast,seed,warning=if(length(warnings)) paste(gsub("[\r\n]+"," ",warnings),collapse=" | ") else "none")
  write_tsv(prepared |> select(target_id,protein_key,locus_tag,display_name,estimate_type,statistic,diff,FDR,
    all_relevant_groups_2of3,group_detection_pattern),file.path(out,paste0(contrast,"__rank_audit.tsv.gz")))
}
ora<-step05_adjust(bind_rows(ora_rows),"p_value") |> arrange(contrast,collection,effect_sign,desc(fold_enrichment),pathway)
gsea<-step05_adjust(bind_rows(gsea_rows),"pval") |> arrange(contrast,collection,desc(NES),pathway)
write_tsv(ora,file.path(out,"ora_all.tsv.gz")); write_tsv(gsea,file.path(out,"gsea_all.tsv.gz"))
for(method in c("ORA","GSEA")) {
  table<-if(method=="ORA") ora else gsea
  significant<-filter(table,is.finite(FDR_study),FDR_study<.1,collection!="UniProt_keywords")
  write_tsv(significant,file.path(report,paste0(tolower(method),"_all_significant_GO_KEGG.tsv")))
}
counts<-bind_rows(list(ORA=ora,GSEA=gsea),.id="method") |> group_by(method,contrast,collection) |>
  summarise(tested=n(),significant=sum(FDR_study<.1,na.rm=TRUE),.groups="drop")
write_tsv(counts,file.path(report,"counts.tsv"))
write_tsv(bind_rows(audit_rows),file.path(report,"universe_audit.tsv"))
write_tsv(bind_rows(logs),file.path(report,"gsea_run_log.tsv"))
writeLines(c("Source October total proteome; observed_only retains estimate_type=='observed'.",
  "No additional >=2/3 measurement filter. Some retained proteins have only1/3 measured values in a relevant group.",
  "UniProt-quality and singleton mapping filters are unchanged. Background, GSEA ranks and ORA lists all exclude source-LOD flags.",
  "All five original annotation collections are tested and included in BH correction, separately for each method across ten contrasts.",
  "GO BP/MF/CC and KEGG displayed without top-N truncation; significance is study-wide FDR<0.1.",
  "Existing 10-500 effective-member size eligibility is unchanged; source protein effects and FDRs are not refitted.",
  paste(names(settings),vapply(settings,as.character,character(1)),sep=" = ")),
  file.path(report,"settings.txt"))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(report,"sessionInfo.txt"))
message("OBSERVED ONLY COMPLETE: ",nrow(ora)," ORA tests, ",nrow(gsea)," GSEA tests.")
