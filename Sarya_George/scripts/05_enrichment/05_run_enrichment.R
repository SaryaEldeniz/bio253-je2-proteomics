# October total proteome: pathway discovery with declared backgrounds and sensitivities.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(!requireNamespace("fgsea",quietly=TRUE)) stop("Run scripts/00_setup/install_enrichment_packages.R first.")
out<-file.path(REPO_ROOT,"results","step05")
published<-file.path(REPO_ROOT,"reports","step05")
for(folder in c(out,file.path(out,"cache"),file.path(out,"contrasts"),published)) dir.create(folder,recursive=TRUE,showWarnings=FALSE)
source_file<-file.path(REPO_ROOT,"results","step04","total","annotated_all.tsv.gz")
terms_file<-file.path(REPO_ROOT,"results","step03","terms_cds.tsv.gz")
map_file<-file.path(REPO_ROOT,"results","step03","id_map.tsv")
stopifnot(file.exists(source_file),file.exists(terms_file),file.exists(map_file))
source_data<-read_tsv(source_file,show_col_types=FALSE,progress=FALSE)
terms<-read_tsv(terms_file,show_col_types=FALSE,progress=FALSE)
plan<-step04_contrast_plan()
policies<-c("primary","replicated","observed_replicated","local_JE2")
settings<-list(min_size=10L,max_size=500L,protein_FDR=.1,protein_abs_log2FC=1,
  pathway_FDR=.1,seed_base=20261112L,eps=1e-10,sample_size=101L,n_perm_simple=1000L,
  score="signed source DEA statistic; no jitter; deterministic ID order at ties",
  family="BH by method and policy across all 10 contrasts and 5 collections; ORA includes both directions")
hash_paths<-c(source_file,terms_file,map_file,file.path(REPO_ROOT,"R","analysis",c("step05_enrichment.R","step05_gsea.R")))
manifest<-tibble(file=substring(hash_paths,nchar(REPO_ROOT)+2),md5=unname(tools::md5sum(hash_paths)))
write_tsv(manifest,file.path(published,"input_hashes.tsv"))
write_tsv(plan,file.path(published,"comparison_plan.tsv"))
saveRDS(settings,file.path(out,"settings.rds"))
all_ora<-all_gsea<-audits<-coverage<-warnings_log<-list()
audit_exclusions<-list()
run_index<-0L
for(policy in policies) for(contrast in plan$contrast) {
  run_index<-run_index+1L
  message("[",run_index,"/40] ",policy," : ",contrast)
  x<-filter(source_data,.data$contrast==.env$contrast)
  prepared<-step05_prepare(x,policy)
  dictionary<-step05_pathways(terms,prepared)
  rank_data<-setNames(prepared$statistic,prepared$target_id)
  positive<-filter(prepared,FDR<settings$protein_FDR,diff>settings$protein_abs_log2FC)$target_id
  negative<-filter(prepared,FDR<settings$protein_FDR,diff< -settings$protein_abs_log2FC)$target_id
  metadata<-filter(plan,.data$contrast==.env$contrast) |> mutate(policy=policy,
    namespace=if_else(policy=="local_JE2","JE2","UniProt"))
  seed<-settings$seed_base+run_index
  fingerprint<-digest::digest(list(manifest,settings,policy,contrast,prepared$target_id,rank_data,
    dictionary$membership,fgsea_version=as.character(packageVersion("fgsea"))),algo="sha256")
  cache_file<-file.path(out,"cache",paste0(policy,"__",contrast,".rds"))
  cache<-if(file.exists(cache_file)) readRDS(cache_file) else NULL
  if(is.null(cache)||!identical(cache$fingerprint,fingerprint)||!identical(cache$seed,seed)) {
    gsea<-run_step05_gsea(rank_data,dictionary$sets,seed=seed,min_size=settings$min_size,
      max_size=settings$max_size,eps=settings$eps,sample_size=settings$sample_size,
      n_perm_simple=settings$n_perm_simple)
    cache<-list(fingerprint=fingerprint,ranks=rank_data,pathways=dictionary$sets,
      metadata=metadata,gsea=gsea,seed=seed)
    saveRDS(cache,cache_file)
  } else message("  Reusing matching deterministic cache")
  warning_text<-attr(cache$gsea,"gsea_warnings")
  warnings_log[[run_index]]<-tibble(policy,contrast,seed,
    warning=if(length(warning_text)) paste(gsub("[\r\n]+"," ",warning_text),collapse=" | ") else "none",
    elapsed_seconds=attr(cache$gsea,"elapsed_seconds"))
  ora<-bind_rows(lapply(c("positive","negative"),function(sign) {
    foreground<-if(sign=="positive") positive else negative
    run_step05_ora(foreground,prepared$target_id,dictionary$sets,
      min_size=settings$min_size,max_size=settings$max_size) |> mutate(effect_sign=sign)
  }))
  ora<-ora |> left_join(dictionary$metadata,by="pathway",relationship="many-to-one") |>
    mutate(policy=policy,contrast=contrast) |> left_join(metadata,by=c("policy","contrast"),relationship="many-to-one") |>
    mutate(direction=if_else(effect_sign=="positive",positive_direction,negative_direction),
      overlap_names=vapply(overlap_ids,step05_member_labels,character(1),mapping=prepared))
  gsea<-cache$gsea |> mutate(leading_edge_ids=vapply(leadingEdge,function(v) paste(v,collapse=";"),character(1))) |>
    select(-leadingEdge) |> left_join(dictionary$metadata,by="pathway",relationship="many-to-one") |>
    mutate(policy=policy,contrast=contrast) |> left_join(metadata,by=c("policy","contrast"),relationship="many-to-one") |>
    mutate(effect_sign=case_when(NES>0~"positive",NES<0~"negative",TRUE~"unavailable"),
      direction=case_when(NES>0~positive_direction,NES<0~negative_direction,TRUE~"Unavailable"),
      leading_edge_names=vapply(leading_edge_ids,step05_member_labels,character(1),mapping=prepared),
      background_size=nrow(prepared),seed=seed)
  all_ora[[run_index]]<-ora; all_gsea[[run_index]]<-gsea
  sizes<-dictionary$membership |> count(pathway,collection,name="effective_size") |>
    mutate(tested=effective_size>=settings$min_size & effective_size<=settings$max_size & effective_size<nrow(prepared))
  coverage[[run_index]]<-step05_collections() |> select(collection) |>
    left_join(dictionary$membership |> group_by(collection) |> summarise(annotated_proteins=n_distinct(target_id),.groups="drop"),by="collection") |>
    left_join(sizes |> group_by(collection) |> summarise(terms_present=n(),terms_tested=sum(tested),.groups="drop"),by="collection") |>
    mutate(policy=policy,contrast=contrast,background_size=nrow(prepared),
      annotated_proteins=coalesce(annotated_proteins,0L),coverage_fraction=annotated_proteins/background_size)
  audits[[run_index]]<-tibble(policy,contrast,source_rows=nrow(x),ordinary_CDS=sum(x$eligible_cds),
    universe=nrow(prepared),positive_hits=length(positive),negative_hits=length(negative),
    source_observed=sum(prepared$estimate_type=="observed"),source_LOD=sum(prepared$estimate_type=="lod_imputed"),
    relevant_groups_2of3=sum(prepared$all_relevant_groups_2of3),rank_ties=sum(duplicated(prepared$statistic)),
    tested_terms=nrow(gsea),failed_GSEA=sum(!is.finite(gsea$pval)),
    floor_GSEA=sum(gsea$pval_at_floor,na.rm=TRUE))
  audit_exclusions[[run_index]]<-x |> select(protein_key,uniprot_acc,entry_class,
    eligible_uniprot_annotation,n_JE2_good_total_DEA,all_relevant_groups_2of3,estimate_type) |>
    mutate(policy=policy,contrast=contrast,in_universe=protein_key %in% prepared$protein_key)
  write_tsv(prepared |> select(target_id,protein_key,locus_tag,display_name,uniprot_acc,statistic,diff,FDR,
    estimate_type,all_relevant_groups_2of3,group_detection_pattern),
    file.path(out,"contrasts",paste0(policy,"__",contrast,"__rank_audit.tsv.gz")))
}
ora<-step05_adjust(bind_rows(all_ora),"p_value") |> arrange(policy,contrast,FDR_study,p_value)
gsea<-step05_adjust(bind_rows(all_gsea),"pval") |> arrange(policy,contrast,FDR_study,pval)
write_tsv(ora,file.path(out,"ora_all.tsv.gz")); write_tsv(gsea,file.path(out,"gsea_all.tsv.gz"))
audit<-bind_rows(audits); coverage<-bind_rows(coverage)
write_tsv(audit,file.path(published,"universe_and_hit_counts.tsv"))
write_tsv(coverage,file.path(published,"annotation_coverage.tsv"))
write_tsv(bind_rows(warnings_log),file.path(published,"gsea_run_log.tsv"))
write_tsv(bind_rows(audit_exclusions),file.path(out,"eligibility_audit.tsv.gz"))
for(method in c("ORA","GSEA")) {
  table<-if(method=="ORA") ora else gsea
  significant<-filter(table,is.finite(FDR_study),FDR_study<settings$pathway_FDR)
  write_tsv(significant,file.path(out,paste0(tolower(method),"_significant_all_policies.tsv")))
  write_tsv(filter(significant,policy=="primary"),file.path(published,paste0(tolower(method),"_primary_significant.tsv")))
  for(contrast in plan$contrast) write_tsv(filter(table,policy=="primary",.data$contrast==.env$contrast),
    file.path(out,"contrasts",paste0("primary__",contrast,"__",tolower(method),".tsv")))
}
summary_counts<-bind_rows(list(ORA=ora,GSEA=gsea),.id="method") |>
  group_by(method,policy,contrast,contrast_label,collection) |>
  summarise(tests=n(),finite_p=sum(is.finite(if(first(method)=="ORA") p_value else pval)),
    significant_study=sum(FDR_study<.1,na.rm=TRUE),significant_local=sum(FDR_local<.1,na.rm=TRUE),.groups="drop")
write_tsv(summary_counts,file.path(published,"pathway_counts.tsv"))
# Side-by-side sensitivity values retain missing tests rather than call them negative.
ora_sens<-ora |> select(contrast,pathway,effect_sign,policy,FDR_study,fold_enrichment,overlap) |>
  pivot_wider(names_from=policy,values_from=c(FDR_study,fold_enrichment,overlap))
gsea_sens<-gsea |> select(contrast,pathway,policy,FDR_study,NES,status) |>
  pivot_wider(names_from=policy,values_from=c(FDR_study,NES,status))
write_tsv(ora_sens,file.path(out,"ora_sensitivity.tsv.gz")); write_tsv(gsea_sens,file.path(out,"gsea_sensitivity.tsv.gz"))
write_tsv(gsea_sens |> filter(FDR_study_primary<.1),file.path(published,"gsea_primary_sensitivity.tsv"))
write_tsv(ora_sens |> filter(FDR_study_primary<.1),file.path(published,"ora_primary_sensitivity.tsv"))
# Long contributing-member tables are for follow-up, not additional independent evidence.
members<-bind_rows(filter(ora,policy=="primary",FDR_study<.1) |>
  transmute(method="ORA",contrast,pathway,collection,term_description,effect_sign,FDR_study,ids=overlap_ids),
  filter(gsea,policy=="primary",FDR_study<.1) |>
    transmute(method="GSEA",contrast,pathway,collection,term_description,effect_sign,FDR_study,ids=leading_edge_ids)) |>
  separate_rows(ids,sep=";") |> filter(!is.na(ids),nzchar(ids))
member_map<-step05_prepare(filter(source_data,contrast==plan$contrast[1]),"primary") |>
  select(ids=target_id,protein_key,locus_tag,display_name,product,uniprot_acc)
members<-members |> left_join(member_map,by="ids",relationship="many-to-one")
write_tsv(members,file.path(out,"contributing_proteins.tsv.gz"))
write_tsv(filter(members,collection=="KEGG"),file.path(published,"kegg_contributing_proteins.tsv"))
# Targeted module/contribution table: no inferred direct phosphorylation or causality.
drivers<-filter(members,collection %in% c("KEGG","GO_BP")) |>
  group_by(contrast,ids,protein_key,locus_tag,display_name,product,uniprot_acc) |>
  summarise(n_terms=n_distinct(pathway),methods=paste(sort(unique(method)),collapse=";"),
    min_pathway_FDR=min(FDR_study),terms=paste(sort(unique(term_description)),collapse=";"),.groups="drop") |>
  left_join(source_data |> select(contrast,protein_key,diff,FDR,estimate_type,n_observed_min_relevant,group_detection_pattern),
    by=c("contrast","protein_key"),relationship="many-to-one")
write_tsv(drivers,file.path(published,"candidate_module_contributors.tsv"))
capture.output(sessionInfo(),file=file.path(published,"sessionInfo.txt"))
writeLines(c("October source export; total proteome only. Phosphosite enrichment is deferred to the localization/parent-protein workflow.",
  paste(names(settings),vapply(settings,as.character,character(1)),sep=" = "),
  "Memberships: supplied course JE2/STRING snapshot transferred through singleton high-quality UniProt ortholog mappings.",
  "No live API membership is mixed into the enrichment snapshot. No source protein-level effects or FDRs are refitted.",
  "Significant means pathway study-wide BH FDR < 0.1; ORA protein thresholds source FDR < 0.1 and abs(log2 effect) > 1."),
  file.path(published,"run_settings.txt"))
message("Completed ",nrow(ora)," ORA tests and ",nrow(gsea)," GSEA tests. Outputs: ",out)
