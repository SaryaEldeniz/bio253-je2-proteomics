# Step 5: declared universes, upper-tail ORA and multiple-testing bookkeeping.
step05_collections <- function() {
  tribble(~category,~collection,
    "KEGG (Kyoto Encyclopedia of Genes and Genomes)","KEGG",
    "Biological Process (Gene Ontology)","GO_BP",
    "Molecular Function (Gene Ontology)","GO_MF",
    "Cellular Component (Gene Ontology)","GO_CC",
    "Annotated Keywords (UniProt)","UniProt_keywords")
}

step05_prepare <- function(x,policy="primary") {
  policy<-match.arg(policy,c("primary","replicated","observed_replicated","local_JE2"))
  if(anyDuplicated(x$protein_key)) stop("Prepare one contrast at a time; protein keys must be unique.")
  x<-filter(x,eligible_cds,is.finite(diff),is.finite(FDR),is.finite(statistic))
  if(policy!="local_JE2") {
    x<-filter(x,eligible_uniprot_annotation,n_JE2_good_total_DEA==1,!is.na(uniprot_acc))
  }
  if(policy %in% c("replicated","observed_replicated")) x<-filter(x,all_relevant_groups_2of3)
  if(policy=="observed_replicated") x<-filter(x,estimate_type=="observed")
  x<-mutate(x,target_id=if(policy=="local_JE2") protein_key else uniprot_acc)
  if(anyDuplicated(x$target_id)) stop("Target collisions must be resolved explicitly before enrichment.")
  x
}

step05_pathways <- function(terms,mapping) {
  dictionary<-terms |> inner_join(step05_collections(),by="category") |>
    transmute(protein_key,collection,term_id=term,term_description) |>
    inner_join(select(mapping,protein_key,target_id),by="protein_key",relationship="many-to-one") |>
    mutate(pathway=paste(collection,term_id,sep="::")) |>
    distinct(pathway,collection,term_id,term_description,target_id)
  metadata<-dictionary |> distinct(pathway,collection,term_id,term_description)
  if(anyDuplicated(metadata$pathway)) stop("Conflicting descriptions for the same pathway key.")
  list(sets=split(dictionary$target_id,dictionary$pathway),metadata=metadata,membership=dictionary)
}

run_step05_ora <- function(foreground,universe,pathways,min_size=10L,max_size=500L) {
  universe<-unique(universe); foreground<-unique(foreground)
  if(anyNA(universe)||anyNA(foreground)||!all(foreground %in% universe)) stop("Foreground must be a subset of a nonmissing universe.")
  sets<-lapply(pathways,function(s) intersect(unique(s),universe))
  sizes<-lengths(sets)
  sets<-sets[sizes>=min_size & sizes<=max_size & sizes<length(universe)]
  N<-length(universe); n<-length(foreground)
  if(!length(sets)) return(tibble(pathway=character(),overlap=integer(),foreground_size=integer(),term_size=integer(),
                                  background_size=integer(),expected_overlap=double(),fold_enrichment=double(),p_value=double(),overlap_ids=character()))
  bind_rows(lapply(names(sets),function(name) {
    members<-sets[[name]]; hit<-intersect(foreground,members); k<-length(hit); K<-length(members)
    tibble(pathway=name,overlap=k,foreground_size=n,term_size=K,background_size=N,
      expected_overlap=n*K/N,fold_enrichment=if(n>0) (k/n)/(K/N) else NA_real_,
      p_value=if(n>0) phyper(k-1,K,N-K,n,lower.tail=FALSE) else 1,
      overlap_ids=paste(sort(hit),collapse=";"))
  }))
}

step05_adjust <- function(x,p_column) {
  # Failed GSEA tests stay NA, but remain in the declared BH family size.
  bh<-function(p) { out<-rep(NA_real_,length(p)); ok<-is.finite(p); out[ok]<-p.adjust(p[ok],"BH",n=length(p)); out }
  x |> group_by(policy,contrast) |> mutate(FDR_local=bh(.data[[p_column]])) |> ungroup() |>
    group_by(policy) |> mutate(FDR_study=bh(.data[[p_column]])) |> ungroup()
}

step05_member_labels <- function(ids,mapping) {
  if(is.na(ids)||!nzchar(ids)) return(NA_character_)
  keys<-strsplit(ids,";",fixed=TRUE)[[1]]
  paste(mapping$display_name[match(keys,mapping$target_id)],collapse=";")
}

read_enrichment <- function(method="ORA",contrast=NULL,policy="primary",significant_only=TRUE,collection=NULL) {
  method<-match.arg(toupper(method),c("ORA","GSEA"))
  path<-file.path(REPO_ROOT,"results","step05",paste0(tolower(method),"_all.tsv.gz"))
  if(!file.exists(path)) stop('First run run_project_step("enrichment").')
  x<-read_tsv(path,show_col_types=FALSE) |> filter(.data$policy==.env$policy)
  if(!is.null(contrast)) x<-filter(x,.data$contrast %in% .env$contrast)
  if(!is.null(collection)) x<-filter(x,.data$collection %in% .env$collection)
  if(significant_only) x<-filter(x,is.finite(FDR_study),FDR_study<.1)
  arrange(x,FDR_study)
}
