# Fresh, exact-input STRING 12.5 networks for poster figures (49, 67 and 75 nodes).
# No automatic partner expansion; the original67 is a distinct input list.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/06_string_network/reports/poster_networks")
raw <- file.path(REPO_ROOT,"data/raw/step06_string/poster_api")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(raw,recursive=TRUE,showWarnings=FALSE)
refresh <- tolower(Sys.getenv("JE2_REFRESH_STRING_POSTER",unset="false"))=="true"
base_url <- "https://version-12-5.string-db.org/api"
lists <- list(`49`=readLines(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/user_selected/combined_selected_plus_anchors.txt")),
  `67`=readLines(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/all_selected_plus_Stk1_Stp1.txt")),
  `75`=readLines(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/complete75_plus_anchors.txt")))
stopifnot(length(lists[[1]])==49L,length(lists[[2]])==67L,length(lists[[3]])==75L,
  all(vapply(lists,function(x)!anyNA(x)&&!anyDuplicated(x)&&all(nzchar(x)),logical(1))))
all_ids <- unique(unlist(lists,use.names=FALSE))
stopifnot(length(all_ids)==75L,length(intersect(lists[[1]],lists[[2]]))==42L,
  all(lists[[1]] %in% lists[[3]]),all(lists[[2]] %in% lists[[3]]))
provenance <- list()
call <- function(name,endpoint,params,filename=NULL) {
  r <- .api_request(paste0(base_url,endpoint),params=params,service="STRING",refresh=refresh)
  provenance[[name]] <<- r$provenance
  if(!is.null(filename)) writeLines(r$body,file.path(raw,filename),useBytes=TRUE)
  r
}
# Map in <=50-ID batches. Do not split network requests: that would omit edges
# between batches. Each complete network fits safely in one GET request.
mapped <- list()
chunks <- split(all_ids,ceiling(seq_along(all_ids)/50))
for(i in seq_along(chunks)) {
  message("STRING mapping batch ",i,"/",length(chunks))
  r <- call(paste0("mapping_",i),"/tsv/get_string_ids",list(identifiers=paste(chunks[[i]],collapse="\r"),
    species=93061,limit=1,echo_query=1,caller_identity="JE2-Project-I-poster"),paste0("mapping_",i,".tsv"))
  z <- .api_tsv(r)
  index <- as.integer(z$queryIndex)+1L
  stopifnot(all(index %in% seq_along(chunks[[i]])),all(z$ncbiTaxonId=="93061"))
  mapped[[i]] <- z |> mutate(uniprot_acc=chunks[[i]][index])
}
mapping <- bind_rows(mapped) |> transmute(uniprot_acc,string_id=stringId,string_name=preferredName,
  string_description=annotation,taxon=ncbiTaxonId)
stopifnot(nrow(mapping)==length(all_ids),!anyDuplicated(mapping$uniprot_acc),!anyDuplicated(mapping$string_id),
  setequal(mapping$uniprot_acc,all_ids))

# Verify every service-selected ID against exact UniProt STRING cross-references.
references <- list()
for(i in seq_along(chunks)) {
  r <- .api_request("https://rest.uniprot.org/uniprotkb/search",params=list(
    query=paste0("accession:",chunks[[i]],collapse=" OR "),format="tsv",
    fields="accession,xref_string",size=500),service="UniProt",refresh=refresh)
  provenance[[paste0("uniprot_xref_",i)]] <- r$provenance
  writeLines(r$body,file.path(raw,paste0("uniprot_crossrefs_",i,".tsv")),useBytes=TRUE)
  z <- .api_tsv(r)
  stopifnot(nrow(z)<500L)
  references[[i]] <- z |> transmute(uniprot_acc=Entry,string_id=strsplit(STRING,";",fixed=TRUE)) |>
    unnest(string_id) |> mutate(string_id=trimws(string_id)) |> filter(nzchar(string_id)) |> distinct()
}
xrefs <- bind_rows(references)
stopifnot(nrow(anti_join(mapping,xrefs,by=c("uniprot_acc","string_id")))==0L)
mapping$mapping_status <- "API mapping verified against exact UniProt STRING cross-reference"
write_csv(mapping,file.path(out,"identifier_mapping.csv"))

# Use curated JE2 identities to avoid ineligible alternate reference mappings.
identity_fields <- c("uniprot_acc","protein_key","locus_tag","display_name","product")
identity <- bind_rows(
  read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/user_selected/network_node_attributes.csv"),show_col_types=FALSE) |> select(all_of(identity_fields)),
  read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/upload_manifest.csv"),show_col_types=FALSE) |> select(all_of(identity_fields)),
  read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/complete75_restored_eight.csv"),show_col_types=FALSE) |> select(all_of(identity_fields))) |> distinct()
stopifnot(!anyDuplicated(identity$uniprot_acc),setequal(identity$uniprot_acc,all_ids))
inventory <- read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/01_module_inventory.csv"),show_col_types=FALSE)
source_data <- read_tsv(file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),show_col_types=FALSE)
codes <- tibble(contrast=c("interaction_stk_vs_anc_at_TSB_vs_SN","interaction_stp_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_stk_at_TSB_vs_SN"),code=c("stk_anc","stp_anc","stp_stk"))
effects <- source_data |> inner_join(codes,by="contrast",relationship="many-to-one") |>
  semi_join(identity,by="protein_key") |>
  mutate(halo_effect=if_else(estimate_type=="observed" & coalesce(all_relevant_groups_2of3,FALSE) &
      is.finite(diff),diff,NA_real_),halo_FDR=if_else(is.finite(halo_effect),FDR,NA_real_)) |>
  select(protein_key,code,diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR) |>
  pivot_wider(names_from=code,values_from=c(diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR))
stopifnot(nrow(effects)==length(all_ids),!anyDuplicated(effects$protein_key))
membership <- bind_rows(lapply(names(lists),function(set) {
  identity |> filter(uniprot_acc %in% lists[[set]]) |>
    inner_join(inventory |> select(protein_key,term_id,module,annotation_review,annotation_note),
      by="protein_key",relationship="one-to-many") |>
    left_join(select(mapping,uniprot_acc,string_id),by="uniprot_acc",relationship="many-to-one") |>
    mutate(set_id=set,.before=1)
}))
write_csv(membership,file.path(out,"pathway_membership.csv"))
compare <- tibble(uniprot_acc=all_ids,in_49=all_ids %in% lists[[1]],in_67=all_ids %in% lists[[2]],
    in_75=all_ids %in% lists[[3]]) |>
  left_join(identity,by="uniprot_acc",relationship="one-to-one")
write_csv(compare,file.path(out,"input_set_comparison.csv"))

summary <- annotation_summary <- list()
module_keys <- c(P="Purine metabolism",R="Oxidative phosphorylation",N="Nitrogen metabolism",T="Aminoacyl-tRNA biosynthesis")
module_code <- function(x) paste(names(module_keys)[module_keys %in% x],collapse="/")
for(set in names(lists)) {
  message("Retrieving complete ",set,"-protein network and annotations")
  ids <- lists[[set]]
  this_map <- mapping[match(ids,mapping$uniprot_acc),]
  params <- list(identifiers=paste(this_map$string_id,collapse="\r"),species=93061,
    network_type="functional",required_score=400,add_nodes=0,caller_identity="JE2-Project-I-poster")
  r <- call(paste0("network_",set),"/tsv/network",params,paste0("network_",set,".tsv"))
  z <- .api_tsv(r)
  stopifnot(all(z$stringId_A %in% this_map$string_id),all(z$stringId_B %in% this_map$string_id))
  edges <- z |> transmute(node_a=pmin(stringId_A,stringId_B),node_b=pmax(stringId_A,stringId_B),
    combined_score=as.numeric(score),neighborhood=as.numeric(nscore),fusion=as.numeric(fscore),
    cooccurrence=as.numeric(pscore),coexpression=as.numeric(ascore),experimental=as.numeric(escore),
    database=as.numeric(dscore),textmining=as.numeric(tscore)) |> distinct()
  stopifnot(!anyDuplicated(edges[c("node_a","node_b")]),all(edges$node_a!=edges$node_b),
    all(edges$combined_score>=.4),all(edges$combined_score<=1))
  write_csv(edges,file.path(out,paste0("edges_",set,".csv")))
  this_membership <- filter(membership,set_id==set) |> group_by(protein_key) |>
    summarise(module_group=paste(unique(module),collapse="; "),annotation_review=any(annotation_review),
      pathway_codes=module_code(module),.groups="drop")
  nodes <- this_map |> left_join(identity,by="uniprot_acc",relationship="one-to-one") |>
    left_join(effects,by="protein_key",relationship="one-to-one") |>
    left_join(this_membership,by="protein_key",relationship="one-to-one") |>
    mutate(node_role=if_else(uniprot_acc %in% c("Q2FZ64","Q2FZ65"),"Biological anchor","Selected protein"),
      module_group=if_else(node_role=="Biological anchor","Regulatory anchors",coalesce(module_group,"Other or unassigned")),
      annotation_review=coalesce(annotation_review,FALSE),
      degree=as.integer(table(factor(c(edges$node_a,edges$node_b),levels=string_id))))
  stopifnot(nrow(nodes)==as.integer(set),!anyNA(nodes$protein_key),!anyDuplicated(nodes$string_id))
  write_csv(nodes,file.path(out,paste0("nodes_",set,".csv")))
  write_csv(this_map,file.path(out,paste0("identifier_mapping_",set,".csv")))
  write_csv(nodes |> select(string_id,uniprot_acc,display_name,degree),file.path(out,paste0("degrees_",set,".csv")))
  annotation <- call(paste0("annotation_",set),"/tsv/functional_annotation",
    list(identifiers=paste(this_map$string_id,collapse="\r"),species=93061,
      caller_identity="JE2-Project-I-poster"),paste0("functional_annotations_",set,".tsv"))
  az <- .api_tsv(annotation)
  # Functional annotations aggregate genes per term; they are not enrichment tests.
  annotation_genes <- strsplit(az$inputGenes,",",fixed=TRUE)
  annotated_ids <- unique(trimws(unlist(annotation_genes)))
  stopifnot(all(annotated_ids %in% this_map$string_id),all(az$ncbiTaxonId=="93061"),
    all(lengths(annotation_genes)==as.integer(az$number_of_genes)))
  annotation_summary[[set]] <- tibble(set_id=set,annotation_rows=nrow(az),annotated_nodes=length(annotated_ids))
  # Retrieve reference FASTA independently; never use it as JE2 phosphosite coordinates.
  fasta <- .api_request("https://rest.uniprot.org/uniprotkb/stream",params=list(
    query=paste0("accession:",ids,collapse=" OR "),format="fasta"),service="UniProt",refresh=refresh)
  provenance[[paste0("reference_fasta_",set)]] <- fasta$provenance
  writeLines(fasta$body,file.path(raw,paste0("reference_sequences_",set,".fa")),useBytes=TRUE)
  fasta_ids <- str_match(strsplit(fasta$body,"\n",fixed=TRUE)[[1]],"^>[^|]+[|]([^|]+)[|]")[,2]
  stopifnot(setequal(na.omit(fasta_ids),ids),sum(!is.na(fasta_ids))==length(ids))
  link <- call(paste0("network_link_",set),"/tsv/get_link",params,paste0("STRING_network_link_",set,".tsv"))
  count_nodes <- nrow(nodes); count_edges <- nrow(edges)
  summary[[set]] <- tibble(set_id=set,isolates=sum(nodes$degree==0),
    edges_ge_07=sum(edges$combined_score>=.7),membership_flags=sum(nodes$annotation_review),
    automatic_added_nodes=0L,nodes=count_nodes,edges=count_edges) |>
    select(set_id,nodes,edges,everything())
}
write_csv(bind_rows(summary),file.path(out,"network_summary.csv"))
write_csv(bind_rows(annotation_summary),file.path(out,"annotation_export_summary.csv"))
write_csv(bind_rows(lapply(names(provenance),function(name) {
  p <- provenance[[name]]
  tibble(query=name,url=p$request_url,retrieved_utc=p$retrieved_utc,status=p$http_status,cache_hit=p$cache_hit)
})),file.path(out,"API_provenance.csv"))
raw_files <- list.files(raw,full.names=TRUE)
raw_files <- raw_files[!grepl("[.]zip$",raw_files)]
write_csv(tibble(file=basename(raw_files),bytes=file.info(raw_files)$size,
  md5=unname(tools::md5sum(raw_files))),file.path(out,"raw_export_inventory.csv"))
write_csv(tibble(key=c("STRING_version","species","network_type","required_score","add_nodes",
  "input_sets","shared_inputs","automatic_partner_expansion","central_anchors"),
  value=c("12.5","93061","functional","400","0","exact49; exact67; exact75","42 (49 vs 67)","none","Stk1; Stp1")),
  file.path(out,"settings.csv"))
print(bind_rows(summary))
message("STRING poster inputs ready; raw network/annotation/FASTA exports stay under ignored data/raw/.")
