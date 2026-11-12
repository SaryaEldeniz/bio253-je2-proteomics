# Small, dated public-ID lookups for interpretation; never changes enrichment memberships.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out<-file.path(REPO_ROOT,"reports","step05")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
gsea<-read_enrichment("GSEA",collection="KEGG")
focus<-c("TSB_vs_SN_at_anc",step04_contrast_plan() |> filter(contrast_type=="interaction") |> pull(contrast))
selected<-gsea |> filter(contrast %in% focus) |> group_by(contrast) |>
  arrange(FDR_study,pathway,.by_group=TRUE) |> slice_head(n=2) |> ungroup()
map<-read_identifier_map()
leading<-selected |> select(contrast,leading_edge_ids) |> separate_rows(leading_edge_ids,sep=";") |>
  filter(nzchar(leading_edge_ids)) |> count(uniprot_acc=leading_edge_ids,sort=TRUE) |>
  arrange(desc(n),uniprot_acc) |> slice_head(n=6)
accessions<-unique(c(lookup_je2(c("Stk1","Stp1"))$uniprot_acc,leading$uniprot_acc))
accessions<-accessions[!is.na(accessions)]
entries<-setNames(lapply(accessions,function(acc) {
  message("Public UniProt lookup: ",acc)
  tryCatch(uniprot_entry(acc),error=function(e) list(error=conditionMessage(e)))
}),accessions)
summaries<-bind_rows(lapply(names(entries),function(acc) {
  entry<-entries[[acc]]
  if(!is.null(entry$error)) return(tibble(requested_accession=acc,lookup_error=entry$error))
  crossrefs<-entry$entry$uniProtKBCrossReferences
  string_ids<-vapply(Filter(function(x) identical(x$database,"STRING"),crossrefs),function(x) x$id,character(1))
  entry$summary |> mutate(lookup_error=NA_character_,uniprot_string_crossrefs=paste(string_ids,collapse=";"),
    retrieved_utc=entry$provenance$retrieved_utc,request_url=entry$provenance$request_url)
})) |> left_join(map |> filter(entry_class=="cds",eligible_uniprot_annotation,in_total_DEA,n_JE2_good_total_DEA==1) |>
  select(requested_accession=uniprot_acc,protein_key,locus_tag,display_name,ortholog_quality),
  by="requested_accession",relationship="one-to-one")
write_tsv(summaries,file.path(out,"live_uniprot_checks.tsv"))
# Check reference-organism STRING matches against UniProt's own STRING cross-reference.
resolved<-tryCatch(string_resolve(accessions,species=93061),error=function(e) {
  message("STRING resolution unavailable: ",conditionMessage(e)); NULL
})
if(!is.null(resolved)) {
  resolved<-resolved |> left_join(select(summaries,identifier=requested_accession,uniprot_string_crossrefs),by="identifier")
  resolved$crossref_verified<-vapply(seq_len(nrow(resolved)),function(i) {
    !is.na(resolved$stringId[i]) && !is.na(resolved$uniprot_string_crossrefs[i]) &&
      resolved$stringId[i] %in% strsplit(resolved$uniprot_string_crossrefs[i],";",fixed=TRUE)[[1]]
  },logical(1))
  resolved$accepted_for_context_network<-resolved$crossref_verified & resolved$ncbiTaxonId=="93061" & resolved$n_candidates_returned==1
  write_tsv(resolved,file.path(out,"live_string_mapping_review.tsv"))
  ids<-unique(resolved$stringId[resolved$accepted_for_context_network %in% TRUE])
  if(length(ids)) {
    network<-tryCatch(string_network(ids,species=93061,network_type="functional",required_score=700),error=function(e) NULL)
    if(!is.null(network)) write_tsv(network,file.path(out,"live_string_context_edges.tsv"))
  }
}
pathway_checks<-bind_rows(lapply(unique(selected$term_id),function(id) {
  # Course terms can be mapXXXXX or organism-prefixed; lookup the reference map.
  query<-sub("^[A-Za-z]+([0-9]{5})$","map\\1",id)
  entry<-tryCatch(kegg_entry(query),error=function(e) list(status="lookup_failed",text=conditionMessage(e)))
  name_line<-grep("^NAME",strsplit(entry$text,"\n",fixed=TRUE)[[1]],value=TRUE)
  tibble(source_term_id=id,query_id=query,status=entry$status,
    returned_name=if(length(name_line)) trimws(sub("^NAME[[:space:]]+","",name_line[1])) else NA_character_,
    retrieved_utc=if(!is.null(entry$provenance)) entry$provenance$retrieved_utc else NA_character_,
    request_url=if(!is.null(entry$provenance)) entry$provenance$request_url else NA_character_)
}))
write_tsv(pathway_checks,file.path(out,"live_kegg_pathway_checks.tsv"))
writeLines(c("# Dated public-resource checks", "",
  "These calls use explicit public protein/pathway IDs only. No sample measurements are submitted.",
  "The frozen course STRING-v12 memberships remain the enrichment input; exact snapshot minor release and GO/KEGG dates are unknown.",
  "Live STRING calls are pinned to version 12.5 and reference organism NCTC8325 (taxon 93061).",
  "Only same-taxon STRING matches also verified by UniProt STRING cross-references enter the small context network.",
  "The live functional network is not a binding assay or a direct kinase-target network. Unmeasured/untested functions are not inferred from absent edges.",
  "The custom JE2 organism IDs in the slides/local snapshot are not sent as reference-organism IDs.",
  "Selection: Stk1/Stp1 plus up to six proteins most frequent in the leading edges of the top two study-significant KEGG terms in each focus contrast; ties by accession.",
  "Inspect TSV lookup status/error fields; public service failures do not change enrichment results.","",
  "Rerun with run_project_step(\"enrichment_resources\"). Responses and HTTP provenance are cached under ignored results/api_cache/."),
  file.path(out,"resource_checks.md"))
capture.output(sessionInfo(),file=file.path(out,"resource_sessionInfo.txt"))
message("Public-ID resource checks complete. Sample plots are made by enrichment_figures.")
