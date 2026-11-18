# One round of reference-organism STRING expansion; no experimental data uploaded.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"reports/step06/string_network")
local <- file.path(REPO_ROOT,"results/step06/string_network")
p <- readRDS(file.path(local,"seed_import.rds"))
species <- "93061"; threshold <- 400L; partner_limit <- 5L; max_added <- 12L
provenance <- list()
partners <- string_partners(p$seeds$string_id,species,network_type="functional",
  required_score=threshold,limit=partner_limit)
provenance[["partner_discovery"]] <- api_provenance(partners)
saveRDS(partners,file.path(local,"partner_query.rds"))
stopifnot(all(partners$ncbiTaxonId==species),all(partners$stringId_A %in% p$seeds$string_id),
  all(partners$score>=threshold/1000))
discovered <- partners |> filter(!stringId_B %in% p$seeds$string_id) |>
  distinct(string_id=stringId_B,string_name=preferredName_B)
stopifnot(!anyDuplicated(discovered$string_id))
all_ids <- unique(c(p$seeds$string_id,discovered$string_id))

# Fetch the induced network across the discovered set so partner ranking does not
# mistake the capped top-five discovery list for each protein's complete degree.
response <- .api_request(paste0(STRING_API_BASE,"/tsv/network"),
  params=list(identifiers=paste(all_ids,collapse="\r"),species=species,
    network_type="functional",required_score=threshold,add_nodes=0,
    caller_identity="JE2-Project-I-course-analysis"),service="STRING")
provenance[["discovered_set_network"]] <- response$provenance
raw_network <- .api_tsv(response)
stopifnot(all(raw_network$stringId_A %in% all_ids),all(raw_network$stringId_B %in% all_ids))
edges <- raw_network |> transmute(node_a=pmin(stringId_A,stringId_B),node_b=pmax(stringId_A,stringId_B),
  combined_score=as.numeric(score),neighborhood=as.numeric(nscore),fusion=as.numeric(fscore),
  cooccurrence=as.numeric(pscore),coexpression=as.numeric(ascore),experimental=as.numeric(escore),
  database=as.numeric(dscore),textmining=as.numeric(tscore),edge_origin="STRING12.5 one-round API") |> distinct()
stopifnot(!anyDuplicated(edges[c("node_a","node_b")]),all(edges$combined_score>=.4))

# Resolve added reference IDs using exact UniProt STRING cross-references. Querying
# names would confuse paralogs; the returned STRING tokens are explicitly matched.
crossrefs <- list()
chunks <- split(discovered$string_id,ceiling(seq_len(nrow(discovered))/20))
for(i in seq_along(chunks)) {
  message("Reference mapping batch ",i,"/",length(chunks))
  r <- .api_request("https://rest.uniprot.org/uniprotkb/search",params=list(
    query=paste(paste0("xref:string-",chunks[[i]]),collapse=" OR "),
    format="tsv",fields="accession,xref_string",size=500),service="UniProt")
  provenance[[paste0("uniprot_xrefs_",i)]] <- r$provenance
  z <- .api_tsv(r)
  stopifnot(nrow(z)<500L)
  if(nrow(z)) crossrefs[[i]] <- z |> transmute(uniprot_acc=Entry,string_id=strsplit(STRING,";",fixed=TRUE)) |>
    unnest(string_id) |> mutate(string_id=trimws(string_id)) |> filter(string_id %in% chunks[[i]]) |> distinct()
}
xref <- bind_rows(crossrefs)
map <- read_identifier_map()
reference <- xref |> left_join(map |> select(uniprot_acc,protein_key,locus_tag,display_name,product,
  eligible_uniprot_annotation,n_JE2_good_total_DEA,in_total_DEA),by="uniprot_acc",relationship="many-to-many")
reference_counts <- reference |> group_by(string_id) |> summarise(n_uniprot=n_distinct(uniprot_acc),
  n_je2=n_distinct(protein_key,na.rm=TRUE),.groups="drop")
accepted_reference <- reference |> inner_join(filter(reference_counts,n_uniprot==1L,n_je2<=1L),by="string_id",relationship="many-to-one") |>
  distinct(string_id,.keep_all=TRUE) |> select(-n_uniprot,-n_je2)
added <- discovered |> left_join(accepted_reference,by="string_id",relationship="one-to-one") |>
  left_join(select(reference_counts,string_id,n_uniprot,n_je2),by="string_id",relationship="one-to-one") |>
  mutate(display_name=coalesce(display_name,string_name),
    mapping_status=case_when(is.na(n_uniprot)~"No exact UniProt cross-reference returned",
      n_uniprot>1 | n_je2>1~"Ambiguous reference/JE2 mapping; no experimental values assigned",
      is.na(protein_key)~"Exact reference cross-reference; no JE2 map",
      TRUE~"Exact UniProt STRING cross-reference and unique JE2 map"),
    node_role="Added partner",measured_in_total=if_else(!is.na(protein_key),coalesce(in_total_DEA,FALSE),NA))
write_csv(reference,file.path(out,"added_reference_mapping_all_candidates.csv"))

# Annotate all discovered candidates from the full course mapping, not just the
# selected proteins. Lack of membership is not proof of absent biological function.
inventory <- read_csv(file.path(REPO_ROOT,"reports/step06/01_module_inventory.csv"),show_col_types=FALSE)
module_map <- inventory |> group_by(protein_key) |> summarise(
  module_group=paste(unique(module),collapse="; "),annotation_review=any(annotation_review),.groups="drop")
added <- added |> left_join(module_map,by="protein_key",relationship="many-to-one") |>
  mutate(module_group=coalesce(module_group,"Other or unassigned"),annotation_review=coalesce(annotation_review,FALSE))
source <- read_tsv(file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),show_col_types=FALSE)
codes <- tibble(contrast=c("interaction_stk_vs_anc_at_TSB_vs_SN","interaction_stp_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_stk_at_TSB_vs_SN"),contrast_code=c("stk_anc","stp_anc","stp_stk"))
effects <- source |> inner_join(codes,by="contrast",relationship="many-to-one") |>
  mutate(halo_effect=if_else(estimate_type=="observed" & coalesce(all_relevant_groups_2of3,FALSE) &
    coalesce(eligible_uniprot_annotation,FALSE) & coalesce(n_JE2_good_total_DEA==1,FALSE) & is.finite(diff),diff,NA_real_),
    halo_FDR=if_else(is.finite(halo_effect),FDR,NA_real_)) |>
  select(protein_key,contrast_code,diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR) |>
  pivot_wider(names_from=contrast_code,values_from=c(diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR))
added <- added |> left_join(effects,by="protein_key",relationship="many-to-one") |>
  mutate(has_quantitative_halo=is.finite(halo_effect_stk_anc)|is.finite(halo_effect_stp_anc)|is.finite(halo_effect_stp_stk))

# Rank ONLY by network context and mapping/measurement availability, never by NES
# or effect magnitude/significance. Degree here is within this discovered network.
seed_ids <- p$seeds$string_id
anchor_ids <- filter(p$seeds,node_role=="Biological anchor")$string_id
candidate_seed_edges <- bind_rows(
  edges |> filter(node_a %in% discovered$string_id,node_b %in% seed_ids) |>
    transmute(string_id=node_a,seed_id=node_b,combined_score,experimental,database,textmining),
  edges |> filter(node_b %in% discovered$string_id,node_a %in% seed_ids) |>
    transmute(string_id=node_b,seed_id=node_a,combined_score,experimental,database,textmining)) |>
  left_join(p$seeds |> select(seed_id=string_id,seed_module=module_group),by="seed_id",relationship="many-to-one")
context <- candidate_seed_edges |> group_by(string_id) |> summarise(n_seed_links=n_distinct(seed_id),
  n_anchor_links=n_distinct(seed_id[seed_id %in% anchor_ids]),
  n_seed_modules=n_distinct(seed_module[!seed_id %in% anchor_ids]),
  best_score=max(combined_score),n_experimental_links=sum(experimental>0),
  seed_names=paste(p$seeds$display_name[match(seed_id,p$seeds$string_id)],collapse="; "),.groups="drop")
added <- added |> left_join(context,by="string_id",relationship="one-to-one") |>
  mutate(priority_class=case_when(n_anchor_links>0~"1: anchor neighbour",
    n_seed_modules>=2~"2: connects selected modules",n_seed_links>=2~"3: connects multiple seeds",
    TRUE~"4: single-seed context")) |>
  arrange(priority_class,desc(has_quantitative_halo),desc(n_seed_modules),desc(n_seed_links),desc(best_score),string_id)

# Compact display budget: first preserve one best contextual partner per module,
# then fill remaining slots by the declared network-context ranking.
coverage <- character()
for(mod in unique(filter(p$seeds,node_role=="Selected seed")$module_group)) {
  linked <- candidate_seed_edges |> filter(seed_module==mod) |> pull(string_id)
  # Prefer a supplied member of the same pathway, when one was discovered;
  # a broadly connected anchor neighbour must not absorb every module slot.
  same_module <- added |> filter(string_id %in% linked,grepl(mod,module_group,fixed=TRUE))
  eligible <- if(nrow(same_module)) same_module$string_id else added$string_id[added$string_id %in% linked]
  if(length(eligible)) coverage <- c(coverage,eligible[1])
}
coverage <- unique(coverage)
chosen <- head(unique(c(coverage,added$string_id)),max_added)
added <- added |> mutate(selected_for_compact=string_id %in% chosen,
  selection_reason=case_when(!selected_for_compact~"Full expansion retained; outside compact display budget",
    string_id %in% coverage~"Best-ranked partner covering a selected module",
    n_anchor_links>0~"Direct anchor neighbour in one-round functional network",
    TRUE~"Network-context rank within compact display budget"))
write_csv(added,file.path(out,"expansion_candidate_decisions.csv"))
write_csv(candidate_seed_edges,file.path(out,"candidate_to_seed_links.csv"))
seed_edges_live <- edges |> filter(node_a %in% seed_ids,node_b %in% seed_ids)
# Compare same-version API seed edges with the student's export; retain export
# baseline regardless of a mismatch and expose differences instead of overwriting.
edge_key <- function(x) paste(x$node_a,x$node_b,sep="|")
differences <- full_join(p$edges |> select(node_a,node_b,export_score=combined_score),
  seed_edges_live |> select(node_a,node_b,api_score=combined_score),by=c("node_a","node_b")) |>
  filter(is.na(export_score)|is.na(api_score)|abs(export_score-api_score)>1e-6)
write_csv(differences,file.path(out,"export_vs_live_seed_edge_differences.csv"))
all_nodes <- bind_rows(p$seeds,added)
stopifnot(!anyDuplicated(all_nodes$string_id))
nodes <- filter(all_nodes,string_id %in% c(seed_ids,chosen))
expanded <- bind_rows(p$edges,edges |> filter(!node_a %in% seed_ids | !node_b %in% seed_ids)) |>
  filter(node_a %in% nodes$string_id,node_b %in% nodes$string_id)
stopifnot(!anyDuplicated(expanded[c("node_a","node_b")]),all(seed_ids %in% nodes$string_id),nrow(nodes)<=61L)
degree <- function(e,ids) as.integer(table(factor(c(e$node_a,e$node_b),levels=ids)))
nodes$degree_expanded <- degree(expanded,nodes$string_id)
write_csv(nodes,file.path(out,"nodes.csv"))
write_csv(expanded,file.path(out,"edges_expanded.csv"))
write_csv(all_nodes,file.path(out,"nodes_full_one_round.csv"))
write_csv(edges,file.path(out,"edges_full_one_round.csv"))
saveRDS(list(nodes=all_nodes,edges=edges,partners=partners,provenance=provenance),file.path(local,"expansion_full.rds"))
prov <- bind_rows(lapply(names(provenance),function(key) {
  r <- provenance[[key]]
  tibble(query=key,url=r$request_url,retrieved_utc=r$retrieved_utc,http_status=r$http_status,cache_hit=r$cache_hit)
}))
write_csv(prov,file.path(out,"API_provenance.csv"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
message(nrow(discovered)," novel partners discovered; ",length(chosen)," retained in compact overview; ",
  nrow(nodes)," total compact nodes; ",nrow(expanded)," edges.")
print(filter(added,selected_for_compact) |> select(display_name,uniprot_acc,module_group,priority_class,
  n_seed_links,n_anchor_links,has_quantitative_halo,selection_reason),n=20,width=170)
