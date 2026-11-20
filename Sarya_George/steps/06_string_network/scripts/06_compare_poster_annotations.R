# Compare source course scope with the fresh STRING reference KEGG memberships.
# This is an annotation comparison, not a new enrichment test or NES correction.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/06_string_network/reports/poster_networks")
raw <- file.path(REPO_ROOT,"data/raw/step06_string/poster_api")
course <- read_csv(file.path(out,"pathway_membership.csv"),show_col_types=FALSE) |>
  mutate(set_id=as.character(set_id))
dictionary <- tibble(term_id=c("map00230","map00190","map00910","map00970"),
  module=c("Purine metabolism","Oxidative phosphorylation","Nitrogen metabolism","Aminoacyl-tRNA biosynthesis"),
  code=c("P","R","N","T"))
live_all <- comparisons <- summaries <- list()
sets <- c("49","67","75")
sets <- sets[file.exists(file.path(out,paste0("nodes_",sets,".csv")))]
for(set in sets) {
  nodes <- read_csv(file.path(out,paste0("nodes_",set,".csv")),show_col_types=FALSE)
  a <- read_tsv(file.path(raw,paste0("functional_annotations_",set,".tsv")),show_col_types=FALSE)
  live <- a |> filter(category=="KEGG",term %in% dictionary$term_id) |>
    transmute(term_id=term,string_id=strsplit(inputGenes,",",fixed=TRUE)) |> unnest(string_id) |>
    mutate(string_id=trimws(string_id)) |>
    left_join(dictionary,by="term_id",relationship="many-to-one") |>
    left_join(nodes |> select(string_id,protein_key,uniprot_acc,locus_tag,display_name),
      by="string_id",relationship="many-to-one") |>
    mutate(set_id=set,annotation_review=FALSE,annotation_source="Fresh STRING12.5 reference KEGG")
  stopifnot(!anyNA(live$protein_key),!anyDuplicated(live[c("term_id","string_id")]))
  live_all[[set]] <- live
  old <- course |> filter(set_id==set) |>
    select(term_id,string_id,uniprot_acc,display_name,course_annotation_review=annotation_review) |>
    mutate(in_course=TRUE)
  cmp <- full_join(old,live |> select(term_id,string_id) |> mutate(in_live=TRUE),by=c("term_id","string_id")) |>
    mutate(set_id=set,in_course=coalesce(in_course,FALSE),in_live=coalesce(in_live,FALSE)) |>
    left_join(dictionary,by="term_id",relationship="many-to-one")
  comparisons[[set]] <- cmp
  summaries[[set]] <- cmp |> group_by(set_id,term_id,module) |>
    summarise(course_members=sum(in_course),live_members=sum(in_live),shared=sum(in_course & in_live),
      course_only=sum(in_course & !in_live),live_only=sum(in_live & !in_course),.groups="drop")
  tags <- live |> group_by(protein_key) |> summarise(
    live_module_group=paste(unique(module),collapse="; "),live_pathway_codes=paste(unique(code),collapse="/"),.groups="drop")
  live_nodes <- nodes |> rename(course_module_group=module_group,course_annotation_review=annotation_review,
    course_pathway_codes=pathway_codes) |>
    left_join(tags,by="protein_key",relationship="one-to-one") |>
    mutate(module_group=if_else(node_role=="Biological anchor","Regulatory anchors",
        coalesce(live_module_group,"Other or unassigned")),
      pathway_codes=coalesce(live_pathway_codes,if_else(node_role=="Biological anchor","","C")),
      annotation_review=FALSE,annotation_group_source="Fresh STRING12.5 reference KEGG")
  write_csv(live_nodes,file.path(out,paste0("nodes_",set,"_live.csv")))
}
write_csv(bind_rows(live_all),file.path(out,"pathway_membership_live.csv"))
write_csv(bind_rows(comparisons),file.path(out,"course_vs_live_KEGG_membership.csv"))
write_csv(bind_rows(summaries),file.path(out,"course_vs_live_KEGG_counts.csv"))
print(bind_rows(summaries),width=160)
message("Live and course annotations preserved separately. No historical pathway results changed.")
