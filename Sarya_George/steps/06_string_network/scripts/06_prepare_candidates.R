# Step 6, parts 1-4: contrast-first candidates and manual STRING upload handoff.
# No NES, GSEA leading edges or pathway significance enter candidate selection.
# No DEA models/FDRs are recalculated and no requests are sent to STRING.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT, "steps/06_string_network/reports")
local_out <- file.path(REPO_ROOT, "results", "step06")
uploads <- file.path(out, "string_uploads")
for (p in c(out, local_out, uploads)) dir.create(p, recursive=TRUE, showWarnings=FALSE)

plan <- step04_contrast_plan() |> filter(contrast_type == "interaction") |>
  mutate(contrast_code = c("stk_anc", "stp_anc", "stp_stk")[match(contrast, c(
    "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
    "interaction_stp_vs_stk_at_TSB_vs_SN"))])
stopifnot(nrow(plan)==3L, !anyNA(plan$contrast_code))
modules <- tibble(term_id=c("map00230","map00190","map00910","map00970"),
  module=c("Purine metabolism","Oxidative phosphorylation","Nitrogen metabolism","Aminoacyl-tRNA biosynthesis"),
  file_stem=c("purine","oxidative_phosphorylation","nitrogen","aminoacyl_tRNA"))
input_files <- c(source=file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),
  terms=file.path(REPO_ROOT,"results/step03/terms_cds.tsv.gz"),
  mapping=file.path(REPO_ROOT,"results/step03/id_map.tsv"),
  abundance=check_derived("abundances_total.tsv.gz"))
source_data <- read_tsv(input_files[["source"]], show_col_types=FALSE)
terms <- read_tsv(input_files[["terms"]], show_col_types=FALSE)
map <- read_identifier_map()
membership <- terms |> filter(category=="KEGG (Kyoto Encyclopedia of Genes and Genomes)") |>
  inner_join(modules, by=c("term"="term_id"), relationship="many-to-one") |>
  transmute(protein_key, term_id=term, module, file_stem) |> distinct()
stopifnot(!anyDuplicated(membership[c("protein_key","term_id")]))

# Module-specific, explicit review flags; an empty flag does not validate membership.
review_file <- file.path(REPO_ROOT,"steps/06_string_network/scripts/module_annotation_review.csv")
reviews <- read_csv(review_file, show_col_types=FALSE)
input_files <- c(input_files, annotation_review=review_file)
stopifnot(!anyDuplicated(reviews[c("term_id","locus_tag")]))
base <- source_data |> filter(contrast %in% plan$contrast) |>
  left_join(select(plan, contrast, contrast_code), by="contrast", relationship="many-to-one") |>
  mutate(mapping_ok=coalesce(eligible_uniprot_annotation,FALSE) &
    coalesce(n_JE2_good_total_DEA==1,FALSE) & !is.na(uniprot_acc),
    observed_model=coalesce(estimate_type=="observed",FALSE),
    fdr_supported=is.finite(FDR) & FDR<.1,
    large_effect_supported=fdr_supported & is.finite(diff) & abs(diff)>1,
    quantitative_eligible=eligible_cds & !is_contaminant & observed_model &
      coalesce(all_relevant_groups_2of3,FALSE) & mapping_ok & is.finite(diff) & is.finite(FDR),
    supported_for_selection=quantitative_eligible & fdr_supported,
    large_for_selection=quantitative_eligible & large_effect_supported)
stopifnot(!anyDuplicated(base[c("protein_key","contrast")]))
write_tsv(plan, file.path(out,"contrast_direction_key.tsv"))

# Retain all annotated proteins, including members absent from the DEA.
inventory <- membership |> left_join(select(map,protein_key,locus_tag,display_name,
  uniprot_acc,product,ortholog_quality,eligible_uniprot_annotation,n_JE2_good_total_DEA,
  in_total_DEA), by="protein_key", relationship="many-to-one") |>
  left_join(reviews, by=c("term_id","locus_tag"), relationship="many-to-one") |>
  mutate(annotation_review=!is.na(annotation_note),
    in_interaction_data=protein_key %in% base$protein_key)
write_csv(inventory, file.path(out,"01_module_inventory.csv"))

all_interactions <- inventory |> select(protein_key,term_id,module,file_stem,annotation_review,annotation_note) |>
  inner_join(base, by="protein_key", relationship="many-to-many") |>
  mutate(selection_status=case_when(
    !eligible_cds | is_contaminant ~ "Not an eligible CDS",
    !observed_model ~ "LOD/other source model: not quantitative selection",
    !mapping_ok ~ "UniProt mapping review",
    !coalesce(all_relevant_groups_2of3,FALSE) ~ "Fewer than 2/3 measured in a relevant group",
    !is.finite(diff) | !is.finite(FDR) ~ "No finite estimate/FDR",
    !fdr_supported ~ "Does not pass protein FDR<0.1",
    annotation_review ~ "Protein contrast supported; module membership review",
    large_effect_supported ~ "Large-effect candidate",
    TRUE ~ "Supported smaller-effect candidate"),
    candidate_in_this_contrast=supported_for_selection & !annotation_review)
write_csv(all_interactions, file.path(out,"02_all_module_interactions.csv"))

# Within-genotype contrasts are TSB-SN in the source. Reverse effect AND CI endpoints
# to show SN-TSB, but keep the supplied FDR unchanged.
environment <- source_data |> filter(contrast %in% paste0("TSB_vs_SN_at_",c("anc","stk","stp"))) |>
  transmute(protein_key,genotype=sub("TSB_vs_SN_at_","",contrast),
    SN_minus_TSB=-diff,response_CI_low=-conf.high,response_CI_high=-conf.low,
    response_FDR=FDR,response_model=estimate_type,
    response_n_min=n_observed_min_relevant,response_detection=group_detection_pattern)
responses_wide <- environment |> select(protein_key,genotype,SN_minus_TSB) |>
  pivot_wider(names_from=genotype,values_from=SN_minus_TSB,names_prefix="R_")
all_interactions <- all_interactions |> left_join(responses_wide, by="protein_key",relationship="many-to-one") |>
  mutate(reconstructed_effect=case_when(contrast_code=="stk_anc"~R_anc-R_stk,
    contrast_code=="stp_anc"~R_anc-R_stp,contrast_code=="stp_stk"~R_stk-R_stp))
stopifnot(all(abs(all_interactions$diff-all_interactions$reconstructed_effect)<1e-8,na.rm=TRUE))
write_csv(all_interactions, file.path(out,"02_all_module_interactions.csv"))
write_csv(environment |> semi_join(membership,by="protein_key") |>
  left_join(select(map,protein_key,locus_tag,display_name,uniprot_acc),by="protein_key",relationship="many-to-one"),
  file.path(out,"03_within_genotype_SN_minus_TSB.csv"))

# Eligibility is 'at least one supported interaction', never best NES or top N.
candidate_summary <- all_interactions |> group_by(term_id,module,file_stem,protein_key,locus_tag,
    display_name,uniprot_acc,product,annotation_review,annotation_note) |>
  summarise(n_supported_interactions=sum(supported_for_selection),
    n_large_effect_interactions=sum(large_for_selection),
    supported_contrasts=paste(contrast_code[supported_for_selection],collapse=";"),
    large_effect_contrasts=paste(contrast_code[large_for_selection],collapse=";"),
    all_six_groups_3of3=all(n_observed_min_relevant==3),
    .groups="drop") |>
  mutate(selected=n_supported_interactions>0 & !annotation_review,
    tier=case_when(annotation_review & n_supported_interactions>0 ~ "Hold: module annotation review",
      n_large_effect_interactions>0 ~ "A: large-effect candidate",
      n_supported_interactions>0 ~ "B: supported smaller-effect candidate",
      TRUE ~ "Not shortlisted"),
    selection_reason=case_when(selected & n_large_effect_interactions>0 ~
      "At least one observed interaction has FDR<0.1, |effect|>1 and >=2/3 measured per relevant group; singleton UniProt",
      selected ~ "At least one observed interaction has FDR<0.1 and >=2/3 measured per relevant group; singleton UniProt; no minimum effect cutoff",
      annotation_review & n_supported_interactions>0 ~ "Contrast supported, but supplied module membership held for review",
      TRUE ~ "No interaction meeting all declared step-6 contrast/support/ID criteria"))
wide <- all_interactions |> select(protein_key,term_id,contrast_code,diff,FDR,conf.low,conf.high,
    estimate_type,n_observed_min_relevant,group_detection_pattern,direction,supported_for_selection) |>
  pivot_wider(names_from=contrast_code,values_from=c(diff,FDR,conf.low,conf.high,estimate_type,
    n_observed_min_relevant,group_detection_pattern,direction,supported_for_selection))
candidates <- candidate_summary |> left_join(wide,by=c("protein_key","term_id"),relationship="one-to-one") |>
  left_join(responses_wide,by="protein_key",relationship="many-to-one") |>
  arrange(match(term_id,modules$term_id),tier,display_name,locus_tag)
write_csv(candidates, file.path(out,"04_candidate_decision_table.csv"))
write_tsv(candidates, file.path(out,"04_candidate_decision_table.tsv"))
write_csv(filter(candidates,selected), file.path(out,"05_selected_candidates.csv"))
write_csv(filter(candidates,annotation_review), file.path(out,"06_annotation_review.csv"))

# Off-module proteins are not lost: provide the full contrast-selected context list.
global <- base |> filter(supported_for_selection) |> mutate(
  in_four_modules=protein_key %in% membership$protein_key)
write_csv(global, file.path(out,"07_global_contrast_supported_context.csv"))

# Mandatory biological anchors are not presented as statistical discoveries.
anchors <- lookup_je2(c("Stk1","Stp1")) |>
  transmute(protein_key,locus_tag,display_name,uniprot_acc,product,
    upload_role="Biological anchor: included regardless of interaction significance")
stopifnot(nrow(anchors)==2L,all(anchors$uniprot_acc %in% c("Q2FZ64","Q2FZ65")))
seeds <- filter(candidates,selected) |> distinct(protein_key,locus_tag,display_name,uniprot_acc,product) |>
  mutate(upload_role="Contrast-selected module protein") |> bind_rows(anchors) |>
  distinct(uniprot_acc,.keep_all=TRUE)
stopifnot(!anyNA(seeds$uniprot_acc),!anyDuplicated(seeds$uniprot_acc))
write_csv(seeds,file.path(uploads,"upload_manifest.csv"))
writeLines(seeds$uniprot_acc,file.path(uploads,"all_selected_plus_Stk1_Stp1.txt"))
writeLines(anchors$uniprot_acc,file.path(uploads,"anchors_Stk1_Stp1.txt"))
held_ids <- filter(candidates,annotation_review,n_supported_interactions>0) |> pull(uniprot_acc) |> unique()
writeLines(unique(c(held_ids,anchors$uniprot_acc)),file.path(uploads,"membership_review_plus_anchors.txt"))
for (i in seq_len(nrow(modules))) {
  ids <- filter(candidates,selected,term_id==modules$term_id[i]) |> pull(uniprot_acc) |> unique()
  writeLines(unique(c(ids,anchors$uniprot_acc)),file.path(uploads,paste0(modules$file_stem[i],"_plus_anchors.txt")))
}
large_ids <- filter(candidates,selected,n_large_effect_interactions>0) |> pull(uniprot_acc) |> unique()
writeLines(unique(c(large_ids,anchors$uniprot_acc)),file.path(uploads,"large_effect_subset_plus_anchors.txt"))

# An editable mapping form for the user's STRING 12.5 export. No guessed STRING IDs.
mapping_review <- seeds |> mutate(requested_STRING_version="12.5",
  intended_reference_taxon=93061L, returned_STRING_ID=NA_character_,
  returned_taxon=NA_character_,returned_name=NA_character_,
  mapping_decision="Pending manual STRING mapping review",mapping_notes=NA_character_)
write_csv(mapping_review,file.path(uploads,"STRING_mapping_review_template.csv"))

# All sample-level values stay in ignored results/, not versioned reports/.
values <- read_abundances("total") |> filter(medium %in% c("TSB","SN")) |>
  semi_join(bind_rows(select(membership,protein_key),select(anchors,protein_key)),by="protein_key") |>
  left_join(select(map,protein_key,display_name,uniprot_acc),by="protein_key",relationship="many-to-one")
stopifnot(!anyDuplicated(values[c("protein_key","sample")]),
  all(is.na(values$normalized_abundance) | is.finite(values$normalized_abundance)))
saveRDS(values,file.path(local_out,"sample_values_TSB_SN.rds"))
saveRDS(list(candidates=candidates,interactions=all_interactions,environment=environment,
  values=values,anchors=anchors,modules=modules,base=base),file.path(local_out,"step06_prepared.rds"))

# Independent joins/detection and contrast algebra checks on actual measurements.
counts <- values |> group_by(protein_key,group) |>
  summarise(n=sum(is.finite(normalized_abundance)),.groups="drop")
for (name in plan$contrast) {
  wanted <- names(contrast_weights(name))[contrast_weights(name)!=0]
  observed <- counts |> filter(as.character(group) %in% wanted) |> group_by(protein_key) |>
    summarise(recounted_min=min(n),.groups="drop")
  test <- base |> filter(contrast==name) |> inner_join(observed,by="protein_key",relationship="one-to-one")
  stopifnot(all(test$n_observed_min_relevant==test$recounted_min))
}
closure <- base |> select(protein_key,contrast_code,diff) |> pivot_wider(names_from=contrast_code,values_from=diff)
stopifnot(all(abs(closure$stp_stk-(closure$stp_anc-closure$stk_anc))<1e-8,na.rm=TRUE))
stopifnot(all(filter(all_interactions,candidate_in_this_contrast)$estimate_type=="observed"),
  all(filter(all_interactions,candidate_in_this_contrast)$n_observed_min_relevant>=2),
  all(filter(all_interactions,candidate_in_this_contrast)$FDR<.1))
module_counts <- inventory |> group_by(term_id,module) |>
  summarise(annotated_members=n_distinct(protein_key),in_interaction_data=sum(in_interaction_data),.groups="drop") |>
  left_join(candidates |> group_by(term_id) |> summarise(
    large_effect=n_distinct(protein_key[selected & n_large_effect_interactions>0]),
    selected=n_distinct(protein_key[selected]),
    annotation_held=n_distinct(protein_key[n_supported_interactions>0 & annotation_review]),.groups="drop"),by="term_id")
stopifnot(sum(module_counts$large_effect)==sum(candidates$selected & candidates$n_large_effect_interactions>0))
write_csv(module_counts,file.path(out,"module_selection_counts.csv"))
public_input_paths <- vapply(unname(input_files),function(path) {
  project_prefix <- paste0(REPO_ROOT,"/")
  if(startsWith(path,project_prefix)) substring(path,nchar(project_prefix)+1L) else
    paste0("course-data/derived/",basename(path))
},character(1))
write_tsv(tibble(file=names(input_files),path=public_input_paths,md5=unname(tools::md5sum(input_files))),
  file.path(out,"input_provenance.tsv"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines(c("Selection uses protein contrasts only; no NES/GSEA inputs are read.",
  "FDR<0.1; source observed; >=2/3 detections in relevant groups; eligible CDS; unique annotation-eligible UniProt.",
  "Tier A also requires |interaction effect|>1; Tier B retains smaller supported effects.",
  "FDR values are supplied per-contrast results, not a new correction across the three selected contrasts.",
  "The three interactions are linearly dependent; multi-contrast support is not independent replication.",
  "Annotation flags hold a module/protein pair for review, not a proof of an incorrect assignment.",
  "Stk1 and Stp1 are biological anchors, not automatically significant candidates.",
  "No STRING calls or expansion performed. User will export version 12.5 after reviewing mappings.",
  "Checks passed: unique joins/IDs, source-observed selections, detection recount, reconstructed effects, contrast closure."),
  file.path(out,"selection_rules_and_checks.txt"))
print(module_counts,width=150)
message(nrow(seeds)-2L," unique contrast-selected proteins plus 2 anchors; ",nrow(seeds)," upload IDs.")
