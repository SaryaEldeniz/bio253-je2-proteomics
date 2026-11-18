# The student's figure-based selection replaces automatic shortlist eligibility.
# Preserve every requested ID, including flagged module memberships.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
input <- file.path(REPO_ROOT,"scripts/06_network/user_selected_uniprot_ids.txt")
ids <- trimws(readLines(input,warn=FALSE))
stopifnot(length(ids)>0L,!anyNA(ids),all(nzchar(ids)),!anyDuplicated(ids))
out <- file.path(REPO_ROOT,"reports/step06/string_uploads/user_selected")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inventory <- read_csv(file.path(REPO_ROOT,"reports/step06/01_module_inventory.csv"),show_col_types=FALSE)
members <- tibble(uniprot_acc=ids,selection_order=seq_along(ids)) |>
  left_join(inventory,by="uniprot_acc",relationship="one-to-many")
stopifnot(!anyNA(members$protein_key))
id_check <- members |> distinct(uniprot_acc,protein_key) |> count(uniprot_acc)
stopifnot(all(id_check$n==1L))
members <- members |> arrange(selection_order)
write_csv(members,file.path(out,"selected_module_membership.csv"))
anchors <- lookup_je2(c("Stk1","Stp1")) |>
  select(protein_key,locus_tag,display_name,uniprot_acc,product)
stopifnot(nrow(anchors)==2L,!anyNA(anchors$uniprot_acc))
manifest <- members |> group_by(selection_order,uniprot_acc,protein_key,locus_tag,display_name,product) |>
  summarise(modules=paste(unique(module),collapse="; "),
    annotation_review=any(annotation_review),
    annotation_notes=paste(unique(na.omit(annotation_note)),collapse="; "),.groups="drop") |>
  arrange(selection_order) |> mutate(upload_role="User selected from sample plots")
write_csv(manifest,file.path(out,"selected_manifest.csv"))
write_csv(filter(manifest,annotation_review),file.path(out,"annotation_review_retained.csv"))
writeLines(ids,file.path(out,"combined_selected.txt"))
writeLines(unique(c(ids,anchors$uniprot_acc)),file.path(out,"combined_selected_plus_anchors.txt"))
stems <- c("purine","oxidative_phosphorylation","nitrogen","aminoacyl_tRNA")
summary <- list()
for (stem in stems) {
  x <- members |> filter(file_stem==stem) |> distinct(uniprot_acc,.keep_all=TRUE)
  writeLines(x$uniprot_acc,file.path(out,paste0(stem,"_selected.txt")))
  writeLines(unique(c(x$uniprot_acc,anchors$uniprot_acc)),file.path(out,paste0(stem,"_selected_plus_anchors.txt")))
  summary[[stem]] <- tibble(module=x$module[1],file_stem=stem,selected=nrow(x),
    with_anchors=length(unique(c(x$uniprot_acc,anchors$uniprot_acc))),
    retained_annotation_flags=sum(x$annotation_review))
}
summary <- bind_rows(summary)
write_csv(summary,file.path(out,"list_counts.csv"))

# Make a corresponding halo table for this exact user selection, including all
# three original estimates for the previously held annotation-review proteins.
nodes <- manifest |> select(-selection_order) |>
  bind_rows(anchors |> mutate(modules=NA_character_,annotation_review=FALSE,annotation_notes="",
    upload_role="Biological anchor"))
stopifnot(nrow(nodes)==length(unique(c(ids,anchors$uniprot_acc))),!anyDuplicated(nodes$uniprot_acc))
codes <- tibble(contrast=c("interaction_stk_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_anc_at_TSB_vs_SN","interaction_stp_vs_stk_at_TSB_vs_SN"),
  contrast_code=c("stk_anc","stp_anc","stp_stk"))
effects <- read_tsv(file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),show_col_types=FALSE) |>
  inner_join(codes,by="contrast",relationship="many-to-one") |>
  semi_join(nodes,by="protein_key") |>
  mutate(halo_effect=if_else(estimate_type=="observed" & coalesce(all_relevant_groups_2of3,FALSE) & is.finite(diff),diff,NA_real_),
    halo_FDR=if_else(is.finite(halo_effect),FDR,NA_real_)) |>
  select(protein_key,contrast_code,diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR)
stopifnot(nrow(effects)==nrow(nodes)*nrow(codes),!anyDuplicated(effects[c("protein_key","contrast_code")]))
nodes <- nodes |> left_join(effects |> pivot_wider(names_from=contrast_code,
  values_from=c(diff,FDR,estimate_type,n_observed_min_relevant,halo_effect,halo_FDR)),
  by="protein_key",relationship="one-to-one")
write_csv(nodes,file.path(out,"network_node_attributes.csv"))
write_csv(nodes |> select(protein_key,locus_tag,display_name,uniprot_acc,upload_role) |>
  mutate(requested_STRING_version="12.5",intended_reference_taxon=93061L,
    returned_STRING_ID=NA_character_,returned_taxon=NA_character_,returned_name=NA_character_,
    mapping_decision="Pending manual STRING mapping review"),
  file.path(out,"STRING_mapping_review_template.csv"))

# Verify exported plain text exactly, including requested ordering and no omissions.
stopifnot(identical(readLines(file.path(out,"combined_selected.txt")),ids),
  setequal(readLines(file.path(out,"combined_selected_plus_anchors.txt")),nodes$uniprot_acc),
  setequal(unique(members$uniprot_acc),ids))
md <- c("# User-selected STRING upload lists","",
  paste0("These lists preserve all ",length(ids)," UniProt IDs selected by the student from the sample plots. No additional NES, FDR, effect-size or automatic-shortlist filter was applied. The previous automatic lists remain unchanged."),"",
  "Files ending _selected.txt contain only the chosen proteins. Files ending _selected_plus_anchors.txt additionally contain Stk1 (Q2FZ64) and Stp1 (Q2FZ65). Plain-text uploads have one accession per line and no header.","",
  "| Module | Selected | With Stk1/Stp1 | Membership flags retained |","|---|---:|---:|---:|")
for(i in seq_len(nrow(summary))) md <- c(md,paste0("| ",summary$module[i]," | ",summary$selected[i],
  " | ",summary$with_anchors[i]," | ",summary$retained_annotation_flags[i]," |"))
md <- c(md,"",paste0("Combined: ",length(ids)," chosen proteins, or ",nrow(nodes)," with both anchors. All IDs were matched uniquely to JE2 entries in the existing inventory."),"",
  "Seven selected proteins retain the previous purine-membership review flag: Q2G010, Q2FYK7, Q2FYG2, Q2FUT0, Q2FV52, Q2FXM9 and Q2FWF8. They are included in the purine and combined lists exactly as requested. This flag concerns pathway assignment, not whether the protein is biologically interesting. See annotation_review_retained.csv for product descriptions and reasons.","",
  paste0("selected_manifest.csv records names, JE2 loci, modules and selection origin. network_node_attributes.csv contains this exact ",nrow(nodes),"-node set with the three original effects/FDRs and halo-ready columns. Source-LOD or <2/3-supported estimates have NA halo values; original estimates remain available with flags."),"",
  "Review STRING's returned organism and mapping. Keep the reference organism namespace consistent; do not mix reference-ortholog STRING IDs with custom JE2 STRG IDs. No STRING requests or network expansion were performed by this script.","",
  "Reproduce from Sarya_George.Rproj:","","```r",
  'source("scripts/06_network/06_export_user_selection.R")',"```", "",
  "Input: scripts/06_network/user_selected_uniprot_ids.txt. Supporting inputs are the step-6 membership inventory and original October annotated total-proteome contrasts.")
writeLines(md,file.path(out,"README.md"))
print(summary,width=130)
message("Exported exact user selection: ",length(ids)," IDs; combined with anchors: ",nrow(nodes)," IDs.")
