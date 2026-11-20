# All leading-edge proteins in the four KEGG modules selected for the summary slide.
# Leading edges are statistical contributors, not a causal or per-protein NES decomposition.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out<-file.path(REPO_ROOT,"steps/05_enrichment/reports/observed_only/four_module_drivers")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
contrasts<-c("interaction_stk_vs_anc_at_TSB_vs_SN","interaction_stp_vs_anc_at_TSB_vs_SN","interaction_stp_vs_stk_at_TSB_vs_SN")
modules<-tibble(term_id=c("map00230","map00190","map00910","map00970"),
  module=c("Purine metabolism","Oxidative phosphorylation","Nitrogen metabolism","Aminoacyl-tRNA biosynthesis"))
gsea<-read_enrichment("GSEA",policy="observed_only",significant_only=FALSE,collection="KEGG") |>
  filter(contrast %in% contrasts) |> inner_join(modules,by="term_id",relationship="many-to-one")
stopifnot(nrow(gsea)==12L,!anyDuplicated(gsea[c("contrast","term_id")]))
source_data<-read_tsv(file.path(REPO_ROOT,"results","step04","total","annotated_all.tsv.gz"),show_col_types=FALSE)
all_rows<-list()
for(contrast in contrasts) {
  cache<-readRDS(file.path(REPO_ROOT,"results","step05","observed_only","cache",paste0(contrast,".rds")))
  proteins<-cache$prepared |>
    transmute(uniprot_acc=target_id,protein_key,locus_tag,protein=display_name,product,
      interaction_log2_effect=diff,protein_FDR=FDR,ranking_statistic=statistic,
      protein_passes_hit_rule=FDR<.1 & abs(diff)>1,estimate_type,
      n_observed_min_relevant,group_detection_pattern)
  stopifnot(all(proteins$estimate_type=="observed"),!anyDuplicated(proteins$uniprot_acc))
  x<-filter(gsea,.data$contrast==.env$contrast)
  for(i in seq_len(nrow(x))) {
    row<-x[i,]
    ids<-if(is.na(row$leading_edge_ids)||!nzchar(row$leading_edge_ids)) character() else strsplit(row$leading_edge_ids,";",fixed=TRUE)[[1]]
    original<-cache$gsea$leadingEdge[[match(row$pathway,cache$gsea$pathway)]]
    stopifnot(setequal(ids,original),!anyDuplicated(ids),all(ids %in% proteins$uniprot_acc),
      all(ids %in% cache$pathways[[row$pathway]]))
    all_rows[[length(all_rows)+1L]]<-tibble(uniprot_acc=ids,leading_edge_order=seq_along(ids)) |>
      left_join(proteins,by="uniprot_acc",relationship="many-to-one") |>
      mutate(module=row$module,term_id=row$term_id,contrast=contrast,
        contrast_label=row$contrast_label,module_NES=row$NES,module_FDR=row$FDR_study,
        module_significant=row$FDR_study<.1,module_direction=row$direction,
        .before=1)
  }
}
drivers<-bind_rows(all_rows)
# The within-genotype response is SN-TSB; the saved environmental contrast is TSB-SN.
responses<-source_data |> filter(contrast %in% c("TSB_vs_SN_at_anc","TSB_vs_SN_at_stk","TSB_vs_SN_at_stp")) |>
  mutate(genotype=sub("TSB_vs_SN_at_","",contrast),response=-diff) |>
  select(protein_key,genotype,response) |> pivot_wider(names_from=genotype,values_from=response,names_prefix="R_")
drivers<-drivers |> left_join(responses,by="protein_key",relationship="many-to-one") |>
  mutate(expected_interaction=case_when(
    contrast==contrasts[1]~R_anc-R_stk,contrast==contrasts[2]~R_anc-R_stp,
    contrast==contrasts[3]~R_stk-R_stp),score_sign_matches_NES=sign(ranking_statistic)==sign(module_NES)) |>
  arrange(match(term_id,modules$term_id),match(contrast,contrasts),leading_edge_order)
# These supplied product descriptions warrant checking against the broad KEGG membership.
# A flag requests annotation review; it does not silently remove/reassign a member.
review_loci<-c("B7H15_14290","B7H15_11645","B7H15_09250","B7H15_04455","B7H15_07765",
  "B7H15_14945","B7H15_06625","B7H15_06670","B7H15_06675")
drivers<-drivers |> mutate(annotation_review=locus_tag %in% review_loci,
  annotation_note=if_else(annotation_review,"Check supplied KEGG membership against product/reference annotation",NA_character_))
stopifnot(!anyNA(drivers$protein_key),all(drivers$estimate_type=="observed"),
  all(abs(drivers$expected_interaction-drivers$interaction_log2_effect)<1e-8),
  !anyDuplicated(drivers[c("contrast","term_id","uniprot_acc")]))
significant<-filter(drivers,module_significant)
write_csv(drivers,file.path(out,"leading_edge_all_12_comparisons.csv"))
write_tsv(drivers,file.path(out,"leading_edge_all_12_comparisons.tsv"))
write_csv(significant,file.path(out,"leading_edge_significant_modules.csv"))
write_tsv(significant,file.path(out,"leading_edge_significant_modules.tsv"))
overview<-drivers |> group_by(module,term_id,contrast,contrast_label,module_NES,module_FDR,module_significant,module_direction) |>
  summarise(leading_edge_proteins=n(),individual_protein_hits=sum(protein_passes_hit_rule),
    protein_names=paste(protein,collapse=", "),uniprot_accessions=paste(uniprot_acc,collapse=", "),.groups="drop") |>
  arrange(match(term_id,modules$term_id),match(contrast,contrasts))
write_tsv(overview,file.path(out,"module_overview.tsv"))
# Compact one-protein-per-module table: interaction labels identify significant GSEA leading edges.
compact<-significant |> mutate(interaction=case_when(contrast==contrasts[1]~"stk/anc",contrast==contrasts[2]~"stp/anc",TRUE~"stp/stk")) |>
  group_by(module,term_id,protein,locus_tag,uniprot_acc,product) |>
  summarise(leading_edge_in=paste(interaction,collapse=", "),.groups="drop") |>
  arrange(match(term_id,modules$term_id),protein,locus_tag)
write_csv(compact,file.path(out,"protein_list_compact.csv"))
write_tsv(compact,file.path(out,"protein_list_compact.tsv"))
format_effect<-function(effect,hit) paste0(sprintf("%+.2f",effect),if_else(hit,"*",""))
display<-significant |> mutate(interaction=case_when(contrast==contrasts[1]~"stk/anc",contrast==contrasts[2]~"stp/anc",TRUE~"stp/stk"),
  effect_cell=format_effect(interaction_log2_effect,protein_passes_hit_rule)) |>
  select(module,term_id,protein,locus_tag,uniprot_acc,product,annotation_review,interaction,effect_cell) |>
  pivot_wider(names_from=interaction,values_from=effect_cell,values_fill="—") |>
  arrange(match(term_id,modules$term_id),protein,locus_tag)
unknown_labels<-c(B7H15_00135="nucleotidase",B7H15_00805="metallophosphatase / nucleotidase",
  B7H15_04035="deoxyribonucleotidase",B7H15_04455="thermonuclease",B7H15_05185="GTP pyrophosphokinase",
  B7H15_05985="purine NTP pyrophosphatase",B7H15_07765="HU",B7H15_13940="GTP pyrophosphokinase",
  B7H15_14555="ribonucleotide reductase",B7H15_14945="RNA polymerase delta",
  B7H15_05630="ubiquinol oxidase I",B7H15_05635="ubiquinol oxidase II",
  B7H15_06670="peptidase M16",B7H15_06675="peptidase M16",
  B7H15_02470="glutamate synthase",B7H15_13270="nitrate reductase gamma",
  B7H15_13280="nitrate reductase beta",B7H15_13285="nitrate reductase alpha",B7H15_13295="nitrite reductase small subunit")
display<-display |> mutate(label=if_else(protein==locus_tag,paste0(protein," (",unname(unknown_labels[locus_tag]),")"),protein),
  label=if_else(protein=="asnC","asnC (asparagine–tRNA ligase)",label),
  label=paste0(label,if_else(annotation_review,"‡","")))
write_tsv(display,file.path(out,"display_table.tsv"))
md<-c("# Complete leading edges for significant module results", "",
  "Interaction effects are native source log2 difference-of-differences. An asterisk means protein FDR<0.1 and abs(effect)>1. A dash means not in that significant comparison's leading edge, not absent/unmeasured. A double-dagger flags supplied pathway membership for annotation review.","")
for(id in modules$term_id) {
  tab<-filter(display,term_id==id)
  interaction_columns<-if(id=="map00230") c("stk/anc","stp/stk") else if(id=="map00190") "stk/anc" else if(id=="map00910") "stp/anc" else "stp/stk"
  md<-c(md,paste0("## ",tab$module[1]," (",nrow(tab)," proteins)"),"",
    paste0("| Protein | UniProt | ",paste(interaction_columns,collapse=" | ")," |"),
    paste0("| --- | --- | ",paste(rep("---:",length(interaction_columns)),collapse=" | ")," |"))
  for(i in seq_len(nrow(tab))) md<-c(md,paste0("| ",tab$label[i]," | ",tab$uniprot_acc[i]," | ",paste(as.character(tab[i,interaction_columns]),collapse=" | ")," |"))
  md<-c(md,"")
}
writeLines(md,file.path(out,"display_table.md"))
lines<-c("# Proteins driving the four selected KEGG GSEA patterns","",
  "These are all GSEA leading-edge proteins for the observed-only analysis. They are core statistical contributors, not proven causal regulators, direct Stk1/Stp1 substrates, or an exact additive decomposition of NES.","",
  "The significant-module tables contain all leading edges for the five module/interaction pairs passing study FDR<0.1. The all-12 table additionally includes the seven nonsignificant module/interaction pairs, explicitly flagged. Individual proteins need not pass the protein-level hit rule to contribute to GSEA.","",
  "Four modules: purine metabolism, oxidative phosphorylation, nitrogen metabolism and aminoacyl-tRNA biosynthesis. No LOD-imputed source models are included.","",
  "Interaction effect = (A_TSB-B_TSB)-(A_SN-B_SN). With R_g=SN_g-TSB_g, the effects are R_anc-R_stk, R_anc-R_stp and R_stk-R_stp. R values below are supplied model estimates, not flux measurements.","",
  "[Detailed table for significant modules](leading_edge_significant_modules.csv) / [Compact protein list](protein_list_compact.csv) / [All twelve module-comparison pairs](leading_edge_all_12_comparisons.csv) / [Module overview](module_overview.tsv)","",
  "| Module | Interaction | NES | Module FDR | Leading-edge proteins | Individually selected proteins |",
  "| --- | --- | ---: | ---: | ---: | ---: |")
for(i in which(overview$module_significant)) {
  r<-overview[i,]
  lines<-c(lines,paste0("| ",r$module," | ",r$contrast_label," | ",sprintf("%.3f",r$module_NES)," | ",
    format(r$module_FDR,digits=3)," | ",r$leading_edge_proteins," | ",r$individual_protein_hits," |"))
}
lines<-c(lines,"","All listed identifiers and product descriptions come from the supplied mapping/source tables. Reference ortholog annotations do not guarantee an identical JE2 sequence or correct assignment of every broad pathway term. Flagged memberships (for example IsaA/SceD in purine metabolism or RbfA/peptidase M16 in oxidative phosphorylation) require review; they remain in the list to reproduce the reported NES.","",
  "Checks: each leading edge matches the cached fgsea output, every contributor belongs to its tested set and observed-only universe, all IDs are uniquely mapped, and each interaction equals its difference-of-differences reconstructed from the within-genotype responses.")
writeLines(lines,file.path(out,"README.md"))
print(overview |> select(module,contrast_label,module_NES,module_FDR,module_significant,leading_edge_proteins,individual_protein_hits,protein_names),n=12,width=190)
message(nrow(significant)," significant module-protein-comparison rows; ",nrow(compact)," unique module-protein rows.")
