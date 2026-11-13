# Rebuild the readable results guide from saved analyses; no tests are rerun.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out<-file.path(REPO_ROOT,"reports","step05")
local<-file.path(REPO_ROOT,"results","step05")
ora<-read_tsv(file.path(local,"ora_all.tsv.gz"),show_col_types=FALSE)
gsea<-read_tsv(file.path(local,"gsea_all.tsv.gz"),show_col_types=FALSE)
audit<-read_tsv(file.path(out,"universe_and_hit_counts.tsv"),show_col_types=FALSE)
coverage<-read_tsv(file.path(out,"annotation_coverage.tsv"),show_col_types=FALSE)
counts<-read_tsv(file.path(out,"pathway_counts.tsv"),show_col_types=FALSE)
sens<-read_tsv(file.path(local,"gsea_sensitivity.tsv.gz"),show_col_types=FALSE)
sens<-sens |> mutate(
  significant_replicated=FDR_study_replicated<.1,
  significant_observed_replicated=FDR_study_observed_replicated<.1,
  significant_local_JE2=FDR_study_local_JE2<.1,
  same_sign_all=sign(NES_primary)==sign(NES_replicated) & sign(NES_primary)==sign(NES_observed_replicated) & sign(NES_primary)==sign(NES_local_JE2),
  significant_same_sign_all=FDR_study_primary<.1 & significant_replicated & significant_observed_replicated & significant_local_JE2 & same_sign_all)
write_tsv(filter(sens,FDR_study_primary<.1),file.path(out,"gsea_primary_sensitivity.tsv"))
gsea_primary<-filter(gsea,policy=="primary") |>
  left_join(select(sens,contrast,pathway,significant_replicated,significant_observed_replicated,
    significant_local_JE2,same_sign_all,significant_same_sign_all),by=c("contrast","pathway"),relationship="one-to-one")
ora_primary<-filter(ora,policy=="primary")
evidence<-full_join(ora_primary |> select(contrast,pathway,collection,term_description,effect_sign,
  ORA_FDR=FDR_study,fold_enrichment,overlap,overlap_names),
  gsea_primary |> select(contrast,pathway,effect_sign,GSEA_FDR=FDR_study,NES,leading_edge_names,significant_same_sign_all),
  by=c("contrast","pathway","effect_sign"),relationship="one-to-one") |>
  filter(ORA_FDR<.1 | GSEA_FDR<.1) |> arrange(contrast,GSEA_FDR,ORA_FDR)
write_tsv(evidence,file.path(out,"primary_method_comparison.tsv"))
md_table<-function(x) {
  if(!nrow(x)) return("No categories meet the declared threshold.")
  x<-mutate(x,across(where(is.numeric),function(v) {
    finite<-v[is.finite(v)]
    if(length(finite) && all(finite==floor(finite))) format(v,scientific=FALSE,trim=TRUE) else format(signif(v,3),trim=TRUE)
  }))
  cells<-lapply(x,function(v) gsub("[|\r\n]"," / ",ifelse(is.na(v),"NA",as.character(v))))
  rows<-do.call(cbind,cells)
  c(paste0("| ",paste(names(x),collapse=" | ")," |"),paste0("| ",paste(rep("---",ncol(x)),collapse=" | ")," |"),
    apply(rows,1,function(row) paste0("| ",paste(row,collapse=" | ")," |")))
}
plan<-step04_contrast_plan()
focus<-c("TSB_vs_SN_at_anc",filter(plan,contrast_type=="interaction")$contrast)
totals<-counts |> filter(policy=="primary") |> group_by(method) |>
  summarise(tests=sum(tests),significant_term_comparisons=sum(significant_study),.groups="drop")
focus_counts<-counts |> filter(policy=="primary",contrast %in% focus) |>
  select(contrast_label,method,collection,significant_study) |>
  pivot_wider(names_from=collection,values_from=significant_study) |> arrange(contrast_label,method)
top_kegg<-gsea_primary |> filter(contrast %in% focus,collection=="KEGG",FDR_study<.1) |>
  group_by(contrast) |> arrange(FDR_study,pathway,.by_group=TRUE) |> slice_head(n=5) |> ungroup() |>
  select(contrast_label,term_description,NES,FDR_study,direction,significant_same_sign_all)
robustness<-gsea_primary |> filter(FDR_study<.1) |> group_by(contrast_label) |>
  summarise(primary_significant=n(),replicated_significant=sum(significant_replicated,na.rm=TRUE),
    observed_replicated_significant=sum(significant_observed_replicated,na.rm=TRUE),
    local_JE2_significant=sum(significant_local_JE2,na.rm=TRUE),
    significant_same_sign_all=sum(significant_same_sign_all,na.rm=TRUE),.groups="drop")
coverage_primary<-coverage |> filter(policy=="primary",contrast=="TSB_vs_SN_at_anc") |>
  select(collection,annotated_proteins,background_size,coverage_fraction,terms_tested)
resource_note<-"Public-resource checks have not been run; run enrichment_resources separately."
if(file.exists(file.path(out,"live_uniprot_checks.tsv"))) {
  live_up<-read_tsv(file.path(out,"live_uniprot_checks.tsv"),show_col_types=FALSE)
  live_kegg<-read_tsv(file.path(out,"live_kegg_pathway_checks.tsv"),show_col_types=FALSE)
  live_string<-read_tsv(file.path(out,"live_string_mapping_review.tsv"),show_col_types=FALSE)
  resource_note<-paste0(nrow(live_up)," UniProt entries were queried; ",sum(is.na(live_up$lookup_error)),
    " returned without lookup errors. ",sum(live_kegg$status=="retrieved")," KEGG maps were retrieved. ",
    sum(live_string$accepted_for_context_network %in% TRUE),
    " STRING reference IDs were verified against UniProt cross-references. The small functional network provides context; it is not independent proof of direct targets.")
}
lines<-c("# Step 5: October total-proteome functional enrichment", "",
  "> Updated display request (2026-11-13): [all significant GO/KEGG results excluding source LOD models](observed_only/README.md). The report below retains the original primary analysis; the new observed-only figures use their own matching background and recalculated enrichment statistics. Sample-evidence PNGs now show TSB and SN only.","",
  paste0("Generated ",format(Sys.time(),"%Y-%m-%d %H:%M %Z")," from the saved October comparisons."),"",
  "ORA and GSEA are complete for ten comparisons, with three separately labelled sensitivity policies. The source protein effects/FDRs were preserved. These analyses identify concentrations of protein abundance changes in functional categories; they do not directly measure pathway activity or phosphorylation.","",
  "## How many tests and findings?", "",md_table(totals),"",
  "Counts are term-comparison results, not independent pathways: GO categories and KEGG maps overlap. Main significance is study-wide BH FDR <0.1, separately for ORA and GSEA across all ten comparisons and five collections. Sensitivity policies have their own correction families.","",
  "## Your ancestor response and three interactions", "",md_table(focus_counts),"",
  "### Leading KEGG findings", "",md_table(top_kegg),"",
  "The table displays the five lowest study FDRs per focus comparison. `significant_same_sign_all` means the term remains significant with the same NES direction in all three sensitivity policies. A false value does not automatically invalidate a result; inspect the exact sensitivity values and changed background sizes. NA means a corresponding sensitivity test/result is unavailable, not a negative result.","",
  "### Check sensitivity before naming a module", "",md_table(robustness),"",
  "`replicated` requires at least 2/3 measurements in every relevant group; `observed_replicated` also requires the source model's observed label; `local_JE2` restores eligible entries excluded by reference mapping. Loss of significance can reflect changed evidence, power, term membership or the background; it is not proof of biological absence.","",
  "## The annotation coverage behind the tests", "",md_table(coverage_primary),"",
  "This is the primary ancestor universe. The complete annotation_coverage.tsv contains each comparison/policy. The common tested-protein background includes proteins without an annotation in a collection; coverage makes this choice visible. Gene sets are intersected with that background before applying the 10-500 member size filter. Interpret broad map names through their contributing proteins, rather than assuming every named process occurs in JE2.","",
  "## What was done, in order", "",
  "1. Read the existing October Step 4 total-proteome table and the Step 3 identifier map/term snapshot. Do not refit source models or replace protein-level FDRs.",
  "2. Remove non-CDS entries and require finite source statistics. For the primary analysis retain high-quality singleton JE2-to-UniProt mappings, so one target accession is counted once. Keep an exclusion audit.",
  "3. Transfer supplied STRING-v12 KEGG, GO BP/MF/CC and UniProt-keyword memberships through those mappings. Deduplicate member-term pairs. Freeze this snapshot and hash its inputs.",
  "4. Define each comparison's entire eligible tested-protein background. Create separate positive/negative ORA lists using source FDR <0.1 and absolute log2 effect >1. Test every eligible term, including zero overlaps, with the upper-tail hypergeometric test.",
  "5. Rank all eligible proteins by the signed source DEA statistic. Run fgseaMultilevel with the declared seeds and serial backend; keep ES, NES, p-value, estimation uncertainty, warnings and leading-edge proteins. The ranking is not restricted to significant proteins.",
  "6. Correct enrichment p-values using the declared local and study-wide BH families. ORA includes both directions in each family; it is adjusted separately from GSEA. The engine's original adjusted values remain separately available.",
  "7. Repeat with the replicated, observed_replicated and local_JE2 policies, always changing foreground, rank vector and background together. Compare each primary finding with its corresponding sensitivity results.",
  "8. Export full pathway tables, directional overlaps, leading edges, sensitivity comparisons and candidate contributors. Plot KEGG/GO findings and NES heatmaps with explicit comparison directions.",
  "9. Query small public-ID lists through the existing UniProt, KEGG and STRING clients for interpretation. Check reference STRING matches against UniProt cross-references. Keep live lookups separate from frozen memberships. Plot measured sample values for selected contributors before interpreting individual proteins.",
  "10. Save the scripts, methods, source hashes, run warnings, package versions, figures and a dated lab notebook entry. Raw inputs and full local caches are not committed.","",
  "## Public-resource checks", "",resource_note,"",
  "See [lookup policy](resource_checks.md), [UniProt entries](live_uniprot_checks.tsv), [STRING mapping review](live_string_mapping_review.tsv), [STRING functional associations](live_string_context_edges.tsv) and [KEGG map checks](live_kegg_pathway_checks.tsv). The API scripts preserve dated responses in local caches. Public service annotations can differ from the frozen course snapshot.","",
  "## Read the direction correctly", "",
  "For TSB_vs_SN, positive NES concentrates toward higher TSB protein values; negative NES toward higher SN values. For interactions define R_g = log2(SN_g)-log2(TSB_g). Then stk_vs_anc is R_anc-R_stk, stp_vs_anc is R_anc-R_stp, and stp_vs_stk is R_stk-R_stp. More positive can mean a larger increase or a smaller decrease. Enrichment of an interaction is not a direct test of pathway activation.","",
  "## Files to open", "",
  "- [Biological interpretation](interpretation.md) discusses the strongest modules and the limitations of narrower candidates.",
  "- [Full methods guide](../../docs/step05_methods.md) explains the formulas, policies and interpretation.",
  "- [ORA primary significant table](ora_primary_significant.tsv) gives enrichment factors and contributing hit IDs/names.",
  "- [GSEA primary significant table](gsea_primary_significant.tsv) gives NES, uncertainty and leading edges.",
  "- [ORA/GSEA comparison](primary_method_comparison.tsv) keeps their evidence separate; it does not combine p-values.",
  "- [GSEA sensitivity](gsea_primary_sensitivity.tsv) and [ORA sensitivity](ora_primary_sensitivity.tsv) show exact policy-specific results.",
  "- [Candidate module contributors](candidate_module_contributors.tsv) are follow-up candidates, not proven kinase substrates.",
  "- [KEGG contributing proteins](kegg_contributing_proteins.tsv) connects terms to protein IDs and annotations.",
  "- [Universe/hit counts](universe_and_hit_counts.tsv), [coverage](annotation_coverage.tsv), [warnings](gsea_run_log.tsv), [input hashes](input_hashes.tsv), and [settings](run_settings.txt) retain the audit trail.",
  "- [Completed result checks](result_validation.txt) independently verify correction families, backgrounds, ranks, signs, leading edges and a real ORA result.",
  "- Full local outputs are results/step05/ora_all.tsv.gz and gsea_all.tsv.gz; per-comparison TSVs, exact ranks and caches are in results/step05/contrasts and cache. Open these in RStudio; they remain ignored by Git.","",
  "## Figures", "",
  "![KEGG across all ten contrasts](../../figures/step05/step05_KEGG_all10_NES_heatmap.png)","",
  "![KEGG interaction GSEA](../../figures/step05/step05_KEGG_interactions_GSEA.png)","",
  "![GO BP interaction GSEA](../../figures/step05/step05_GO_BP_interactions_GSEA.png)","",
  "![Measured sample values](../../figures/step05/step05_candidate_sample_values.png)","",
  "![Nitrogen module sample values](../../figures/step05/step05_nitrogen_sample_values.png)","",
  "## Important limits", "",
  "The local annotation snapshot is STRING v12 (custom organism STRG0A35SHG); its exact minor release and underlying GO/KEGG dates are not supplied. Slides show a separate STRING12.5 example organism. Live reference-organism lookups are explicitly separate and do not silently refresh memberships.","",
  "The supplied observed/LOD label is model-wide; a flagged model may have complete measurements in this particular contrast. Conversely observed does not guarantee all 3/3 values. Missing measurements remain NA. The sample plots and replication sensitivities expose support.","",
  "These are protein-abundance analyses with n=3 biological replicates per group. Preranked GSEA does not permute sample labels. Neither ORA nor GSEA proves flux, direct physical binding, phosphorylation dependence or direct Stk1 substrates. The reduced Stk1 abundance in the stp strain complicates pure phosphatase-loss interpretations. The phosphosite/Venn/sequence-logo work follows later after site localization and parent-protein checks.","",
  "## RStudio commands", "", "```r", 'source("run_step.R")',
  'run_project_step("enrichment")', 'run_project_step("enrichment_figures")',
  'run_project_step("enrichment_resources")  # small public-ID queries; network required',
  'run_project_step("enrichment_report")',
  'View(read_enrichment("GSEA", contrast="interaction_stk_vs_anc_at_TSB_vs_SN", collection="KEGG"))',
  'View(read_enrichment("ORA", significant_only=FALSE))',"```","",
  "## Method references", "",
  "- [Official fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html).",
  "- [R hypergeometric test](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html).",
  "- [R multiple-testing correction](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).",
  "- Course 06_web_resources_2026.pdf: resource selection; local memberships and ortholog mapping are the supplied Week 2 snapshot.")
writeLines(lines,file.path(out,"README.md"))
message("Readable Step 5 report rebuilt: reports/step05/README.md")
