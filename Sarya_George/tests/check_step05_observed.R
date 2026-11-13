# Read-only verification of the observed-only enrichment and uncapped displays.
# Run from Sarya_George: Rscript tests/check_step05_observed.R
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

validate_step05_observed <- function() {
  report<-file.path(REPO_ROOT,"reports","step05","observed_only")
  output<-file.path(REPO_ROOT,"results","step05","observed_only")
  required<-c(file.path(output,c("ora_all.tsv.gz","gsea_all.tsv.gz")),
    file.path(report,c("universe_audit.tsv","gsea_run_log.tsv")))
  if(!all(file.exists(required))) {
    message("Observed-only analysis/figures are not complete. Validation was not run.")
    message("Missing: ",paste(basename(required[!file.exists(required)]),collapse=", "))
    return(invisible(FALSE))
  }
  log<-c("Observed-only Step 5 verification",format(Sys.time(),"%Y-%m-%d %H:%M %Z"),"")
  on.exit(writeLines(log,file.path(report,"validation.txt")),add=TRUE)
  check<-function(description,condition) {
    if(!isTRUE(condition)) {
      log<<-c(log,paste("FAIL:",description)); stop(description)
    }
    log<<-c(log,paste("PASS:",description)); message("PASS: ",description)
  }
  eq<-function(a,b) isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-10,check.attributes=FALSE))
  split_ids<-function(s) if(is.na(s)||!nzchar(s)) character() else strsplit(s,";",fixed=TRUE)[[1]]
  read<-function(path) read_tsv(path,show_col_types=FALSE,progress=FALSE)
  ora<-read(required[1]); gsea<-read(required[2]); audit<-read(required[3]); run_log<-read(required[4])
  source_data<-read(file.path(REPO_ROOT,"results","step04","total","annotated_all.tsv.gz"))
  plan<-step04_contrast_plan(); collections<-c("KEGG","GO_BP","GO_MF","GO_CC")
  check("Ten complete observed-only runs and single-line warning records",
    nrow(audit)==10 && nrow(run_log)==10 && !anyDuplicated(audit$contrast) &&
      !anyDuplicated(run_log$contrast) && setequal(audit$contrast,plan$contrast) && setequal(run_log$contrast,plan$contrast))
  check("Both methods preserve the observed-only UniProt policy and five-collection testing family",
    all(ora$policy=="observed_only") && all(gsea$policy=="observed_only") &&
      all(ora$namespace=="UniProt") && all(gsea$namespace=="UniProt") &&
      setequal(ora$collection,c(collections,"UniProt_keywords")) && setequal(gsea$collection,c(collections,"UniProt_keywords")))
  check("Unique enrichment tests",
    !anyDuplicated(ora[c("contrast","pathway","effect_sign")]) && !anyDuplicated(gsea[c("contrast","pathway")]))
  for(method in c("ORA","GSEA")) {
    table<-if(method=="ORA") ora else gsea; p<-table[[if(method=="ORA") "p_value" else "pval"]]
    ok<-is.finite(p); expected<-rep(NA_real_,length(p)); expected[ok]<-p.adjust(p[ok],"BH",n=length(p))
    check(paste(method,"study-wide BH independently recalculated across all tests"),eq(table$FDR_study,expected))
    for(contrast in plan$contrast) {
      rows<-which(table$contrast==contrast); q<-rep(NA_real_,length(rows)); finite<-is.finite(p[rows])
      q[finite]<-p.adjust(p[rows][finite],"BH",n=length(rows))
      check(paste(method,contrast,"local BH"),eq(table$FDR_local[rows],q))
    }
  }
  caches<-list()
  for(contrast in plan$contrast) {
    cache<-readRDS(file.path(output,"cache",paste0(contrast,".rds"))); caches[[contrast]]<-cache
    prepared<-cache$prepared
    expected<-source_data |> filter(.data$contrast==.env$contrast,eligible_cds,
      is.finite(diff),is.finite(FDR),is.finite(statistic),eligible_uniprot_annotation,
      n_JE2_good_total_DEA==1,!is.na(uniprot_acc),estimate_type=="observed")
    row<-audit[audit$contrast==contrast,]
    check(paste(contrast,"exact source-observed universe, with no additional replication restriction"),
      !anyDuplicated(prepared$protein_key) && !anyDuplicated(prepared$target_id) &&
        setequal(prepared$protein_key,expected$protein_key) && all(prepared$estimate_type=="observed") &&
        nrow(prepared)==row$universe && row$source_LOD==0 && row$source_observed==nrow(prepared) &&
        row$sparse_observed_retained==sum(!expected$all_relevant_groups_2of3))
    idx<-match(prepared$protein_key,expected$protein_key)
    check(paste(contrast,"source effects/FDR/statistics and full ranks preserved"),
      eq(prepared$diff,expected$diff[idx]) && eq(prepared$FDR,expected$FDR[idx]) &&
        eq(prepared$statistic,expected$statistic[idx]) &&
        setequal(names(cache$ranks),prepared$target_id) &&
        eq(unname(cache$ranks[prepared$target_id]),prepared$statistic))
    check(paste(contrast,"source-observed sparse proteins are retained when present"),
      all(expected$protein_key[!expected$all_relevant_groups_2of3] %in% prepared$protein_key))
    o<-ora[ora$contrast==contrast,]; g<-gsea[gsea$contrast==contrast,]
    foreground<-list(positive=prepared$target_id[prepared$FDR<.1 & prepared$diff>1],
      negative=prepared$target_id[prepared$FDR<.1 & prepared$diff< -1])
    check(paste(contrast,"matched ORA lists and backgrounds"),
      row$positive_hits==length(foreground$positive) && row$negative_hits==length(foreground$negative) &&
        all(o$background_size==nrow(prepared)) && all(g$background_size==nrow(prepared)) &&
        nrow(o)==2*nrow(g) && nrow(g)==row$tested_terms)
    ora_members_ok<-vapply(seq_len(nrow(o)),function(i) {
      ids<-split_ids(o$overlap_ids[i]); pathway<-intersect(cache$pathways[[o$pathway[i]]],prepared$target_id)
      expected_hit<-intersect(foreground[[o$effect_sign[i]]],pathway)
      setequal(ids,expected_hit) && length(ids)==o$overlap[i] && length(pathway)==o$term_size[i]
    },logical(1))
    gsea_members_ok<-vapply(seq_len(nrow(g)),function(i) {
      ids<-split_ids(g$leading_edge_ids[i]); pathway<-intersect(cache$pathways[[g$pathway[i]]],prepared$target_id)
      cached<-match(g$pathway[i],cache$gsea$pathway)
      all(ids %in% pathway) && all(ids %in% prepared$target_id) && length(pathway)==g$size[i] &&
        !is.na(cached) && setequal(ids,cache$gsea$leadingEdge[[cached]]) &&
        eq(g$NES[i],cache$gsea$NES[cached]) && eq(g$pval[i],cache$gsea$pval[cached])
    },logical(1))
    check(paste(contrast,"all overlaps and leading edges consist of eligible observed proteins"),
      all(ora_members_ok) && all(gsea_members_ok))
  }
  chosen<-ora |> filter(is.finite(FDR_study),FDR_study<.1) |> arrange(p_value) |> slice_head(n=1)
  check("A real significant ORA row is available for independent verification",nrow(chosen)==1)
  N<-chosen$background_size; n<-chosen$foreground_size; K<-chosen$term_size; k<-chosen$overlap
  contingency<-matrix(c(k,n-k,K-k,N-K-n+k),nrow=2,byrow=TRUE)
  check("Independent one-sided Fisher exact test matches a real ORA result",
    all(contingency>=0) && eq(chosen$p_value,fisher.test(contingency,alternative="greater")$p.value) &&
      eq(chosen$fold_enrichment,(k/n)/(K/N)))
  figure_paths<-file.path(report,c("figure_manifest.tsv","all_plotted_significant.tsv"))
  if(!all(file.exists(figure_paths))) {
    log<-c(log,"","Numerical checks passed. Figure coverage checks await completed rendering.")
    message(tail(log,1)); return(invisible(FALSE))
  }
  manifest<-read(figure_paths[1]); plotted<-read(figure_paths[2])
  expected_groups<-expand_grid(method=c("ORA","GSEA"),collection=collections,contrast=plan$contrast)
  group_key<-function(x) paste(x$method,x$collection,x$contrast,sep="::")
  check("All 80 plot groups are documented, including empty groups",
    nrow(manifest)==80 && !anyDuplicated(group_key(manifest)) && setequal(group_key(manifest),group_key(expected_groups)))
  check("Every plot declares observed-only policy and no top-N truncation",
    all(manifest$policy=="observed_only") && all(manifest$top_limit=="none"))
  expected_plotted<-bind_rows(list(ORA=ora,GSEA=gsea),.id="method") |>
    filter(collection %in% collections,is.finite(FDR_study),FDR_study<.1)
  result_key<-function(x) paste(x$method,x$contrast,x$collection,x$pathway,x$effect_sign,sep="::")
  check("Combined plotted-row table contains every significant GO/KEGG result exactly once",
    !anyDuplicated(result_key(plotted)) && setequal(result_key(plotted),result_key(expected_plotted)))
  resolve<-function(path) if(grepl("^([A-Za-z]:|/)",path)) path else file.path(REPO_ROOT,path)
  for(i in seq_len(nrow(manifest))) {
    m<-manifest[i,]; table<-read(resolve(m$table))
    expected<-expected_plotted |> filter(method==m$method,collection==m$collection,contrast==m$contrast)
    coverage_ok<-nrow(table)==m$significant_rows && setequal(result_key(table),result_key(expected))
    if(m$method=="GSEA") {
      order_ok<-all(diff(table$NES)<=0) && identical(as.integer(table$draw_order),seq_len(nrow(table)))
    } else {
      order_ok<-all(vapply(split(table,table$effect_sign),function(part) {
        all(diff(part$fold_enrichment)<=0) && identical(as.integer(part$draw_order),seq_len(nrow(part)))
      },logical(1)))
    }
    check(paste(m$method,m$collection,m$contrast,"complete membership, descending rows and PNG/PDF files"),
      coverage_ok && order_ok && file.exists(resolve(m$png)) && file.exists(resolve(m$pdf)))
  }
  log<-c(log,"","All observed-only analysis and uncapped plotting checks passed.")
  message(tail(log,1))
  invisible(TRUE)
}
validate_step05_observed()
