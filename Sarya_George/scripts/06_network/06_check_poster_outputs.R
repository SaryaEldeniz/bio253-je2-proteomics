# Independent source/output reconciliations for the poster networks.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"reports/step06/poster_networks")
fig <- file.path(REPO_ROOT,"figures/step06/poster_networks")
raw <- file.path(REPO_ROOT,"data/raw/step06_string/poster_api")
original <- read_tsv(file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),show_col_types=FALSE)
codes <- c(stk_anc="interaction_stk_vs_anc_at_TSB_vs_SN",stp_anc="interaction_stp_vs_anc_at_TSB_vs_SN",
  stp_stk="interaction_stp_vs_stk_at_TSB_vs_SN")
comparison <- read_csv(file.path(out,"input_set_comparison.csv"),show_col_types=FALSE)
sets <- c("49","67")
if(file.exists(file.path(out,"nodes_75.csv"))) sets <- c(sets,"75")
has75 <- "75" %in% sets
stopifnot(sum(comparison$in_49)==49,sum(comparison$in_67)==67,
  sum(comparison$in_49 & comparison$in_67)==42,
  sum(comparison$in_49 | comparison$in_67)==74,
  nrow(comparison)==if(has75) 75L else 74L)
if(has75) stopifnot("in_75" %in% names(comparison),sum(comparison$in_75)==75L)
network_summary <- read_csv(file.path(out,"network_summary.csv"),show_col_types=FALSE)
live_membership <- read_csv(file.path(out,"pathway_membership_live.csv"),show_col_types=FALSE)
input_files <- c(`49`="user_selected/combined_selected_plus_anchors.txt",
  `67`="all_selected_plus_Stk1_Stp1.txt",`75`="complete75_plus_anchors.txt")
for(set in sets) {
  n <- read_csv(file.path(out,paste0("nodes_",set,".csv")),show_col_types=FALSE)
  live <- read_csv(file.path(out,paste0("nodes_",set,"_live.csv")),show_col_types=FALSE)
  e <- read_csv(file.path(out,paste0("edges_",set,".csv")),show_col_types=FALSE)
  stopifnot(nrow(n)==as.integer(set),nrow(live)==as.integer(set),
    !anyDuplicated(n$string_id),!anyDuplicated(e[c("node_a","node_b")]),
    all(c(e$node_a,e$node_b) %in% n$string_id),all(e$combined_score>=.4),
    setequal(n$string_id,live$string_id),sum(n$node_role=="Biological anchor")==2L)
  # Input, summary and degree checks also catch an isolate lost during joins.
  ids <- readLines(file.path(REPO_ROOT,"reports/step06/string_uploads",input_files[[set]]))
  stopifnot(length(ids)==as.integer(set),!anyDuplicated(ids),setequal(n$uniprot_acc,ids))
  s <- network_summary |> filter(set_id==set)
  stopifnot(nrow(s)==1L,s$nodes==nrow(n),s$edges==nrow(e),
    s$isolates==sum(n$degree==0),s$automatic_added_nodes==0L,
    s$edges_ge_07==sum(e$combined_score>=.7),
    all(n$degree==as.integer(table(factor(c(e$node_a,e$node_b),levels=n$string_id)))))
  if(set=="75") {
    full_ids <- readLines(file.path(REPO_ROOT,"reports/step06/string_uploads/complete75_plus_anchors.txt"))
    held <- c("Q2G010","Q2FYK7","Q2FYG2","Q2FUT0","Q2FV52","Q2FXM9","Q2FWF8","Q2G0W1")
    plot_index <- read_csv(file.path(REPO_ROOT,"reports/step06/sample_plot_index.csv"),show_col_types=FALSE)
    plotted_keys <- unique(unlist(strsplit(plot_index$protein_keys,";",fixed=TRUE)))
    old67 <- read_csv(file.path(out,"nodes_67.csv"),show_col_types=FALSE)
    stopifnot(length(full_ids)==75L,!anyDuplicated(full_ids),
      setequal(n$uniprot_acc,full_ids),all(held %in% n$uniprot_acc),
      setequal(setdiff(n$uniprot_acc,old67$uniprot_acc),held),
      length(plotted_keys)==75L,setequal(n$protein_key,plotted_keys))
  }
  for(code in names(codes)) {
    z <- original |> filter(contrast==codes[[code]])
    j <- match(n$protein_key,z$protein_key)
    stopifnot(!anyNA(j),isTRUE(all.equal(n[[paste0("diff_",code)]],z$diff[j])),
      isTRUE(all.equal(n[[paste0("FDR_",code)]],z$FDR[j])))
    good <- is.finite(n[[paste0("halo_effect_",code)]])
    expected <- z$estimate_type[j]=="observed" & coalesce(z$all_relevant_groups_2of3[j],FALSE) &
      is.finite(z$diff[j])
    stopifnot(all(z$estimate_type[j][good]=="observed"),all(z$n_observed_min_relevant[j][good]>=2),
      identical(good,expected),
      all(n[[paste0("halo_effect_",code)]][good]==z$diff[j][good]),
      isTRUE(all.equal(n[[paste0("halo_FDR_",code)]][good],z$FDR[j][good])))
    # Applying current pathway labels must not change the contrast data.
    k <- match(n$string_id,live$string_id)
    for(field in c("diff_","FDR_","halo_effect_","halo_FDR_"))
      stopifnot(isTRUE(all.equal(n[[paste0(field,code)]],live[[paste0(field,code)]][k])))
  }
  a <- read_tsv(file.path(raw,paste0("functional_annotations_",set,".tsv")),show_col_types=FALSE)
  genes <- strsplit(a$inputGenes,",",fixed=TRUE)
  stopifnot(setequal(unique(unlist(genes)),n$string_id),all(lengths(genes)==a$number_of_genes))
  # Only current STRING KEGG assignments may create pathway outlines.
  m <- live_membership |> filter(set_id==set)
  dictionary <- c(P="map00230",R="map00190",N="map00910",T="map00970")
  for(code in names(dictionary)) {
    rows <- which(a$category=="KEGG" & a$term==dictionary[[code]])
    api_members <- unique(trimws(unlist(genes[rows])))
    display_members <- live$string_id[grepl(paste0("(^|/)",code,"(/|$)"),live$pathway_codes)]
    stopifnot(setequal(api_members,display_members),
      setequal(api_members,m$string_id[m$code==code]))
  }
  for(ext in c("png","pdf","svg")) {
    p <- file.path(fig,paste0("STRING_",set,"_poster_halos.",ext))
    stopifnot(file.exists(p),file.info(p)$size>1000)
  }
}
old <- read_csv(file.path(REPO_ROOT,"reports/step06/string_network/edges_seed_export.csv"),show_col_types=FALSE)
fresh <- read_csv(file.path(out,"edges_49.csv"),show_col_types=FALSE)
joined <- full_join(select(old,node_a,node_b,old=combined_score),
  select(fresh,node_a,node_b,fresh=combined_score),by=c("node_a","node_b"))
stopifnot(nrow(joined)==204L,!anyNA(joined),all(joined$old==joined$fresh))
geometry <- read_csv(file.path(out,"poster_geometry_checks.csv"),show_col_types=FALSE)
stopifnot(nrow(geometry)==length(sets),
  setequal(sub("^STRING_([0-9]+)_poster_halos$","\\1",geometry$set),sets),
  all(geometry$minimum_centre_distance>geometry$outer_node_diameter),
  all(geometry$pinned_anchors))
stopifnot(all(c("anchors_inside_any_outline","minimum_anchor_boundary_margin") %in% names(geometry)),
  all(geometry$anchors_inside_any_outline==0),all(geometry$minimum_anchor_boundary_margin>0))
writeLines(c(if(has75)
  "PASS: exact49/exact67/exact75 inputs; 49/67 common42/union74; all eight additions and all 75 sample-plot proteins retained."
  else "PASS: exact49/exact67 inputs; common42/union74; all nodes and isolates retained.",
  "PASS: fresh49 edge set/scores exactly match original204-edge export.",
  "PASS: original protein effects/FDRs preserved; finite halos are observed and >=2/3 supported.",
  "PASS: current annotation inputGenes scopes/counts match each input set.",
  "PASS: poster exports nonempty; node separation exceeds node diameter; anchors pinned and outside pathway outlines."),
  file.path(out,"validation_checks.txt"))
message("Poster output/source reconciliations passed.")
