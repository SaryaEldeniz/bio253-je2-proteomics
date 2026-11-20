# Restore all eight manually held candidates to the exploratory poster network.
# A manual question about pathway membership is not evidence against inclusion.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads")
old <- readLines(file.path(out,"all_selected_plus_Stk1_Stp1.txt"))
added_ids <- c("Q2G010","Q2FYK7","Q2FYG2","Q2FUT0","Q2FV52","Q2FXM9","Q2FWF8","Q2G0W1")
decisions <- read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/04_candidate_decision_table.csv"),show_col_types=FALSE)
added <- decisions |> filter(uniprot_acc %in% added_ids) |>
  mutate(review_origin="Analyst-authored review note; not an official STRING/lecturer warning",
    inclusion_reason="Restored for exploratory network analysis; manual pathway question does not exclude protein")
stopifnot(nrow(added)==8L,!anyDuplicated(added$uniprot_acc),
  all(added$annotation_review),all(added$n_supported_interactions>0),
  !any(added$uniprot_acc %in% old))
ids <- c(old,added_ids)
plot_index <- read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/sample_plot_index.csv"),show_col_types=FALSE)
plotted_keys <- unique(unlist(strsplit(plot_index$protein_keys,";",fixed=TRUE)))
identity <- bind_rows(
  read_csv(file.path(out,"upload_manifest.csv"),show_col_types=FALSE) |> select(uniprot_acc,protein_key),
  added |> select(uniprot_acc,protein_key)) |> distinct()
stopifnot(length(ids)==75L,!anyDuplicated(ids),nrow(identity)==75L,
  setequal(identity$protein_key,plotted_keys),setequal(ids,identity$uniprot_acc))
writeLines(ids,file.path(out,"complete75_plus_anchors.txt"))
write_csv(added,file.path(out,"complete75_restored_eight.csv"))
message("Complete 75-protein input matches the sample-plot inventory exactly; all eight restored.")
