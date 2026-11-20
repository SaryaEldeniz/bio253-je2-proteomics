source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
run_project_step <- function(step) {
  scripts<-c(inventory="steps/01_design/scripts/00_inventory_inputs.R",design="steps/01_design/scripts/01_validate_design.R",
    identifiers="steps/03_identifiers_resources/scripts/03_build_identifier_map.R",api="steps/03_identifiers_resources/scripts/03_query_web_resources.R",
    analysis="steps/04_total_proteome/scripts/04_assemble_comparisons.R",reconcile="steps/04_total_proteome/scripts/04_reconcile_counts.R",
    figures="steps/04_total_proteome/scripts/04_render_figures.R",enrichment="steps/05_enrichment/scripts/05_run_enrichment.R",
    enrichment_figures="steps/05_enrichment/scripts/05_render_enrichment.R",
    enrichment_resources="steps/05_enrichment/scripts/05_check_pathway_resources.R",
    enrichment_report="steps/05_enrichment/scripts/05_write_report.R",
    enrichment_observed="steps/05_enrichment/scripts/05_run_observed_enrichment.R",
    enrichment_all_figures="steps/05_enrichment/scripts/05_render_all_observed.R")
  if(length(step)!=1 || !step %in% names(scripts)) stop("Choose: ",paste(names(scripts),collapse=", "))
  source(file.path(REPO_ROOT,scripts[[step]]))
}
if(sys.nframe()==0) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)!=1) stop("Usage: Rscript Sarya_George/run_step.R identifiers (or another step).")
  run_project_step(args[1])
}
