source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
run_project_step <- function(step) {
  scripts<-c(inventory="01_data/00_inventory_inputs.R",design="01_data/01_validate_design.R",
    identifiers="02_identifiers/03_build_identifier_map.R",api="03_web_resources/03_query_web_resources.R",
    analysis="04_analysis/04_assemble_comparisons.R",reconcile="04_analysis/04_reconcile_counts.R",
    figures="05_figures/04_render_figures.R",enrichment="05_enrichment/05_run_enrichment.R",
    enrichment_figures="05_figures/05_render_enrichment.R",
    enrichment_resources="03_web_resources/05_check_pathway_resources.R",
    enrichment_report="05_enrichment/05_write_report.R",
    enrichment_observed="05_enrichment/05_run_observed_enrichment.R",
    enrichment_all_figures="05_figures/05_render_all_observed.R")
  if(length(step)!=1 || !step %in% names(scripts)) stop("Choose: ",paste(names(scripts),collapse=", "))
  source(file.path(REPO_ROOT,"scripts",scripts[[step]]))
}
if(sys.nframe()==0) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)!=1) stop("Usage: Rscript Sarya_George/run_step.R identifiers (or another step).")
  run_project_step(args[1])
}
