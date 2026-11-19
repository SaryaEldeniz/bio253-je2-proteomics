# Archive native STRING images alongside tabular exports, independently of the
# custom poster layout. Website clustering is not substituted for pathway groups.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"reports/step06/poster_networks")
raw <- file.path(REPO_ROOT,"data/raw/step06_string/poster_api")
dir.create(raw,recursive=TRUE,showWarnings=FALSE)
log <- list()
sets <- trimws(strsplit(Sys.getenv("JE2_POSTER_SETS",unset="49,67,75"),",",fixed=TRUE)[[1]])
sets <- sets[file.exists(file.path(out,paste0("nodes_",sets,".csv")))]
for(set in sets) {
  nodes <- read_csv(file.path(out,paste0("nodes_",set,".csv")),show_col_types=FALSE)
  params <- list(identifiers=paste(nodes$string_id,collapse="\r"),species=93061,
    required_score=400,network_type="functional",network_flavor="evidence",
    add_color_nodes=0,add_white_nodes=0,hide_disconnected_nodes=0,
    network_clustering_algorithm="leiden",network_clustering_parameter_leiden=.7,
    cluster_aware_layout=0,network_clustering_edge_visibility="faded",show_cluster_names=0,
    caller_identity="JE2-Project-I-poster")
  encode <- function(x) utils::URLencode(as.character(x),reserved=TRUE)
  url <- paste0("https://version-12-5.string-db.org/api/image/network?",
    paste0(encode(names(params)),"=",vapply(params,encode,character(1)),collapse="&"))
  path <- file.path(raw,paste0("STRING_native_",set,".png"))
  message("Retrieving native STRING ",set,"-node image")
  .api_throttle("STRING")
  handle <- curl::new_handle(timeout=60,connecttimeout=10,followlocation=TRUE,
    useragent="JE2-Project-I-poster/1.0",accept_encoding="")
  r <- tryCatch(curl::curl_fetch_memory(url,handle=handle),error=identity)
  good <- !inherits(r,"error") && r$status_code==200L && length(r$content)>8L &&
    identical(as.integer(r$content[1:8]),c(137L,80L,78L,71L,13L,10L,26L,10L))
  if(good) writeBin(r$content,path)
  log[[set]] <- tibble(set_id=set,url=url,retrieved_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    status=if(inherits(r,"error")) NA_integer_ else r$status_code,
    success=good,file=if(good) basename(path) else NA_character_,
    note=if(inherits(r,"error")) conditionMessage(r) else if(!good) "Not a successful PNG response" else
      "Native STRING evidence view; Leiden resolution0.7; no added nodes; not the custom pathway-group poster")
}
log_path <- file.path(out,"native_image_export_provenance.csv")
old_log <- if(file.exists(log_path)) read_csv(log_path,show_col_types=FALSE) |>
  mutate(set_id=as.character(set_id)) |> filter(!set_id %in% sets) else tibble()
write_csv(bind_rows(old_log,bind_rows(log)),log_path)
files <- list.files(raw,full.names=TRUE)
files <- files[!grepl("[.]zip$",files)]
write_csv(tibble(file=basename(files),bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),
  file.path(out,"raw_export_inventory.csv"))
print(bind_rows(log) |> select(set_id,status,success,file))
