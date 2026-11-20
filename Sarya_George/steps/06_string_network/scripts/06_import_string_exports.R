# Import the student's 49-node STRING 12.5 exports; preserve raw downloads locally.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/06_string_network/reports/string_network")
local <- file.path(REPO_ROOT,"results/step06/string_network")
raw <- file.path(REPO_ROOT,"data/raw/step06_string/seed_export")
for(p in c(out,local,raw)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
download_dir <- Sys.getenv("JE2_STRING_DOWNLOAD_DIR",unset=file.path(Sys.getenv("USERPROFILE"),"Downloads"))
files <- c("string_interactions_short_49_nodes.tsv","string_functional_annotations_49_nodes.tsv",
  "string_protein_annotations_49_nodes.tsv","string_protein_sequences_49_nodes.fa",
  "string_coordinates_49_nodes.tsv","string_node_degrees_49_nodes.tsv","string_interactions_psimi_49_nodes.xml")
for(f in files) {
  source_file <- file.path(download_dir,f)
  if(!file.exists(source_file) && !file.exists(file.path(raw,f))) stop("Missing STRING input: ",f)
  if(file.exists(source_file)) stopifnot(file.copy(source_file,file.path(raw,f),overwrite=TRUE))
}
read_export <- function(f) read_tsv(file.path(raw,f),show_col_types=FALSE,name_repair="check_unique")
proteins <- read_export(files[3])
# This export has URL/description headings swapped. Detect content, not assumptions.
swapped <- all(grepl("^https?://",proteins$annotation)) &&
  all(!grepl("^https?://",proteins$domain_summary_url))
if(swapped) {
  temp <- proteins$annotation
  proteins$annotation <- proteins$domain_summary_url
  proteins$domain_summary_url <- temp
}
aliases <- proteins |> select(string_id=identifier,aliases=other_names_and_aliases) |>
  mutate(alias=strsplit(aliases,",",fixed=TRUE)) |> unnest(alias) |>
  mutate(alias=trimws(alias)) |> distinct(string_id,alias)
previous <- read_csv(file.path(REPO_ROOT,"steps/06_string_network/reports/string_uploads/user_selected/network_node_attributes.csv"),show_col_types=FALSE)
mapping <- previous |> select(uniprot_acc,protein_key,locus_tag,display_name,upload_role) |>
  inner_join(aliases,by=c("uniprot_acc"="alias"),relationship="one-to-one")
stopifnot(nrow(mapping)==49L,!anyDuplicated(mapping$string_id),!anyDuplicated(mapping$uniprot_acc),
  setequal(mapping$uniprot_acc,previous$uniprot_acc),all(startsWith(mapping$string_id,"93061.")))
seeds <- previous |> left_join(select(mapping,uniprot_acc,string_id),by="uniprot_acc",relationship="one-to-one") |>
  left_join(proteins |> transmute(string_id=identifier,string_name=.data[["#node"]],
    string_description=annotation,domain_summary_url),by="string_id",relationship="one-to-one") |>
  mutate(node_role=if_else(upload_role=="Biological anchor","Biological anchor","Selected seed"),
    module_group=if_else(node_role=="Biological anchor","Regulatory anchors",modules),
    measured_in_total=TRUE,mapping_status="Exact UniProt alias token in student STRING export")
edges_raw <- read_export(files[1])
edges <- edges_raw |> transmute(node_a=pmin(node1_string_id,node2_string_id),
  node_b=pmax(node1_string_id,node2_string_id),combined_score,
  neighborhood=neighborhood_on_chromosome,fusion=gene_fusion,cooccurrence=phylogenetic_cooccurrence,
  coexpression,experimental=experimentally_determined_interaction,database=database_annotated,
  textmining=automated_textmining,edge_origin="Student export") |> distinct()
stopifnot(!anyDuplicated(edges[c("node_a","node_b")]),all(edges$node_a!=edges$node_b),
  all(c(edges$node_a,edges$node_b) %in% seeds$string_id),all(edges$combined_score>=.4 & edges$combined_score<=1))
coordinates <- read_export(files[5])
degrees <- read_export(files[6])
sequence_ids <- readLines(file.path(raw,files[4]),warn=FALSE)
sequence_ids <- sub(" .*","",substring(sequence_ids[startsWith(sequence_ids,">")],2))
stopifnot(setequal(coordinates$identifier,seeds$string_id),setequal(degrees$identifier,seeds$string_id),
  setequal(sequence_ids,seeds$string_id))
calc_degree <- function(e,ids) as.integer(table(factor(c(e$node_a,e$node_b),levels=ids)))
seeds$degree_seed_export <- calc_degree(edges,seeds$string_id)
seeds$degree_seed_high_confidence <- calc_degree(filter(edges,combined_score>=.7),seeds$string_id)
stopifnot(all(seeds$degree_seed_export==degrees$node_degree[match(seeds$string_id,degrees$identifier)]))
annotations <- read_export(files[2])
stopifnot(all(annotations$identifier %in% seeds$string_id))
annotation_counts <- annotations |> group_by(category) |>
  summarise(node_term_rows=n(),unique_terms=n_distinct(.data[["term ID"]]),
    annotated_nodes=n_distinct(identifier),.groups="drop")
write_csv(mapping,file.path(out,"seed_identifier_mapping.csv"))
write_csv(seeds,file.path(out,"nodes_seed.csv"))
write_csv(edges,file.path(out,"edges_seed_export.csv"))
write_csv(filter(edges,combined_score>=.7),file.path(out,"edges_seed_high_confidence.csv"))
write_csv(annotation_counts,file.path(out,"functional_annotation_inventory.csv"))
write_csv(tibble(file=files,md5=unname(tools::md5sum(file.path(raw,files)))),file.path(out,"export_input_hashes.csv"))
write_csv(tibble(key=c("STRING_version","species","network_type","main_threshold","sensitivity_threshold",
  "expansion_partners_per_seed","max_added_nodes","export_annotation_headers_swapped",
  "seed_export_network_type_UI","seed_export_edge_display","seed_export_evidence_sources",
  "seed_export_first_shell","seed_export_second_shell","seed_export_layout",
  "seed_export_clustering","seed_export_clustering_resolution","seed_export_use_clusters_in_layout",
  "seed_export_cluster_edge_display"),
  value=c("12.5 (student declaration)","93061","functional","0.4","0.7","5","12",as.character(swapped),
    "Full STRING network","Evidence","All shown sources enabled including evidence transfer",
    "none / query proteins only","none","Force-directed","Leiden","0.7","FALSE","faded lines")),
  file.path(out,"network_settings.csv"))
saveRDS(list(seeds=seeds,edges=edges,aliases=aliases,annotations=annotations,coordinates=coordinates),
  file.path(local,"seed_import.rds"))
writeLines(c("All 49 requested IDs uniquely matched through exact alias tokens; all in taxon93061.",
  paste0(nrow(edges)," edges at >=0.4; ",sum(edges$combined_score>=.7)," at >=0.7."),
  paste0(sum(seeds$degree_seed_export==0)," isolates at >=0.4; ",sum(seeds$degree_seed_high_confidence==0)," at >=0.7."),
  paste0("Annotation description/URL columns corrected: ",swapped),
  "Functional annotation rows are memberships, not ORA/GSEA results; no new enrichment FDR is present.",
  "Reference FASTA sequences are not JE2 phosphosite-coordinate sequences.",
  "PSI-MI XML was archived as supporting evidence, not treated as a STRING version declaration."),
  file.path(out,"import_checks.txt"))
message("Imported49nodes; ",nrow(edges)," edges; alias matches/degree/FASTA checks passed.")
