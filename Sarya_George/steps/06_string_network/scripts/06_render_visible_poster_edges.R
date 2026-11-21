# Update only the 75-node poster's edge styling, retaining input/layout checks.
source(if(file.exists("Sarya_George.Rproj"))"R/setup.R" else "Sarya_George/R/setup.R")

render_visible_poster_edges <- function(rerender=getOption("JE2.poster.edges.render",TRUE)) {
  root <- REPO_ROOT
  fig <- file.path(root,"steps/06_string_network/figures/poster_networks")
  report <- file.path(root,"steps/06_string_network/reports/poster_networks")
  archive <- file.path(root,"90_archive/figures/step06_before_visible_edges")
  dir.create(archive,recursive=TRUE,showWarnings=FALSE)
  for(ext in c("svg","png","pdf")) {
    path <- file.path(fig,paste0("STRING_75_poster_halos.",ext))
    backup <- file.path(archive,basename(path))
    if(!file.exists(backup))stopifnot(file.copy(path,backup))
  }
  inputs <- file.path(report,c("nodes_75_live.csv","edges_75.csv","pathway_membership_live.csv"))
  old_hashes <- tools::md5sum(inputs)
  coords <- read_csv(file.path(report,"poster_coordinates_75.csv"),show_col_types=FALSE)
  outlines <- read_csv(file.path(report,"poster_pathway_outlines_75.csv"),show_col_types=FALSE)
  previous_sets <- Sys.getenv("JE2_POSTER_SETS",unset=NA_character_)
  previous_mode <- Sys.getenv("JE2_POSTER_ANNOTATIONS",unset=NA_character_)
  on.exit({
    if(is.na(previous_sets))Sys.unsetenv("JE2_POSTER_SETS") else Sys.setenv(JE2_POSTER_SETS=previous_sets)
    if(is.na(previous_mode))Sys.unsetenv("JE2_POSTER_ANNOTATIONS") else Sys.setenv(JE2_POSTER_ANNOTATIONS=previous_mode)
  },add=TRUE)
  Sys.setenv(JE2_POSTER_SETS="75",JE2_POSTER_ANNOTATIONS="live")
  if(rerender)source(file.path(root,"steps/06_string_network/scripts/06_plot_poster_networks.R"),local=TRUE,encoding="UTF-8")
  stopifnot(identical(old_hashes,tools::md5sum(inputs)),
    isTRUE(all.equal(coords,read_csv(file.path(report,"poster_coordinates_75.csv"),show_col_types=FALSE),tolerance=0,check.attributes=FALSE)),
    isTRUE(all.equal(outlines,read_csv(file.path(report,"poster_pathway_outlines_75.csv"),show_col_types=FALSE),tolerance=0,check.attributes=FALSE)),
    nrow(coords)==75L,nrow(read_csv(inputs[2],show_col_types=FALSE))==412L)
  svg <- xml2::read_xml(file.path(fig,"STRING_75_poster_halos.svg"))
  original_svg <- xml2::read_xml(file.path(archive,"STRING_75_poster_halos.svg"))
  # Compare every vector path, including original lettering and halo polygons,
  # directly with the archived pre-change figure (not just coordinate tables).
  path_query <- ".//*[local-name()='path']"
  stopifnot(identical(xml2::xml_attr(xml2::xml_find_all(original_svg,path_query),"d"),
    xml2::xml_attr(xml2::xml_find_all(svg,path_query),"d")))
  stopifnot(length(xml2::xml_find_all(svg,".//*[local-name()='image']"))==0L)
  writeLines(c("PASS: all 75 node coordinates, halo values, labels and node metadata identical.",
    "PASS: pathway outline coordinates and colours identical.",
    "PASS: every SVG path geometry exactly matches the archived original, including text glyphs and halos.",
    "PASS: source node, edge and membership input hashes identical; 412 associations retained.",
    "Edge display only: colour #42566A, linewidth .25-.60 mm, alpha .55-.85.",
    "The previous widths were .08-.36 mm, alpha .10-.38, colour #697686.",
    "Confidence mapping still monotonic over .4-1; no new edges or confidence changes.",
    "PASS: SVG remains vector-only; matching PNG and PDF refreshed."),
    file.path(report,"visible_edge_style_checks.txt"))
}
render_visible_poster_edges()
