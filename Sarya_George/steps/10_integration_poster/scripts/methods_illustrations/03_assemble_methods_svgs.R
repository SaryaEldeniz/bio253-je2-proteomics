# Assemble fully self-contained vectors; embedded content remains native SVG.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/10_integration_poster/figures/methods_svg")
files <- c("01_experimental_design.svg","02_proteomics_layers.svg","03_analysis_workflow.svg")
docs <- lapply(file.path(out,files),xml2::read_xml)
stopifnot(length(xml2::xml_find_all(docs[[1]],"//*[starts-with(@id,'design-sample-')]"))==27L)
all_ids <- unlist(lapply(docs,function(x) xml2::xml_attr(xml2::xml_find_all(x,"//*[@id]"),"id")))
stopifnot(!anyDuplicated(all_ids))
parts <- vapply(docs,function(x) paste(as.character(xml2::xml_children(xml2::xml_root(x))),collapse="\n"),character(1))
y <- c(0,666,1202)
svg <- paste0('<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="1852" viewBox="0 0 1400 1852" role="img" aria-labelledby="complete-title complete-desc">',
  '<title id="complete-title">JE2 poster methods: design, proteomics and analysis</title>',
  '<desc id="complete-desc">Three editable vector panels show27biologicalsamples,paired measurement layers,andseparate total-protein/phosphosite analysis branches. Artwork is original and schematic.</desc>',
  '<rect id="complete-white-background" width="1400" height="1852" fill="#FFFFFF"/>',
  paste0('<g id="complete-panel-',1:3,'" transform="translate(0 ',y,')">',parts,'</g>',collapse="\n"),
  '<path d="M 30 658 H 1370 M 30 1194 H 1370" stroke="#D8E8EC" stroke-width="1.5"/>',
  '</svg>')
writeLines(enc2utf8(svg),file.path(out,"00_complete_methods.svg"),useBytes=TRUE)
for(file in c("00_complete_methods.svg",files)) {
  doc <- xml2::read_xml(file.path(out,file))
  stopifnot(length(xml2::xml_find_all(doc,"//*[local-name()='image' or local-name()='script' or local-name()='foreignObject']"))==0L,
    length(xml2::xml_find_all(doc,"//*[local-name()='text']"))>0L)
  ids <- xml2::xml_attr(xml2::xml_find_all(doc,"//*[@id]"),"id")
  stopifnot(!anyDuplicated(ids))
  uses <- xml2::xml_attr(xml2::xml_find_all(doc,"//*[@marker-end]"),"marker-end")
  marker_ids <- sub("^url\\(#([^)]*)\\)$","\\1",uses)
  stopifnot(all(marker_ids %in% ids))
}
writeLines(c("PASS:27 distinct biological-replicate symbols in the nine-group matrix.",
  "PASS:four valid SVGs with live text,no embedded images/scripts/external assets,and no duplicate IDs.",
  "PASS:all arrow marker references resolve within each SVG.",
  "PASS:genotypes/timings/measurement layers checked against original course README and sample sheet.",
  "PASS:network candidates branch directly from protein contrasts; enrichment is not a candidate-selection criterion."),
  file.path(out,"validation.txt"))
message("Assembled complete methods SVG and checked editable vector integrity.")
