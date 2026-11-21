# Convert the existing PDF atlases without recalculating any biological results.
# Inkscape's Poppler importer preserves vector art and embeds existing images.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

export_step09_atlas_svgs <- function() {
  folder <- file.path(REPO_ROOT, "steps", "09_sequence_logos", "reports", "structures_and_logos")
  page_dir <- file.path(folder, "svg_pages")
  qa_dir <- file.path(REPO_ROOT, "results", "step09", "svg_exports")
  for (directory in c(page_dir, qa_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  inkscape <- Sys.getenv("INKSCAPE_BIN", "")
  if (!nzchar(inkscape)) inkscape <- unname(Sys.which("inkscape.com"))
  if (!nzchar(inkscape)) inkscape <- unname(Sys.which("inkscape"))
  if (!nzchar(inkscape) && file.exists("C:/Program Files/Inkscape/bin/inkscape.com")) {
    inkscape <- "C:/Program Files/Inkscape/bin/inkscape.com"
  }
  if (!nzchar(inkscape) || !file.exists(inkscape)) stop("Set INKSCAPE_BIN to an Inkscape 1.4+ executable.")
  if (.Platform$OS.type == "windows") inkscape <- sub("\\.com$", ".exe", inkscape)
  atlases <- tibble(
    base = c("JE2_Phosphosite_Protein_Models", "JE2_Phosphosites_3D_Structures",
             "JE2_Sequence_Logos_Plus_Minus_10_AA"),
    pages = c(4L, 6L, 2L)
  )
  source_pdfs <- file.path(folder, paste0(atlases$base, ".pdf"))
  stopifnot(all(file.exists(source_pdfs)))
  original_hashes <- tools::md5sum(source_pdfs)
  manifest <- list()
  all_outputs <- character()
  for (i in seq_len(nrow(atlases))) {
    base <- atlases$base[i]
    page_files <- character()
    for (page in seq_len(atlases$pages[i])) {
      svg <- file.path(page_dir, sprintf("%s_page%02d.svg", base, page))
      log <- tempfile(sprintf("%s_page%02d_", base, page), tmpdir = qa_dir, fileext = ".log")
      message("Converting ", base, " page ", page, "/", atlases$pages[i])
      reuse <- isTRUE(getOption("step09.svg.reuse", FALSE)) && file.exists(svg)
      if (reuse) {
        status <- 0L
      } else if (.Platform$OS.type == "windows") {
        helper <- file.path(REPO_ROOT, "steps", "09_sequence_logos", "scripts", "phosphosite_figures",
                            "convert_pdf_page_windows.ps1")
        status <- system2("powershell.exe", c("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", shQuote(helper),
          "-Inkscape", shQuote(inkscape), "-InputPdf", shQuote(source_pdfs[i]),
          "-OutputSvg", shQuote(svg), "-Page", page, "-LogPath", shQuote(log)),
          wait = FALSE, stdout = paste0(log, ".launcher.stdout"), stderr = paste0(log, ".launcher.stderr"))
        # Windows graphics programs can keep inherited process handles open.
        # The helper writes its unique completion log only after validating and
        # copying a complete SVG and closing its own Inkscape process.
        started <- Sys.time()
        while (!file.exists(log) && !file.exists(paste0(log, ".failed")) &&
               as.numeric(difftime(Sys.time(), started, units = "secs")) < 55) Sys.sleep(0.2)
        if (file.exists(paste0(log, ".failed"))) stop(paste(readLines(paste0(log, ".failed")), collapse = "\n"))
        if (!file.exists(log)) stop("Conversion did not complete; inspect ", log, ".launcher.stderr")
        status <- 0L
      } else {
        output <- system2(inkscape, c(shQuote(source_pdfs[i]), paste0("--pages=", page),
          "--pdf-poppler", "--export-type=svg", "--export-plain-svg", "--export-text-to-path",
          paste0("--export-filename=", shQuote(svg))), stdout = TRUE, stderr = TRUE, timeout = 60)
        writeLines(output, log, useBytes = TRUE)
        status <- attr(output, "status")
        if (is.null(status)) status <- 0L
      }
      if (status != 0L || !file.exists(svg)) stop("PDF conversion failed; inspect ", log)
      # Prefix all IDs and local references so pages can safely share one SVG.
      text <- readr::read_file(svg)
      prefix <- sprintf("p%02d_", page)
      existing_ids <- xml2::xml_attr(xml2::xml_find_all(xml2::read_xml(svg), "//*[@id]"), "id")
      if (!length(existing_ids) || !all(startsWith(existing_ids, prefix))) {
        text <- gsub('(\\s+id=")([^"]+)(")', paste0("\\1", prefix, "\\2\\3"), text, perl = TRUE)
        text <- gsub('url\\(#([^)]*)\\)', paste0("url(#", prefix, "\\1)"), text, perl = TRUE)
        text <- gsub('((?:xlink:)?href=")#([^"]+)"', paste0("\\1#", prefix, '\\2"'), text, perl = TRUE)
        writeLines(text, svg, useBytes = TRUE)
      }
      document <- xml2::read_xml(svg)
      ids <- xml2::xml_attr(xml2::xml_find_all(document, "//*[@id]"), "id")
      stopifnot(!anyDuplicated(ids))
      images <- xml2::xml_find_all(document, "//*[local-name()='image']")
      paths <- xml2::xml_find_all(document, "//*[local-name()='path']")
      # These two source PDFs are pure vector drawings; do not rasterise them.
      if (base != "JE2_Phosphosites_3D_Structures") stopifnot(length(images) == 0L)
      manifest[[length(manifest) + 1L]] <- tibble(
        source_pdf = basename(source_pdfs[i]), source_pdf_md5 = unname(original_hashes[i]),
        page = page, svg_file = file.path("svg_pages", basename(svg)),
        vector_paths = length(paths), embedded_images = length(images),
        rendering = if (length(images)) "Vector labels with source raster structure images" else "Vector artwork"
      )
      page_files <- c(page_files, svg)
      message("Exported ", basename(svg))
    }
    first <- xml2::read_xml(page_files[1])
    bounds <- as.numeric(strsplit(xml2::xml_attr(xml2::xml_root(first), "viewBox"), "[[:space:]]+")[[1]])
    stopifnot(length(bounds) == 4L, all(is.finite(bounds)), all(bounds[1:2] == 0))
    width <- bounds[3]
    height <- bounds[4]
    # SVG has no PDF-style pages: the same-named atlas stacks pages vertically.
    # The individual page SVGs are easier to crop and arrange on the poster.
    combined <- xml2::read_xml(sprintf(
      '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="%spt" height="%spt" viewBox="0 0 %s %s"/>',
      width, height * length(page_files), width, height * length(page_files)))
    xml2::xml_add_child(combined, "title", base)
    xml2::xml_add_child(combined, "desc", paste("All pages from", basename(source_pdfs[i]), "stacked vertically; original PDF unchanged."))
    for (page in seq_along(page_files)) {
      document <- xml2::read_xml(page_files[page])
      nested <- xml2::xml_root(document)
      xml2::xml_set_attrs(nested, c(x = "0", y = as.character((page - 1) * height),
                                   width = as.character(width), height = as.character(height)))
      xml2::xml_add_child(xml2::xml_root(combined), nested)
    }
    combined_ids <- xml2::xml_attr(xml2::xml_find_all(combined, "//*[@id]"), "id")
    stopifnot(!anyDuplicated(combined_ids))
    atlas_svg <- file.path(folder, paste0(base, ".svg"))
    xml2::write_xml(combined, atlas_svg)
    all_outputs <- c(all_outputs, atlas_svg, page_files)
  }
  stopifnot(identical(original_hashes, tools::md5sum(source_pdfs)))
  manifest <- bind_rows(manifest) |> mutate(source_pdf_unchanged = TRUE)
  write_tsv(manifest, file.path(folder, "svg_conversion_manifest.tsv"))
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(file.path(folder, "JE2_structures_and_logos_SVGs.zip"),
              files = c(substring(all_outputs, nchar(folder) + 2L), "svg_conversion_manifest.tsv"),
              root = folder, mode = "mirror")
  }
  message("Created 12 page SVGs and 3 complete atlas SVGs. Original PDFs are unchanged.")
  invisible(manifest)
}

step09_svg_exports <- export_step09_atlas_svgs()
