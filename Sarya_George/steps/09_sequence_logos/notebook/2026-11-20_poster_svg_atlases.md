# 2026-11-20: PDF atlas SVG exports and poster consultation

Question: which results deserve space on an A0 poster, and can the three existing sequence/structure PDF atlases be supplied as SVGs?

The source atlases contain 4 linear-map pages, 6 structure pages and 2 sequence-logo pages. Added an R entry script to convert them through Inkscape's Poppler importer, export twelve individual page SVGs, safely prefix SVG IDs and combine each atlas as a vertical SVG. A Windows helper handles hidden graphics-process cleanup after verifying a complete SVG. Added a ZIP bundle and a conversion manifest recording source PDF MD5 hashes, page numbers, vector path counts and embedded image counts. The original PDFs are preserved. No analysis, sequence selection, model retrieval or biological annotation was rerun.

The linear maps and logos remain vector. The source structure PDF contains PyMOL raster renders; those remain embedded images in the converted SVG. Text/page artwork remains vector and typography is outlined for reliable poster import. Local QA renders compare each source PDF page with the corresponding SVG page. QA files remain under ignored `results/step09/svg_exports/`.

Poster consultation distinguishes 14 medium-specific comparison hits from 13 distinct phosphosites on 12 proteins. CsbD S5 occurs in both media. IsaA S194, Pyk T537 and QoxB T552 meet their relevant TSB simple-comparison criteria; none of the nine corresponding interaction estimates meets the FDR threshold. Sequence logos are descriptive frequency summaries of 13 pooled windows or 3-4 windows per group, rather than validated kinase motifs. Structure thumbnails provide mapped reference-model context.

Prepared a poster content plan and a curated supplementary index in Step 10. The index links derived figures, statistical summaries and code; raw course measurements remain local. Publication and QR printing still require checking the final online supplement URL.

AI assistance: Codex reviewed local project evidence, wrote the format-conversion/reproduction code and documentation, and checked the SVG renderings. It did not change the statistical estimates or claim experimental validation of structural effects.
