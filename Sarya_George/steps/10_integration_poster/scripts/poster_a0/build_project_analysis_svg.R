# Poster adaptation of 03_analysis_workflow.svg; original artwork is preserved.
# At 529.33 mm wide, 24 SVG units print at 25.7 pt (above the lecture minimum).
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
ink <- "#17394D"; blue <- "#247BA0"; teal <- "#078D92"; muted <- "#52717C"
out <- file.path(REPO_ROOT, "steps/10_integration_poster/figures/methods_svg/03_analysis_workflow_poster.svg")
esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); gsub("<", "&lt;", x, fixed = TRUE) }
label <- function(x, y, text, size = 24, weight = 400, color = ink, anchor = "start") {
  sprintf('<text x="%s" y="%s" font-family="Arial" font-size="%s" font-weight="%s" fill="%s" text-anchor="%s">%s</text>',
          x, y, size, weight, color, anchor, esc(text))
}
box <- function(x,y,w,h,fill) sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="16" fill="%s"/>',x,y,w,h,fill)
line <- function(d,color=teal,width=2.5,arrow=FALSE) {
  path <- sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linecap="round"/>',d,color,width)
  if (!arrow) return(path)
  # Explicit triangle arrowheads survive PowerPoint SVG import more reliably
  # than SVG marker references. The analysis connections are unchanged.
  xy <- as.numeric(regmatches(d,gregexpr("[0-9]+",d))[[1]])
  segment <- tail(xy,4); start <- segment[1:2]; end <- segment[3:4]
  direction <- (end-start)/sqrt(sum((end-start)^2)); normal <- c(-direction[2],direction[1])
  a <- end-10*direction+5*normal; b <- end-10*direction-5*normal
  paste0(path,sprintf('<polygon points="%s,%s %s,%s %s,%s" fill="%s"/>',end[1],end[2],a[1],a[2],b[1],b[2],color))
}
circle <- function(x,y,r,fill,stroke=blue) sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s" stroke="%s" stroke-width="2"/>',x,y,r,fill,stroke)
z <- c('<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="460" viewBox="0 0 1400 460" role="img" aria-labelledby="title desc">',
 '<title id="title">Project analysis</title>',
 '<desc id="desc">Poster adaptation of the existing project-analysis diagram. Identifier mapping feeds parallel total-protein and phosphosite analyses. STRING candidates come directly from protein contrasts, independently of enrichment NES. Venn sets use within-medium mutant comparisons; phosphosite interactions are tested separately. Integration prioritises candidates without claiming direct substrates or phosphorylation occupancy.</desc>',
 '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10Z" fill="#078D92"/></marker></defs>',
 label(20,34,"Project analysis",34,700),label(1380,34,"TSB vs long-exposure SN",24,400,muted,"end"),
 box(20,80,245,325,"#EAF0F7"),
 box(310,80,310,140,"#EDF7F8"),box(665,80,300,140,"#EDF7F8"),box(1010,80,370,140,"#EDF7F8"),
 box(310,245,310,160,"#EAF0F7"),box(665,245,300,160,"#EAF0F7"),box(1010,245,370,160,"#EAF0F7"),
 line("M265 150L300 150",arrow=TRUE),line("M265 325L300 325",arrow=TRUE),
 line("M620 150L655 150",arrow=TRUE),line("M620 325L655 325",arrow=TRUE),line("M965 325L1000 325",arrow=TRUE),
 line("M1195 220L1195 235",arrow=TRUE),
 line("M465 80L465 66L1195 66L1195 77",arrow=TRUE),
 label(830,60,"Contrast-based selection",24,700,teal,"middle"),
 label(142,116,"October exports",28,700,ink,"middle"),
 '<ellipse cx="142" cy="146" rx="34" ry="10" fill="#DDEEF5" stroke="#247BA0" stroke-width="2"/>',
 '<path d="M108 146L108 188C108 202 176 202 176 188L176 146" fill="#DDEEF5" stroke="#247BA0" stroke-width="2"/>',
 '<ellipse cx="142" cy="146" rx="34" ry="10" fill="#DDEEF5" stroke="#247BA0" stroke-width="2"/>',
 label(142,242,"JE2 to UniProt",26,700,blue,"middle"),
 label(142,282,"Retain sequences",24,400,muted,"middle"),label(142,314,"+ site positions",24,400,muted,"middle"),
 label(142,372,"Identifier mapping",24,700,blue,"middle"),
 label(330,115,"Total proteins",28,700),label(330,154,"TSB vs SN interactions",24),label(330,188,"Three genotype contrasts",24),
 label(685,115,"Functional patterns",28,700),label(685,154,"GO + KEGG",24,700),label(685,188,"ORA + ranked GSEA",24),
 label(1030,115,"Association network",28,700),label(1030,154,"Contrast-selected proteins",24),
 label(1030,184,"STRING 12.5 / score >=0.4",24),label(1030,211,"Stk1/Stp1 / three halos",24),
 label(330,280,"Phosphosites",28,700),label(330,318,"Mutant vs ancestor",24),label(330,349,"within TSB and SN",24),
 label(330,386,"Separate interaction tests",24,400,muted),
 label(685,280,"Overlap + sequence",28,700),label(685,320,"Venn sets / residues",24),label(685,352,"JE2 sequence logos",24),
 circle(726,382,15,"#DDEEF5"),circle(747,382,15,"#BAE2DA",teal),
 label(1030,280,"Integrate evidence",28,700),label(1030,318,"Parent-protein abundance",24),label(1030,349,"+ phosphosite signal",24),
 label(1030,386,"Prioritise candidates",24,700,teal),
 label(20,447,"I = (A_TSB - B_TSB) - (A_SN - B_SN), using log2 signals",24,700),
 label(1380,447,"Edges: functional associations",24,400,muted,"end"),'</svg>')
writeLines(enc2utf8(paste(z,collapse="\n")), out, useBytes=TRUE)
doc <- xml2::read_xml(out)
fonts <- as.numeric(xml2::xml_attr(xml2::xml_find_all(doc,"//*[local-name()='text']"),"font-size"))
stopifnot(min(fonts) >= 24)
message("Created the poster-sized Project analysis SVG; minimum printed text at 529.33mm = ",round(min(fonts)*529.333/1400*72/25.4,1),"pt.")
