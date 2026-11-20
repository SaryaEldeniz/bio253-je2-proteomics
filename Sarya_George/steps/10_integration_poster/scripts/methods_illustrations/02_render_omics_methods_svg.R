# Original editable vector artwork; icons are drawn here rather than copied
# from BioRender. Source design: course README and reports/step01_design.md.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/10_integration_poster/figures/methods_svg")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ink <- "#17394D"; secondary <- "#52717C"; blue <- "#247BA0"
teal <- "#078D92"; mint <- "#86D1C4"; pale <- "#EDF7F8"
esc <- function(x) gsub(">","&gt;",gsub("<","&lt;",gsub("&","&amp;",x,fixed=TRUE),fixed=TRUE),fixed=TRUE)
elements <- character()
add <- function(x) elements <<- c(elements,x)
text <- function(x,y,label,size=18,colour=ink,weight=400,anchor="middle") {
  add(sprintf('<text x="%s" y="%s" font-size="%s" fill="%s" font-weight="%s" text-anchor="%s">%s</text>',x,y,size,colour,weight,anchor,esc(label)))
}
path <- function(d,stroke=ink,width=2,fill="none",extra="") add(sprintf('<path d="%s" fill="%s" stroke="%s" stroke-width="%s" stroke-linecap="round" stroke-linejoin="round" %s/>',d,fill,stroke,width,extra))
line <- function(x1,y1,x2,y2,colour=ink,width=2,extra="") path(sprintf("M %s %s L %s %s",x1,y1,x2,y2),colour,width,extra=extra)
circle <- function(x,y,r,fill,stroke=ink,width=2,extra="") add(sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s" stroke="%s" stroke-width="%s" %s/>',x,y,r,fill,stroke,width,extra))
rect <- function(x,y,w,h,fill,stroke=ink,rx=10,width=2,extra="") add(sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s" stroke="%s" stroke-width="%s" %s/>',x,y,w,h,rx,fill,stroke,width,extra))
arrow <- function(x1,y1,x2,y2,colour=secondary,width=2.5) line(x1,y1,x2,y2,colour,width,'marker-end="url(#omics-arrow)"')
gstart <- function(x,y,scale=1) add(sprintf('<g transform="translate(%s,%s) scale(%s)">',x,y,scale))
gend <- function() add('</g>')

defs <- paste0('<defs>',
  '<marker id="omics-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M1 1 L9 5 L1 9" fill="none" stroke="',secondary,'" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></marker>',
  '<linearGradient id="omics-tube-glass" x1="0" x2="1"><stop offset="0" stop-color="#DDEEF1"/><stop offset="0.4" stop-color="#FFFFFF"/><stop offset="1" stop-color="#CEE6EA"/></linearGradient>',
  '<linearGradient id="omics-liquid" x1="0" x2="1"><stop offset="0" stop-color="#37A7AC"/><stop offset="0.45" stop-color="#86D1C4"/><stop offset="1" stop-color="#2A959D"/></linearGradient>',
  '</defs>')

draw_cells <- function(x,y) {
  gstart(x,y)
  coords <- rbind(c(-25,-25),c(9,-40),c(35,-10),c(-7,0),c(-37,18),c(28,30),c(-5,40))
  for(i in seq_len(nrow(coords))) {
    circle(coords[i,1],coords[i,2],23,if(i%%2) "#9DDAD1" else "#B9DFE8",teal,2.2)
    circle(coords[i,1]-6,coords[i,2]-6,13,"none","#EAF8F6",2)
    path(sprintf("M %s %s Q %s %s %s %s",coords[i,1]-17,coords[i,2]+4,coords[i,1],coords[i,2]+10,coords[i,1]+17,coords[i,2]+4),teal,1.3)
  }
  gend()
}
draw_tube <- function(x,y) {
  gstart(x,y)
  rect(-32,-61,64,13,"#D6E8EE",blue,3,2.2)
  path("M -26 -48 L -26 24 Q -26 49 0 63 Q 26 49 26 24 L 26 -48 Z",blue,2.3,"url(#omics-tube-glass)")
  path("M -24 0 L -24 23 Q -24 47 0 59 Q 24 47 24 23 L 24 0 Q 0 9 -24 0 Z",teal,1,"url(#omics-liquid)")
  path("M -14 -41 L -14 -4", "#FFFFFF",4)
  for(i in 0:3) line(13,-30+i*13,21,-30+i*13,blue,1.2)
  path("M -8 26 C -22 9 8 4 12 20 C 17 35 -12 31 -2 45",ink,2)
  gend()
}
draw_peptide <- function(x,y,phospho=FALSE,scale=1) {
  gstart(x,y,scale)
  xs <- c(-45,-23,0,23,45); ys <- c(4,-7,3,-7,4)
  for(i in 1:4) line(xs[i],ys[i],xs[i+1],ys[i+1],teal,4)
  for(i in 1:5) circle(xs[i],ys[i],9,c("#B6DDE8",mint,"#D6EBE5",blue,"#B5D8E4")[i],ink,1.2)
  if(phospho) {
    line(0,-5,5,-23,teal,1.7)
    circle(7,-30,12,"#FFFFFF",teal,1.8)
    text(7,-25,"P",15,teal,700)
  }
  gend()
}
draw_imac <- function(x,y,scale=1) {
  gstart(x,y,scale)
  path("M -25 -57 L 25 -57 L 25 35 L 8 58 L 8 73 L -8 73 L -8 58 L -25 35 Z",blue,2.2,"url(#omics-tube-glass)")
  rect(-32,-65,64,11,"#D7EBEE",blue,3,1.8)
  for(row in 0:3) for(col in 0:2) circle(-14+col*14,-35+row*16,6,if((row+col)%%2) mint else "#49A7AC",teal,1)
  line(-16,31,16,31,teal,2)
  path("M -4 82 Q -14 99 0 102 Q 14 99 4 82",teal,1.4,mint)
  circle(35,-14,12,"#FFFFFF",teal,1.8);text(35,-9,"P",15,teal,700)
  line(23,-14,17,-14,teal,1.6)
  gend()
}
draw_ms <- function(x,y) {
  gstart(x,y)
  # Generic liquid-chromatography module and mass analyser, not an instrument model.
  rect(-86,-49,57,100,"#E6F0F5",blue,8,2.2)
  rect(-78,-39,41,29,"#FFFFFF",blue,4,1.2)
  line(-69,-28,-47,-28,teal,2)
  for(i in 0:2) rect(-77,1+i*15,39,9,"#F7FBFC","#A9C7D1",2,1)
  circle(-58,43,3,mint,teal,1)
  rect(-21,-49,126,100,"#EEF5F7",blue,8,2.2)
  rect(-10,-38,104,55,"#FFFFFF","#A3C4CF",4,1.2)
  # Chromatogram above a set of mass-to-charge sticks.
  path("M -1 -12 L 12 -12 Q 18 -12 21 -23 Q 24 -35 29 -12 L 42 -12 Q 46 -12 49 -30 Q 54 -43 57 -12 L 82 -12",blue,1.9)
  line(0,9,83,9,"#BDD3DB",1)
  for(i in seq_along(c(6,10,15,21,27,34,42,48,56,64,73,79))) {
    xx<-c(6,10,15,21,27,34,42,48,56,64,73,79)[i]
    hh<-c(5,11,6,17,9,23,6,18,10,6,15,5)[i]
    line(xx,8,xx,8-hh,teal,1.7)
  }
  circle(81,34,4,mint,teal,1)
  line(-30,-25,-21,-25,blue,2)
  path("M -29 -25 C -23 -25 -26 -34 -21 -34",blue,1.7)
  gend()
}
draw_abundance <- function(x,y) {
  gstart(x,y)
  rect(-76,-45,152,100,"#FFFFFF","#C6DCE2",9,1.7)
  line(-54,32,-54,-22,secondary,1.5);line(-54,32,57,32,secondary,1.5)
  heights<-c(24,43,32,52,40)
  for(i in 1:5) rect(-41+(i-1)*20,32-heights[i],12,heights[i],c(blue,teal,mint,blue,teal)[i],"none",3,0)
  gend()
}
draw_sites <- function(x,y) {
  gstart(x,y)
  rect(-76,-53,152,107,"#FFFFFF","#C6DCE2",9,1.7)
  draw_peptide(0,-7,TRUE,.88)
  text(0,35,"Search: S / T / Y / H",14,ink,500)
  gend()
}

text(30,43,"Two complementary proteomics layers",26,ink,700,"start")
draw_cells(91,257)
text(91,351,"Harvested JE2",18,ink,600)
text(91,376,"cells",18,secondary)
arrow(159,257,218,257)
draw_tube(273,252)
text(273,351,"Protein extraction",18,ink,600)
text(273,376,"Peptide preparation",18,secondary)
arrow(315,257,403,257)
circle(422,257,5,teal,"none",0)
path("M 427 257 L 454 257 Q 470 257 470 241 L 470 177 Q 470 164 487 164 L 520 164",secondary,2.4,extra='marker-end="url(#omics-arrow)"')
path("M 427 257 L 454 257 Q 470 257 470 273 L 470 351 Q 470 364 487 364 L 520 364",secondary,2.4,extra='marker-end="url(#omics-arrow)"')

text(600,95,"Total proteome",23,blue,700)
text(600,121,"Unenriched peptides",17,secondary)
draw_peptide(600,165,FALSE,1.12)
draw_peptide(600,204,FALSE,.80)
arrow(673,177,758,177)
draw_ms(860,177)
text(870,250,"LC-MS/MS",18,ink,600)
arrow(982,177,1080,177)
draw_abundance(1190,172)
text(1190,105,"Protein abundance",21,blue,700)
text(1190,250,"Total-protein measurements",17,secondary)

text(600,284,"Phosphoproteome",23,teal,700)
text(600,310,"Fe-IMAC enrichment",17,secondary)
draw_imac(600,381,.8)
arrow(664,378,758,378)
draw_ms(860,378)
text(870,451,"LC-MS/MS",18,ink,600)
arrow(982,378,1080,378)
draw_sites(1190,378)
text(1190,300,"Phosphosite signals",21,teal,700)
text(1190,451,"Phosphorylated peptides",17,secondary)

line(30,477,1368,477,"#CCDEE4",1)
text(700,508,"27 biological samples; each sample measured in both layers",18,secondary,500)
svg <- c('<?xml version="1.0" encoding="UTF-8"?>',
  '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="520" viewBox="0 0 1400 520" role="img" aria-labelledby="omics-title omics-desc">',
  '<title id="omics-title">JE2 total-proteome and phosphoproteome acquisition</title>',
  '<desc id="omics-desc">The same27 biological samples feed total-protein and Fe-IMAC-enriched phosphosite measurements. Both branches use LC-MS/MS. Original editable vector artwork.</desc>',
  defs,'<g font-family="Arial,Helvetica,sans-serif">',elements,'</g></svg>')
writeLines(svg,file.path(out,"02_proteomics_layers.svg"),useBytes=TRUE)
stopifnot(!any(grepl('<image|<script|href=',svg)),length(grep('id="omics-',svg))>=3)
message("Saved editable original SVG:02_proteomics_layers.svg(1400x520).")
