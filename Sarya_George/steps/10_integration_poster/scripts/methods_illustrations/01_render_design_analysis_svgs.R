# Original editable vector artwork for the poster methods, using verified
# course design and completed project analyses. No BioRender assets are used.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/10_integration_poster/figures/methods_svg")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ink <- "#17394D"; secondary <- "#52717C"; blue <- "#247BA0"
teal <- "#078D92"; mint <- "#86D1C4"; pale <- "#EDF7F8"; pale_blue <- "#EAF0F7"
esc <- function(x) {
  x <- gsub("&","&amp;",x,fixed=TRUE)
  x <- gsub("<","&lt;",x,fixed=TRUE); gsub(">","&gt;",x,fixed=TRUE)
}
txt <- function(x,y,text,size=22,fill=ink,weight=400,anchor="start",italic=FALSE) paste0(
  '<text x="',x,'" y="',y,'" font-family="Arial, Helvetica, sans-serif" font-size="',size,
  '" font-weight="',weight,'" text-anchor="',anchor,'" fill="',fill,'"',
  if(italic) ' font-style="italic"' else '', '>',esc(text),'</text>')
rect <- function(x,y,w,h,fill=pale,stroke="none",r=14,sw=2) sprintf(
  '<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s" stroke="%s" stroke-width="%s"/>',x,y,w,h,r,fill,stroke,sw)
circle <- function(x,y,r,fill=mint,stroke=teal,sw=2) sprintf(
  '<circle cx="%s" cy="%s" r="%s" fill="%s" stroke="%s" stroke-width="%s"/>',x,y,r,fill,stroke,sw)
path <- function(d,stroke=teal,fill="none",sw=3,extra="") paste0(
  '<path d="',d,'" fill="',fill,'" stroke="',stroke,'" stroke-width="',sw,'" stroke-linecap="round" stroke-linejoin="round" ',extra,'/>')
group <- function(body,x=0,y=0,s=1) paste0('<g transform="translate(',x,' ',y,') scale(',s,')">',paste(body,collapse=""),'</g>')
arrow <- function(x1,y1,x2,y2,prefix,colour=teal,sw=3) path(sprintf("M %s %s L %s %s",x1,y1,x2,y2),colour,sw=sw,
  extra=paste0('marker-end="url(#',prefix,'-arrow)"'))
defs <- function(prefix) paste0('<defs><marker id="',prefix,'-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 Z" fill="',teal,'"/></marker></defs>')
write_svg <- function(file,w,h,body,prefix,title,description) {
  svg <- paste0('<svg xmlns="http://www.w3.org/2000/svg" width="',w,'" height="',h,'" viewBox="0 0 ',w,' ',h,'" role="img" aria-labelledby="',prefix,'-title ',prefix,'-desc">',
    '<title id="',prefix,'-title">',esc(title),'</title><desc id="',prefix,'-desc">',esc(description),'</desc>',
    defs(prefix),paste(body,collapse=""),'</svg>')
  writeLines(enc2utf8(svg),file.path(out,file),useBytes=TRUE)
  stopifnot(length(xml2::xml_find_all(xml2::read_xml(file.path(out,file)),"//*[local-name()='image']"))==0L)
}
cocci <- function(x,y,s=1,colour=mint) group(c(
  circle(0,0,12,colour),circle(23,-2,12,colour),circle(12,19,12,colour),
  circle(-10,21,12,colour),circle(28,21,12,colour),circle(4,40,12,colour),circle(25,43,12,colour),
  path("M -4 -7 Q 0 -10 5 -7",stroke="#FFFFFF",sw=2),path("M 20 14 Q 26 10 31 15",stroke="#FFFFFF",sw=2)),x,y,s)
tube <- function(x,y,s=1,fill=mint) group(c(
  path("M 0 5 L 0 45 Q 0 50 6 57 L 13 65 Q 17 68 21 65 L 29 57 Q 35 50 35 45 L 35 5 Z",ink,"#FFFFFF",2),
  path("M 3 31 L 32 31 L 32 45 Q 31 49 26 54 L 18 63 L 10 55 Q 3 49 3 45 Z",stroke="none",fill=fill),
  rect(-3,0,41,10,fill=pale_blue,stroke=ink,r=3,sw=2),
  path("M 7 15 L 7 27",stroke="#C7DFE5",sw=3)),x,y,s)
clock <- function(x,y,r=23) paste0(circle(x,y,r,"#FFFFFF",teal,2),path(sprintf("M %s %s L %s %s L %s %s",x,y-r*.58,x,y,x+r*.4,y),ink,sw=2.5))
flask <- function(x,y,s=1,rods=TRUE) group(c(
  path("M 30 0 L 67 0 L 67 15 L 61 15 L 61 54 L 91 113 Q 99 130 80 132 L 17 132 Q 0 130 8 113 L 38 54 L 38 15 L 30 15 Z",ink,"#FFFFFF",3),
  path("M 28 83 L 72 83 L 88 115 Q 94 126 79 126 L 18 126 Q 5 125 12 113 Z",stroke="none",fill="#D7EFEC"),
  if(rods) paste0('<g transform="rotate(-25 46 100)">',rect(27,94,32,12,mint,teal,6,1.5),'</g>',
    '<g transform="rotate(25 61 113)">',rect(47,108,27,10,"#7CC6D2",blue,5,1.5),'</g>') else "",
  path("M 22 113 Q 18 119 24 121",stroke="#FFFFFF",sw=3)),x,y,s)
gene <- function(x,y,mutant="anc") {
  # Insertion mark is on the disrupted gene, not on the bacterium itself.
  c(path(sprintf("M %s %s L %s %s",x,y+8,x+156,y+8),"#A2B7C1",sw=2),
    rect(x+3,y-4,68,24,"#DDEEF5",blue,4,1.5),rect(x+79,y-4,70,24,"#D4EEEA",teal,4,1.5),
    txt(x+37,y+13,"stp1",13,ink,anchor="middle"),txt(x+114,y+13,"stk1",13,ink,anchor="middle"),
    if(mutant!="anc") {
      xx <- x+if(mutant=="stp") 37 else 114
      path(sprintf("M %s %s L %s %s L %s %s Z",xx-10,y-23,xx+10,y-23,xx,y-7),ink,teal,1.5)
    } else "")
}
samples <- read_tsv(file.path(DATA_ROOT,"derived/samples.tsv"),show_col_types=FALSE)
design <- samples |> count(genotype,medium,name="n")
stopifnot(nrow(samples)==27L,n_distinct(samples$sample)==27L,nrow(design)==9L,all(design$n==3L))

# PANEL 1: a sterile-supernatant source plus a true nine-group/27-replicate matrix.
p <- "design"
d <- c(txt(30,42,"Experimental design",30,weight=700),txt(1370,40,"S. aureus JE2",21,secondary,anchor="end",italic=TRUE),
  flask(80,77,.78),txt(115,199,"P. aeruginosa",15,anchor="middle",italic=TRUE),
  txt(115,219,"PAO1 spent medium",15,anchor="middle"),arrow(166,129,194,129,p),
  rect(210,102,97,18,"#DDEEF5",ink,8,2),path("M 222 122 L 294 122 L 280 143 L 235 143 Z",teal,"#D7EFEC",2),
  path("M 257 145 L 257 155",teal,sw=3),circle(257,165,5,teal,stroke="none"),
  txt(257,204,"Sterile filtration",18,anchor="middle"),arrow(320,129,371,129,p),tube(393,92,1.15),
  txt(462,117,"Cell-free PA supernatant",23,weight=700),txt(462,145,"PASN: secreted products + altered nutrients",18,secondary),
  rect(945,79,425,119,pale,r=16),flask(966,92,.58,FALSE),
  txt(1045,111,"Long exposure (SN)",22,weight=700),txt(1045,140,"70% TSB + 30% PASN",21,teal,700),
  txt(1045,168,"Approximately 15 h",18,secondary))
column_x <- c(595,865,1135); groups <- c("anc","stk","stp")
for(j in seq_along(groups)) {
  x <- column_x[j]
  d <- c(d,rect(x-113,229,226,100,pale_blue,r=15),cocci(x-84,244,.57),
    txt(x-20,254,c("JE2 ancestor","JE2 stk mutant","JE2 stp mutant")[j],18.5,weight=700),
    txt(x+30,279,c("wild type","Kinase | stk1::Tn","Phosphatase | stp1::Tn")[j],13.5,secondary,anchor="middle"),gene(x-63,306,groups[j]))
}
row_y <- c(350,434,518)
for(i in 1:3) {
  y <- row_y[i]; m <- c("TSB","SNpulse","SN")[i]
  fills <- c("#EAF0F7","#F0F8F7","#D9EFEC")
  d <- c(d,rect(30,y,427,72,fills[i],r=14),clock(65,y+35,21),
    txt(104,y+28,c("TSB | control","SNpulse | early exposure","SN | long exposure")[i],22,weight=700),
    txt(104,y+54,c("Control broth","Growth in TSB, then 15 min PASN","Approximately 15 h in PASN mixture")[i],17,secondary))
  for(j in 1:3) {
    x <- column_x[j]
    d <- c(d,rect(x-113,y,226,72,fills[i],r=14))
    for(k in 1:3) d <- c(d,paste0('<g id="design-sample-',groups[j],'-',m,'-',k,'">',
      group(tube(0,0,.69,c("#A6CADC","#B8DFD7","#65B9AF")[i]),x-66+(k-1)*54,y+11),'</g>'))
  }
}
d <- c(d,txt(30,627,"3 genotypes x 3 conditions x 3 biological replicates = 27 samples",23,weight=700),
  txt(1370,627,"Each matrix tube = one replicate",17,secondary,anchor="end"))
write_svg("01_experimental_design.svg",1400,650,d,p,"JE2 experimental design",
  "PAO1 spent medium is sterile-filtered before JE2 exposure. Nine genotype-by-condition groups each have three biological samples. Pulse is15minutes; long exposure isapproximately15hours.")

# PANEL 3: distinct analysis branches; Venn selection is not an interaction hit list.
p <- "analysis"
database_icon <- function(x,y,s=1) group(c(
  '<ellipse cx="43" cy="12" rx="40" ry="12" fill="#DDEEF5" stroke="#247BA0" stroke-width="2"/>',
  path("M 3 12 L 3 76 C 3 92 83 92 83 76 L 83 12",blue,"#EDF7F8",2),
  '<ellipse cx="43" cy="12" rx="40" ry="12" fill="#DDEEF5" stroke="#247BA0" stroke-width="2"/>',
  path("M 3 43 C 3 59 83 59 83 43",blue,sw=1.5)),x,y,s)
contrast_icon <- function(x,y,s=1) group(c(path("M 45 0 L 45 92","#B2C6CD",sw=2),
  path("M 5 18 L 66 18",blue,sw=4),circle(31,18,6,blue,blue),
  path("M 28 45 L 86 45",teal,sw=4),circle(61,45,6,teal,teal),
  path("M 11 73 L 74 73",blue,sw=4),circle(39,73,6,blue,blue)),x,y,s)
enrich_icon <- function(x,y,s=1) group(c(
  path("M 0 83 L 113 83 M 0 0 L 0 83","#A2B7C1",sw=2),
  rect(10,49,13,34,mint,r=3),rect(31,31,13,52,teal,r=3),rect(52,40,13,43,blue,r=3),
  path("M 2 70 L 16 47 L 29 52 L 42 23 L 57 13 L 74 22 L 90 44 L 110 55",teal,sw=3)),x,y,s)
venn_logo <- function(x,y,s=1) group(c(circle(28,27,25,"#DDEEF5",blue,2),circle(60,27,25,"#BAE2DA",teal,2),
  txt(6,94,"A",27,blue,700),txt(27,94,"T",36,teal,700),txt(52,94,"S",47,ink,700),txt(82,94,"P",26,blue,700)),x,y,s)
network <- function(x,y,s=1) {
  xx<-c(24,89,56,144,12,121); yy<-c(25,22,84,72,119,132)
  z<-c()
  for(e in list(c(1,2),c(1,3),c(2,3),c(2,4),c(3,5),c(3,6),c(4,6))) z<-c(z,path(sprintf("M %s %s L %s %s",xx[e[1]],yy[e[1]],xx[e[2]],yy[e[2]]),"#9EBCC7",sw=2))
  for(i in seq_along(xx)) z<-c(z,circle(xx[i],yy[i],19,"none",if(i%%2) blue else teal,4),circle(xx[i],yy[i],13,"#FFFFFF",ink,1.5))
  z<-c(z,txt(xx[1],yy[1]+4,"Stk",9,ink,700,"middle"),txt(xx[2],yy[2]+4,"Stp",9,ink,700,"middle"))
  group(z,x,y,s)
}
integrate <- function(x,y,s=1) group(c(
  path("M 6 4 C 23 -8 22 25 39 12 S 61 -1 72 18",blue,sw=8),
  circle(40,10,9,mint,teal),txt(40,14,"P",10,ink,700,"middle"),
  path("M 10 45 L 72 45 M 40 46 L 40 65",secondary,sw=2),
  circle(9,84,17,pale_blue,blue),circle(41,84,17,pale,teal),circle(74,84,17,"#D9EFEC",teal),
  txt(41,117,"protein + site",16,secondary,anchor="middle")),x,y,s)
a <- c(txt(30,42,"Project analysis",30,weight=700),txt(1370,40,"TSB vs long-exposure SN is the main focus",19,secondary,anchor="end"),
  rect(30,95,235,480,pale_blue,r=20),database_icon(105,133,.9),
  txt(147,256,"October exports",22,weight=700,anchor="middle"),
  txt(147,291,"JE2 protein IDs",20,anchor="middle"),txt(147,323,"to UniProt",20,anchor="middle"),
  path("M 95 349 L 201 349",secondary,sw=2),
  txt(147,394,"Retain JE2",19,secondary,anchor="middle"),txt(147,421,"sequences and",19,secondary,anchor="middle"),
  txt(147,448,"site positions",19,secondary,anchor="middle"),
  txt(147,532,"Identifier mapping",18,blue,700,"middle"))
# Main branch panels are directly labelled, not identified by colour alone.
for(y in c(95,345)) a <- c(a,arrow(271,y+107,300,y+107,p),arrow(625,y+107,655,y+107,p))
a <- c(a,arrow(970,452,1001,452,p),
  path("M 466 94 L 466 77 Q 466 69 475 69 L 1182 69 Q 1191 69 1191 78 L 1191 91",teal,sw=2,
    extra='marker-end="url(#analysis-arrow)"'),
  txt(830,61,"Contrast-based selection",15,teal,700,"middle"),arrow(1191,324,1191,340,p,sw=2))
a <- c(a,rect(311,95,310,225,pale,r=18),contrast_icon(330,143,.73),
  txt(427,132,"Total proteins",23,weight=700),txt(427,168,"TSB vs SN interactions",18,secondary),
  txt(427,202,"stk vs anc",20),txt(427,233,"stp vs anc",20),txt(427,264,"stp vs stk",20),
  rect(666,95,300,225,pale,r=18),enrich_icon(686,153,.76),
  txt(816,132,"Functional patterns",22,weight=700,anchor="middle"),txt(791,222,"GO + KEGG",22,weight=700),
  txt(791,252,"ORA + ranked GSEA",17,secondary),txt(816,293,"Observed source estimates",17,secondary,anchor="middle"),
  rect(1011,95,359,225,pale,r=18),network(1027,144,.75),
  txt(1191,132,"Association network",22,weight=700,anchor="middle"),txt(1174,196,"STRING 12.5",21,weight=700),
  txt(1174,226,"Confidence >=0.4",17,secondary),txt(1174,255,"Stk1 / Stp1 anchors",17,secondary),txt(1191,295,"Contrast-selected proteins + 3 halos",17,secondary,anchor="middle"),
  rect(311,345,310,230,pale_blue,r=18),contrast_icon(330,392,.73),
  txt(427,383,"Phosphosites",23,weight=700),txt(427,420,"Mutant vs ancestor",18),
  txt(427,450,"within TSB and SN",18),txt(466,533,"Separate interaction tests",17,secondary,anchor="middle"),
  rect(666,345,300,230,pale_blue,r=18),venn_logo(686,407,.93),
  txt(816,383,"Overlap + sequence",22,weight=700,anchor="middle"),txt(800,447,"Venn sets",21,weight=700),
  txt(800,477,"Residue profiles",18,secondary),txt(800,507,"JE2 sequence logos",17,secondary),
  rect(1011,345,359,230,pale_blue,r=18),integrate(1037,415,1),
  txt(1191,383,"Integrate evidence",22,weight=700,anchor="middle"),txt(1177,438,"Parent protein",20,weight=700),
  txt(1177,469,"+ phosphosite signal",18),txt(1177,514,"Prioritise candidates",18,teal,700),
  txt(30,621,"Interaction = (A_TSB - B_TSB) - (A_SN - B_SN), using log2 signals",20,weight=700),
  txt(1370,621,"Network edges = functional associations",17,secondary,anchor="end"))
write_svg("03_analysis_workflow.svg",1400,650,a,p,"JE2 protein and phosphosite analysis workflow",
  "Separate protein interaction and within-medium phosphosite comparisons feed functional enrichment, association networks, Venn sets and sequence logos. Protein-site matching prioritises candidates; direct substrates or occupancy are not established.")
message("Created experimental-design and analysis SVG panels; sample matrix verified against27unique source samples.")
