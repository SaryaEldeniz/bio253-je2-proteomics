# Poster network figures: exact supplied protein sets, fixed central anchors,
# pathway outlines, redundant membership labels and the original three effects.
# Regions are supplied pathway annotations, not Leiden communities. Node position
# is editorial; it does not represent biological distance or regulatory direction.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
input <- file.path(REPO_ROOT,"reports/step06/poster_networks")
output <- file.path(REPO_ROOT,"figures/step06/poster_networks")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
annotation_mode <- match.arg(Sys.getenv("JE2_POSTER_ANNOTATIONS",unset="live"),c("live","course"))
input_suffix <- if(annotation_mode=="live") "_live" else ""
# A selected set can be rerendered without replacing the earlier figures. By
# default include each of the current comparison lists whose inputs are ready.
requested_sets <- Sys.getenv("JE2_POSTER_SETS",unset="")
sets <- if(nzchar(trimws(requested_sets))) {
  unique(trimws(strsplit(requested_sets,",",fixed=TRUE)[[1]]))
} else {
  c("49","67","75") |> keep(\(set)
    file.exists(file.path(input,paste0("nodes_",set,input_suffix,".csv"))) &&
    file.exists(file.path(input,paste0("edges_",set,".csv"))))
}
stopifnot(length(sets)>0L,all(grepl("^[1-9][0-9]*$",sets)))
for(set in sets) {
  needed <- file.path(input,c(paste0("nodes_",set,input_suffix,".csv"),paste0("edges_",set,".csv")))
  if(!all(file.exists(needed))) stop("Missing poster inputs: ",paste(needed[!file.exists(needed)],collapse=", "))
}
memberships <- read_csv(file.path(input,paste0("pathway_membership",input_suffix,".csv")),show_col_types=FALSE)
settings <- read_csv(file.path(input,"settings.csv"),show_col_types=FALSE)
pathways <- tibble(module=c("Purine metabolism","Oxidative phosphorylation","Nitrogen metabolism","Aminoacyl-tRNA biosynthesis"),
  code=c("P","R","N","T"),short=c("Purine metabolism","Respiration","Nitrogen metabolism","tRNA charging"),
  colour=c("#0072B2","#E69F00","#333333","#AA3377"),line=c("solid","longdash","dotdash","dotted"))
codes <- c("stk_anc","stp_anc","stp_stk")
effect_palette <- grDevices::colorRampPalette(c("#0072B2","#FCFCFC","#D55E00"))(401)
effect_colour <- function(x) {
  ans <- rep("#B8BDC4",length(x)); good <- is.finite(x)
  ans[good] <- effect_palette[1+round((pmax(-2,pmin(2,x[good]))+2)/4*400)]
  ans
}
circle_points <- function(x,y,r,n=73L) {
  a <- seq(0,2*pi,length.out=n)
  tibble(x=x+r*cos(a),y=y+r*sin(a))
}
ring_points <- function(x,y,inner,outer) {
  theta <- seq(0,2*pi,length.out=73L)
  tibble(x=c(x+outer*cos(theta),x+inner*cos(rev(theta))),
    y=c(y+outer*sin(theta),y+inner*sin(rev(theta))))
}

# Pathway-constrained packing: preserves disconnected proteins next to their
# annotation group instead of putting isolated proteins in a remote extra row.
# All positions are explicit and exported so poster edits remain reproducible.
pack_group <- function(tab,edges,centre,ncol) {
  n <- nrow(tab)
  if(!n) return(tibble(string_id=character(),x=double(),y=double()))
  degree <- table(factor(c(edges$node_a,edges$node_b),levels=tab$string_id))
  ord <- order(-as.integer(degree),tab$display_name,tab$uniprot_acc)
  ncol <- min(ncol,n); nrow_grid <- ceiling(n/ncol)
  grid <- expand_grid(row=seq_len(nrow_grid),column=seq_len(ncol)) |>
    mutate(x=(column-(ncol+1)/2)*1.90,y=((nrow_grid+1)/2-row)*1.72,
      radius=sqrt(x^2+y^2)) |> arrange(radius,row,column) |> slice_head(n=n)
  # Centre the actual occupied grid, including partially filled final rows.
  grid$x <- grid$x-mean(range(grid$x)); grid$y <- grid$y-mean(range(grid$y))
  tibble(string_id=tab$string_id[ord],x=grid$x+centre[1],y=grid$y+centre[2])
}
prepare_positions <- function(nodes,edges,membership) {
  membership <- membership |> filter(string_id %in% nodes$string_id)
  clean <- membership |> filter(!coalesce(annotation_review,FALSE)) |> select(string_id,module) |> distinct() |>
    inner_join(pathways,by="module",relationship="many-to-one")
  code_map <- clean |> group_by(string_id) |> summarise(
    membership_codes=paste(unique(code),collapse="/"),n_pathways=n_distinct(module),
    region=first(module),.groups="drop")
  nodes <- nodes |> left_join(code_map,by="string_id",relationship="one-to-one") |>
    mutate(is_anchor=node_role=="Biological anchor" | uniprot_acc %in% c("Q2FZ64","Q2FZ65"),
      region=case_when(is_anchor~"Anchors",coalesce(annotation_review,FALSE) & is.na(region)~"Assignment under review",
        is.na(region)~if_else(annotation_mode=="live","Additional selected proteins","Other / unassigned"),TRUE~region),
      membership_codes=case_when(is_anchor~"",region=="Assignment under review"~"?",
        region=="Additional selected proteins"~"C",
        is.na(membership_codes)~"-",TRUE~membership_codes))
  centres <- list("Purine metabolism"=c(-8,-.5),"Oxidative phosphorylation"=c(7.6,-3),
    "Nitrogen metabolism"=c(.6,-3.9),"Aminoacyl-tRNA biosynthesis"=c(7.3,3.4),
    "Assignment under review"=c(-8,4.7),"Other / unassigned"=c(.8,7.4))
  columns <- c("Purine metabolism"=6L,"Oxidative phosphorylation"=4L,
    "Nitrogen metabolism"=3L,"Aminoacyl-tRNA biosynthesis"=5L,
    "Assignment under review"=4L,"Other / unassigned"=4L)
  if(annotation_mode=="live") {
    n_purine <- sum(nodes$region=="Purine metabolism")
    n_context <- sum(nodes$region=="Additional selected proteins")
    centres[["Purine metabolism"]] <- c(-8,2)
    columns[["Purine metabolism"]] <- if(n_purine<=12L) 5L else 6L
    centres[["Additional selected proteins"]] <- c(-8,if(n_purine<=12L) -3.2 else -4.0)
    columns[["Additional selected proteins"]] <- if(n_context>15L) 6L else 5L
    # Extra selected proteins need space for their own enclosing outline and
    # heading. Keep earlier layouts unchanged when the original gap is enough.
    purine_half_height <- max(0L,ceiling(n_purine/columns[["Purine metabolism"]])-1L)*.86
    context_half_height <- max(0L,ceiling(n_context/columns[["Additional selected proteins"]])-1L)*.86
    context_top_limit <- 2-purine_half_height-2.56
    centres[["Additional selected proteins"]][2] <- min(
      centres[["Additional selected proteins"]][2],context_top_limit-context_half_height)
  }
  positions <- list()
  for(group in names(centres)) {
    tab <- filter(nodes,region==group)
    positions[[group]] <- pack_group(tab,edges,centres[[group]],columns[[group]])
  }
  # Central positions are fixed for biological prominence, not inferred hubness.
  anchor_pos <- nodes |> filter(is_anchor) |> transmute(string_id,
    x=if_else(uniprot_acc=="Q2FZ64",-.98,.98),y=.4)
  positions <- bind_rows(c(positions,list(anchors=anchor_pos)))
  stopifnot(nrow(positions)==nrow(nodes),!anyDuplicated(positions$string_id))
  nodes <- nodes |> left_join(positions,by="string_id",relationship="one-to-one")
  # A true multi-pathway node lies between its pathway centres and keeps all its
  # membership codes. Current inputs do not need this adjustment, but no duplicate
  # node is invented if memberships overlap in a future dataset.
  for(id in nodes$string_id[coalesce(nodes$n_pathways>1L,FALSE)]) {
    regions <- clean$module[clean$string_id==id] |> unique()
    xy <- colMeans(do.call(rbind,centres[regions]))
    nodes$x[nodes$string_id==id] <- xy[1]; nodes$y[nodes$string_id==id] <- xy[2]
  }
  # Keep labels inside centre discs. Full names remain in exported node tables.
  nodes <- nodes |> mutate(label=case_when(is_anchor & uniprot_acc=="Q2FZ64"~"Stk1",
    is_anchor & uniprot_acc=="Q2FZ65"~"Stp1",
    !is.na(display_name) & !grepl("^(B7H15_|ABD|SAOUHSC_)",display_name) & nchar(display_name)<=7 ~ display_name,
    !is.na(uniprot_acc)~uniprot_acc,TRUE~substr(string_id,nchar(string_id)-5,nchar(string_id))),
    label=label)
  nodes
}
make_outlines <- function(nodes,membership) {
  polys <- labels <- list()
  for(i in seq_len(nrow(pathways))) {
    row <- pathways[i,]
    ids <- membership |> filter(module==row$module,!coalesce(annotation_review,FALSE)) |> pull(string_id) |> unique()
    tab <- filter(nodes,string_id %in% ids,!is_anchor)
    if(!nrow(tab)) next
    cloud <- bind_rows(lapply(seq_len(nrow(tab)),function(k)circle_points(tab$x[k],tab$y[k],.97,25)))
    boundary <- cloud[chull(cloud$x,cloud$y),] |> mutate(region=row$module,colour=row$colour,line=row$line)
    polys[[row$module]] <- boundary
    # Use a labelled heading just above the region; its colour is redundant with
    # P/R/N/T badges and the boundary line type.
    labels[[row$module]] <- tibble(x=mean(range(tab$x)),y=max(tab$y)+1.24,
      label=paste0(row$code," | ",row$short),colour=row$colour)
  }
  for(group in c("Assignment under review","Other / unassigned","Additional selected proteins")) {
    tab <- filter(nodes,region==group)
    if(!nrow(tab)) next
    cloud <- bind_rows(lapply(seq_len(nrow(tab)),function(k)circle_points(tab$x[k],tab$y[k],.97,25)))
    polys[[group]] <- cloud[chull(cloud$x,cloud$y),] |> mutate(region=group,colour="#6B6B6B",line="dashed")
    labels[[group]] <- tibble(x=mean(range(tab$x)),y=max(tab$y)+1.24,
      label=paste0(if(group=="Additional selected proteins") "C | " else "? | ",group),colour="#505050")
  }
  list(polys=bind_rows(polys),labels=bind_rows(labels))
}
poster_plot <- function(nodes,edges,outlines,set_label) {
  centre <- rings <- dots <- list()
  for(i in seq_len(nrow(nodes))) {
    row <- nodes[i,]
    centre[[i]] <- circle_points(row$x,row$y,.34) |> mutate(id=i,fill=if_else(row$is_anchor,"#172B4D","white"))
    for(j in seq_along(codes)) {
      effect <- row[[paste0("halo_effect_",codes[j])]]
      q <- row[[paste0("halo_FDR_",codes[j])]]
      rings[[length(rings)+1L]] <- ring_points(row$x,row$y,.38+(j-1)*.12,.47+(j-1)*.12) |>
        mutate(id=paste(i,j),fill=effect_colour(effect))
      if(is.finite(effect) && is.finite(q) && q<.1) dots[[length(dots)+1L]] <-
        tibble(x=row$x,y=row$y+.425+(j-1)*.12)
    }
  }
  e <- edges |> left_join(select(nodes,string_id,x,y),by=c("node_a"="string_id"),relationship="many-to-one") |>
    left_join(select(nodes,string_id,xend=x,yend=y),by=c("node_b"="string_id"),relationship="many-to-one")
  boundary <- outlines$polys; headings <- outlines$labels
  limits_x <- range(c(boundary$x,nodes$x))+c(-.25,.25)
  limits_y <- range(c(boundary$y,headings$y))+c(-.20,.35)
  p <- ggplot() +
    geom_polygon(data=boundary,aes(x,y,group=region,fill=colour),alpha=.027,colour=NA) +
    geom_segment(data=e,aes(x,y,xend=xend,yend=yend,linewidth=combined_score,alpha=combined_score),colour="#697686") +
    scale_linewidth_continuous(range=c(.08,.36),limits=c(.4,1),guide="none") +
    scale_alpha_continuous(range=c(.10,.38),limits=c(.4,1),guide="none") +
    geom_polygon(data=boundary,aes(x,y,group=region,colour=colour,linetype=line),fill=NA,linewidth=.65) +
    scale_colour_identity() + scale_linetype_identity() +
    geom_polygon(data=bind_rows(rings),aes(x,y,group=id,fill=fill),colour=NA) +
    geom_polygon(data=bind_rows(centre),aes(x,y,group=id,fill=fill),colour="#26364A",linewidth=.25) +
    scale_fill_identity() +
    geom_text(data=nodes,aes(x,y=y+.055,label=label),size=2.45,fontface="bold",
      colour=ifelse(nodes$is_anchor,"white","#182A3B")) +
    geom_text(data=filter(nodes,!is_anchor),aes(x,y=y-.21,label=membership_codes),size=1.9,fontface="bold",colour="#505050") +
    geom_text(data=headings,aes(x,y,label=label,colour=colour),size=4.2,fontface="bold") +
    coord_equal(xlim=limits_x,ylim=limits_y,expand=FALSE,clip="off") +
    labs(title=paste0("PASN response network | ",set_label," proteins"),
      subtitle=paste0(nrow(edges)," functional associations | confidence >=0.4 | ",
        if(annotation_mode=="live") "STRING KEGG annotation groups" else "Course pathway annotation groups")) +
    theme_void(base_size=11) + theme(plot.title=element_text(face="bold",size=18),
      plot.subtitle=element_text(size=10,colour="#4B5A6C"),plot.margin=margin(4,5,0,5))
  if(length(dots)) p <- p+geom_point(data=bind_rows(dots),aes(x,y),size=.55,colour="#142536")
  list(plot=p,xlim=limits_x,ylim=limits_y,boundary=boundary)
}
boundary_safety <- function(nodes,boundaries) {
  anchors <- filter(nodes,is_anchor)
  distances <- numeric(); inside <- logical()
  for(region in unique(boundaries$region)) {
    poly <- filter(boundaries,.data$region==.env$region)
    following <- c(2:nrow(poly),1L)
    for(k in seq_len(nrow(anchors))) {
      px <- anchors$x[k]; py <- anchors$y[k]
      ax <- poly$x; ay <- poly$y; bx <- ax[following]; by <- ay[following]
      t <- pmax(0,pmin(1,((px-ax)*(bx-ax)+(py-ay)*(by-ay))/pmax((bx-ax)^2+(by-ay)^2,1e-12)))
      distances <- c(distances,min(sqrt((px-ax-t*(bx-ax))^2+(py-ay-t*(by-ay))^2)))
      crossed <- (ay>py)!=(by>py)
      crossing_x <- ax+(bx-ax)*(py-ay)/ifelse(abs(by-ay)<1e-12,1e-12,by-ay)
      inside <- c(inside,sum(crossed & px<crossing_x)%%2L==1L)
    }
  }
  tibble(anchors_inside_any_outline=sum(inside),minimum_anchor_boundary_margin=min(distances)-.71)
}
draw_compact_legend <- function() {
  # One compact horizontal footer replaces the former three sprawling blocks.
  grid::grid.rect(gp=grid::gpar(fill="#F8F9FA",col=NA))
  grid::grid.text("RINGS",x=.014,y=.80,just="left",gp=grid::gpar(fontsize=8,fontface="bold"))
  for(j in 3:1) grid::grid.circle(x=.035,y=.39,r=grid::unit(1.6+j*1.2,"mm"),gp=grid::gpar(fill="#C5CAD0",col="white",lwd=2))
  grid::grid.circle(x=.035,y=.39,r=grid::unit(1.7,"mm"),gp=grid::gpar(fill="#172B4D",col="white"))
  grid::grid.text("Inner stk/anc | middle stp/anc | outer stp/stk",x=.061,y=.52,just="left",gp=grid::gpar(fontsize=8.5))
  grid::grid.text("Dot: protein FDR <0.1 | grey: no suitable estimate",x=.061,y=.24,just="left",gp=grid::gpar(fontsize=8))
  grid::grid.text("INTERACTION EFFECT (log2)",x=.377,y=.80,just="left",gp=grid::gpar(fontsize=8,fontface="bold"))
  for(j in 1:101) grid::grid.rect(x=.382+(j-1)/100*.19,y=.49,width=.19/100,height=.14,
    gp=grid::gpar(fill=effect_colour(-2+(j-1)/100*4),col=NA))
  for(j in 0:4) grid::grid.text(j-2,x=.382+j*.0475,y=.24,gp=grid::gpar(fontsize=7.5))
  grid::grid.text(if(annotation_mode=="live") "STRING KEGG GROUPS" else "PATHWAY OUTLINES",x=.602,y=.80,just="left",gp=grid::gpar(fontsize=8,fontface="bold"))
  for(j in seq_len(nrow(pathways))) {
    x <- .607+((j-1)%%2)*.19; y <- .52-floor((j-1)/2)*.28
    grid::grid.lines(x=c(x,x+.025),y=c(y,y),gp=grid::gpar(col=pathways$colour[j],lty=pathways$line[j],lwd=1.8))
    grid::grid.text(paste0(pathways$code[j]," ",pathways$short[j]),x=x+.03,y=y,just="left",gp=grid::gpar(fontsize=8))
  }
}
save_poster <- function(result,nodes,edges,stem) {
  width <- 16
  body_ratio <- diff(result$xlim)/diff(result$ylim)
  # Tight aspect ratio follows the occupied drawing; footer is only 0.85 inch.
  body_height <- width/body_ratio+.68
  height <- body_height+1.10
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1,heights=grid::unit(c(body_height,.85,.25),"in"))))
    print(result$plot,newpage=FALSE,vp=grid::viewport(layout.pos.row=1))
    grid::pushViewport(grid::viewport(layout.pos.row=2)); draw_compact_legend(); grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row=3))
    grid::grid.text(paste0(if(annotation_mode=="live")
      "C: additional selected proteins outside these four KEGG maps. Anchors are centred editorially; position is not biological distance."
      else "Central anchors and pathway grouping are editorial, not biological distance. ? marks pathway membership under review.",
      "\nColours clip beyond +/-2; original values are retained. Edges are functional associations, not proof of binding or phosphorylation."),
      x=.01,y=.5,just="left",gp=grid::gpar(fontsize=7.3,col="#4B5A6C",lineheight=1.1))
    grid::popViewport(2)
  }
  grDevices::png(file.path(output,paste0(stem,".png")),width=width,height=height,units="in",res=300,bg="white")
  draw();grDevices::dev.off()
  grDevices::pdf(file.path(output,paste0(stem,".pdf")),width=width,height=height,onefile=TRUE)
  draw();grDevices::dev.off()
  grDevices::svg(file.path(output,paste0(stem,".svg")),width=width,height=height,bg="white")
  draw();grDevices::dev.off()
  d <- as.matrix(dist(as.matrix(select(nodes,x,y)))); diag(d) <- Inf
  safety <- boundary_safety(nodes,result$boundary)
  stopifnot(min(d)>1.42,safety$anchors_inside_any_outline==0L,safety$minimum_anchor_boundary_margin>0)
  bind_cols(tibble(set=stem,annotation_mode,n_nodes=nrow(nodes),n_edges=nrow(edges),width_inches=width,height_inches=height,
    minimum_centre_distance=min(d),outer_node_diameter=1.42,
    footer_fraction=1.10/height,pinned_anchors=all(abs(filter(nodes,is_anchor)$y-.4)<1e-10)),safety)
}

metrics <- list()
for(set in sets) {
  nodes <- read_csv(file.path(input,paste0("nodes_",set,input_suffix,".csv")),show_col_types=FALSE)
  edges <- read_csv(file.path(input,paste0("edges_",set,".csv")),show_col_types=FALSE)
  stopifnot(nrow(nodes)==as.integer(set),!anyDuplicated(nodes$string_id),
    all(c(edges$node_a,edges$node_b) %in% nodes$string_id))
  membership <- memberships |> filter(as.character(set_id)==set)
  positions <- prepare_positions(nodes,edges,membership)
  boundaries <- make_outlines(positions,membership)
  result <- poster_plot(positions,edges,boundaries,set)
  mode_suffix <- if(annotation_mode=="course") "_course_scope" else ""
  metrics[[set]] <- save_poster(result,positions,edges,paste0("STRING_",set,"_poster_halos",mode_suffix))
  write_csv(positions,file.path(input,paste0("poster_coordinates_",set,mode_suffix,".csv")))
  write_csv(boundaries$polys,file.path(input,paste0("poster_pathway_outlines_",set,mode_suffix,".csv")))
}
metric_path <- file.path(input,paste0("poster_geometry_checks",if(annotation_mode=="course") "_course_scope" else "",".csv"))
new_metrics <- bind_rows(metrics)
old_metrics <- if(file.exists(metric_path)) read_csv(metric_path,show_col_types=FALSE) else new_metrics[0,]
# Retain the checks for untouched sets when rendering only one new input list.
retained_metrics <- anti_join(old_metrics,select(new_metrics,set,annotation_mode),by=c("set","annotation_mode"))
write_csv(bind_rows(retained_metrics,new_metrics),metric_path)
writeLines(c("# Poster network figures","",
  "Exact 49-, 67- and 75-protein input sets can be shown separately; these are not the former 61-node one-round expansion. The 75-protein set restores all eight previously held proteins to the original 67-protein list. All supplied nodes and edges are retained, including isolates. Stk1/Stp1 are pinned symmetrically at the centre for biological prominence, not inferred hub status.","",
  "Pathway outlines use blue, orange, charcoal and purple, plus distinct line types, labels and P/R/N/T node badges. Node centres are neutral white, with dark anchors. The default uses freshly retrieved STRING KEGG annotation groups, not Leiden clustering. Additional selected proteins outside those four maps appear in a grey enclosure with C badges; this is not a claim that they lack other functions or pathways. Historical course review flags remain in the source metadata and are not automatically applied to fresh KEGG membership. Actual multiple memberships, if present, have multiple badges and participate in overlapping outlines without duplicate nodes.","",
  "Ring order is inner stk/anc, middle stp/anc, outer stp/stk. Signed original protein interaction effects use a shared ±2 colour scale; grey means unavailable/unsuitable quantitative estimate and white means zero. Black dots mark protein FDR<0.1. Position does not show biological distance, edges are functional associations, and neither establishes direct phosphorylation.","",
  "PNG is 300 dpi; PDF and SVG preserve vectors. Canvas aspect follows occupied coordinates and a compact single-row footer. Poster coordinates, outline geometry and minimum-distance checks are exported beside source tables. Set JE2_POSTER_SETS=75 to render only the new 75-protein figure; existing figures and geometry rows for other lists are preserved. Without this variable, the script renders each available 49-, 67- and 75-protein input. To reproduce optional original-course grouping, set JE2_POSTER_ANNOTATIONS=course before running; those files receive a _course_scope suffix and do not replace live-group figures."),file.path(output,"README.md"))
message("Poster networks exported for exact inputs: ",paste(sets,collapse=", "),".")
