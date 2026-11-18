# Reproducible STRING networks with three concentric interaction-effect rings.
# Uses only base R, grid and the project's existing ggplot2/tidyverse packages.
# Node positions show association topology, not biological distance or causality.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
input <- file.path(REPO_ROOT,"reports/step06/string_network")
out <- file.path(REPO_ROOT,"figures/step06/string_network")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
nodes <- read_csv(file.path(input,"nodes.csv"),show_col_types=FALSE)
seed_edges <- read_csv(file.path(input,"edges_seed_export.csv"),show_col_types=FALSE)
high_edges <- read_csv(file.path(input,"edges_seed_high_confidence.csv"),show_col_types=FALSE)
expanded_edges <- read_csv(file.path(input,"edges_expanded.csv"),show_col_types=FALSE)
settings <- read_csv(file.path(input,"network_settings.csv"),show_col_types=FALSE)
threshold <- as.numeric(settings$value[settings$key=="main_threshold"])
if(length(threshold)!=1L || !is.finite(threshold)) threshold <- .4
stopifnot(!anyDuplicated(nodes$string_id),all(c("node_a","node_b","combined_score") %in% names(expanded_edges)))
codes <- c("stk_anc","stp_anc","stp_stk")
roles <- c("Selected seed","Biological anchor","Added partner")
module_colours <- c("Purine metabolism"="#DAE5F4", "Oxidative phosphorylation"="#D9EAD3",
  "Nitrogen metabolism"="#F8E0B6", "Aminoacyl-tRNA biosynthesis"="#E9DDF2",
  "Other or unassigned"="#E4E7EA", "Regulatory anchors"="#145A70",
  "Multiple selected modules"="#F1D6D6")
colour_ramp <- grDevices::colorRampPalette(c("#0072B2","#FAFAFA","#D55E00"))(401)
effect_colours <- function(x) {
  answer <- rep("#BFC4C9",length(x)); finite <- is.finite(x)
  answer[finite] <- colour_ramp[1L+round((pmax(-2,pmin(2,x[finite]))+2)/4*400)]
  answer
}

# Connected components and a deterministic force layout are implemented directly
# to avoid installing new graph packages in the student's course environment.
component_ids <- function(ids,edges) {
  adjacency <- matrix(FALSE,length(ids),length(ids),dimnames=list(ids,ids))
  if(nrow(edges)) {
    a <- match(edges$node_a,ids); b <- match(edges$node_b,ids)
    ok <- !is.na(a) & !is.na(b) & a!=b
    adjacency[cbind(a[ok],b[ok])] <- TRUE
    adjacency[cbind(b[ok],a[ok])] <- TRUE
  }
  membership <- integer(length(ids)); component <- 0L
  for(i in seq_along(ids)) if(membership[i]==0L) {
    component <- component+1L; frontier <- i; membership[i] <- component
    while(length(frontier)) {
      neighbours <- which(colSums(adjacency[frontier,,drop=FALSE])>0 & membership==0L)
      membership[neighbours] <- component; frontier <- neighbours
    }
  }
  list(membership=membership,adjacency=adjacency)
}
force_component <- function(adjacency,seed=20261117L) {
  n <- nrow(adjacency)
  if(n==1L) return(matrix(c(0,0),ncol=2))
  if(n==2L) return(cbind(c(-1.2,1.2),c(0,0)))
  set.seed(seed)
  angle <- seq(0,2*pi,length.out=n+1L)[seq_len(n)]
  xy <- cbind(cos(angle),sin(angle))*sqrt(n)+matrix(runif(n*2,-.1,.1),ncol=2)
  k <- 1.8
  for(iteration in seq_len(800L)) {
    dx <- outer(xy[,1],xy[,1],"-"); dy <- outer(xy[,2],xy[,2],"-")
    distance <- sqrt(dx^2+dy^2); diag(distance) <- Inf
    distance <- pmax(distance,.05)
    repulsion <- k^2/distance^2
    attraction <- adjacency*distance/k
    attraction[!is.finite(attraction)] <- 0
    force <- cbind(rowSums(dx*(repulsion-attraction)),rowSums(dy*(repulsion-attraction)))-xy*.025
    magnitude <- sqrt(rowSums(force^2))
    temperature <- .45*(1-iteration/801)^1.3+.004
    step <- pmin(magnitude,temperature)/pmax(magnitude,1e-10)
    xy <- xy+force*step
    xy <- sweep(xy,2,colMeans(xy))
  }
  # Sparse tails and densely connected cliques can make a raw force layout very
  # elongated. Compress outliers monotonically, then separate colliding centres
  # rather than scaling the entire drawing by one unusually close pair.
  scale <- apply(abs(xy),2,median)
  xy <- asinh(sweep(xy,2,pmax(scale,.1),"/"))
  extent <- c(max(8,sqrt(n)*3.0),max(7,sqrt(n)*2.5))
  for(axis in 1:2) {
    xy[,axis] <- (xy[,axis]-mean(range(xy[,axis])))/max(diff(range(xy[,axis])),.1)*extent[axis]
  }
  for(iteration in seq_len(1200L)) {
    dx <- outer(xy[,1],xy[,1],"-"); dy <- outer(xy[,2],xy[,2],"-")
    distance <- pmax(sqrt(dx^2+dy^2),1e-5); diag(distance) <- Inf
    overlap <- pmax(0,2.20-distance)/distance
    xy <- xy+.22*cbind(rowSums(dx*overlap),rowSums(dy*overlap))
  }
  xy
}
layout_all <- function(tab,edges) {
  ids <- tab$string_id
  components <- component_ids(ids,edges)
  groups <- split(seq_along(ids),components$membership)
  groups <- groups[order(-lengths(groups),vapply(groups,function(i)min(ids[i]),character(1)))]
  positions <- vector("list",length(groups))
  cursor_x <- cursor_y <- row_height <- 0
  wrap_width <- max(20,sqrt(length(ids))*3.1)
  for(i in seq_along(groups)) {
    indices <- groups[[i]]
    xy <- force_component(components$adjacency[indices,indices,drop=FALSE],seed=20261117L+i)
    xy <- sweep(xy,2,apply(xy,2,min))
    width <- diff(range(xy[,1]))+2.2; height <- diff(range(xy[,2]))+2.6
    if(cursor_x>0 && cursor_x+width>wrap_width) {
      cursor_x <- 0; cursor_y <- cursor_y+row_height+1.7; row_height <- 0
    }
    positions[[i]] <- tibble(string_id=ids[indices],x=xy[,1]+cursor_x+1.1,
      y=-(xy[,2]+cursor_y+1.3),component_id=components$membership[indices])
    cursor_x <- cursor_x+width+1.2; row_height <- max(row_height,height)
  }
  bind_rows(positions)
}
positions <- layout_all(nodes,expanded_edges)
write_csv(positions,file.path(input,"network_plot_coordinates.csv"))
nodes <- nodes |> left_join(positions,by="string_id",relationship="one-to-one") |>
  mutate(module_group=if_else(module_group %in% names(module_colours),module_group,"Other or unassigned"),
    centre_colour=unname(module_colours[module_group]),
    label=case_when(!is.na(display_name) & nzchar(display_name) & !grepl("^(B7H15_|ABD|SAOUHSC_)",display_name) ~ display_name,
      !is.na(uniprot_acc) & nzchar(uniprot_acc) ~ uniprot_acc,
      TRUE ~ sub("^93061\\.","",string_id)),
    label=if_else(nchar(label)>8L,paste0(substr(label,1L,ceiling(nchar(label)/2)),"\n",
      substring(label,ceiling(nchar(label)/2)+1L)),label),
    label=if_else(coalesce(annotation_review,FALSE),paste0(label,"\u2020"),label))

# Annular polygons encode effect only by colour; width and area never vary by FC.
annulus <- function(x,y,r_inner,r_outer) {
  theta <- seq(0,2*pi,length.out=65L)
  tibble(x=c(x+r_outer*cos(theta),x+r_inner*cos(rev(theta))),
    y=c(y+r_outer*sin(theta),y+r_inner*sin(rev(theta))))
}
circle <- function(x,y,r) {
  theta <- seq(0,2*pi,length.out=65L)
  tibble(x=x+r*cos(theta),y=y+r*sin(theta))
}
draw_legend <- function() {
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,3,widths=grid::unit(c(.37,.29,.34),"null"))))
  grid::pushViewport(grid::viewport(layout.pos.col=1))
  grid::grid.text("THREE INTERACTION RINGS",x=.03,y=.94,just="left",gp=grid::gpar(fontface="bold",fontsize=9))
  grid::grid.circle(x=.12,y=.47,r=grid::unit(10,"mm"),gp=grid::gpar(fill="#DCE1E5",col="white"))
  grid::grid.circle(x=.12,y=.47,r=grid::unit(7.9,"mm"),gp=grid::gpar(fill="white",col="white"))
  grid::grid.circle(x=.12,y=.47,r=grid::unit(7.3,"mm"),gp=grid::gpar(fill="#DCE1E5",col="white"))
  grid::grid.circle(x=.12,y=.47,r=grid::unit(5.4,"mm"),gp=grid::gpar(fill="white",col="white"))
  grid::grid.circle(x=.12,y=.47,r=grid::unit(4.8,"mm"),gp=grid::gpar(fill="#DCE1E5",col="white"))
  grid::grid.circle(x=.12,y=.47,r=grid::unit(3.1,"mm"),gp=grid::gpar(fill="#52677F",col="white"))
  grid::grid.text("Inner: stk / anc\nMiddle: stp / anc\nOuter: stp / stk",x=.25,y=.52,just="left",gp=grid::gpar(fontsize=9,lineheight=1.4))
  grid::grid.text("Black dot on a ring: protein FDR < 0.1",x=.03,y=.08,just="left",gp=grid::gpar(fontsize=8))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.col=2))
  grid::grid.text("SIGNED EFFECT (log2)",x=.03,y=.94,just="left",gp=grid::gpar(fontface="bold",fontsize=9))
  for(j in seq_len(101L)) grid::grid.rect(x=.08+(j-1)/100*.78,y=.62,width=.78/100,height=.13,
    gp=grid::gpar(fill=effect_colours(-2+(j-1)/100*4),col=NA))
  for(j in 0:4) grid::grid.text(as.character(j-2),x=.08+j*.195,y=.43,gp=grid::gpar(fontsize=8))
  grid::grid.rect(x=.1,y=.23,width=.06,height=.08,gp=grid::gpar(fill="#BFC4C9",col=NA))
  grid::grid.text("Unavailable / unsuitable estimate",x=.17,y=.23,just="left",gp=grid::gpar(fontsize=7.5))
  grid::grid.text("Values outside +/-2 are colour-clipped.",x=.03,y=.08,just="left",gp=grid::gpar(fontsize=7.5))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.col=3))
  grid::grid.text("NODE CENTRES / STATUS",x=.03,y=.94,just="left",gp=grid::gpar(fontface="bold",fontsize=9))
  mods <- names(module_colours)
  short <- c("Purine","Respiration","Nitrogen","tRNA charging","Other / unassigned","Stk1 / Stp1","Multiple modules")
  for(j in seq_along(mods)) {
    column <- if(j<=4) 0 else .48; row <- if(j<=4) j else j-4
    grid::grid.circle(x=.06+column,y=.81-(row-1)*.145,r=grid::unit(1.8,"mm"),gp=grid::gpar(fill=module_colours[[mods[j]]],col="#52677F"))
    grid::grid.text(short[j],x=.10+column,y=.81-(row-1)*.145,just="left",gp=grid::gpar(fontsize=7.5))
  }
  grid::grid.text("Dashed centre outline: added partner\n\u2020: pathway membership under review",x=.03,y=.12,just="left",gp=grid::gpar(fontsize=7.5,lineheight=1.2))
  grid::popViewport(2)
}
network_plot <- function(tab,edges,title,subtitle) {
  stopifnot(all(edges$node_a %in% tab$string_id),all(edges$node_b %in% tab$string_id))
  rings <- list(); centres <- list(); dots <- list()
  for(i in seq_len(nrow(tab))) {
    for(j in 1:3) {
      value <- tab[[paste0("halo_effect_",codes[j])]][i]
      fdr <- tab[[paste0("halo_FDR_",codes[j])]][i]
      rings[[length(rings)+1L]] <- annulus(tab$x[i],tab$y[i],.40+(j-1)*.13,.51+(j-1)*.13) |>
        mutate(id=paste(i,j,sep="_"),fill=effect_colours(value))
      if(is.finite(value) && is.finite(fdr) && fdr<.1) {
        theta <- pi/2 # same radial axis makes nested ring order explicit
        dots[[length(dots)+1L]] <- tibble(x=tab$x[i]+(.455+(j-1)*.13)*cos(theta),
          y=tab$y[i]+(.455+(j-1)*.13)*sin(theta))
      }
    }
    centres[[i]] <- circle(tab$x[i],tab$y[i],.37) |>
      mutate(id=i,fill=tab$centre_colour[i],line=if_else(tab$node_role[i]=="Added partner","dashed","solid"))
  }
  edge_positions <- edges |> left_join(select(tab,string_id,x,y),by=c("node_a"="string_id"),relationship="many-to-one") |>
    left_join(select(tab,string_id,xend=x,yend=y),by=c("node_b"="string_id"),relationship="many-to-one")
  p <- ggplot() +
    geom_segment(data=edge_positions,aes(x,y,xend=xend,yend=yend,linewidth=combined_score,alpha=combined_score),colour="#738699") +
    scale_linewidth_continuous(range=c(.12,.5),limits=c(0,1),guide="none") +
    scale_alpha_continuous(range=c(.08,.50),limits=c(0,1),guide="none") +
    geom_polygon(data=bind_rows(rings),aes(x,y,group=id,fill=fill),colour=NA) +
    geom_polygon(data=bind_rows(centres),aes(x,y,group=id,fill=fill,linetype=line),colour="#33475B",linewidth=.4) +
    scale_fill_identity() + scale_linetype_identity() +
    geom_text(data=tab,aes(x=x,y=y,label=label),size=2.25,fontface="bold",
      colour=ifelse(tab$node_role=="Biological anchor","white","#203449")) +
    coord_equal(xlim=range(tab$x)+c(-1.2,1.2),ylim=range(tab$y)+c(-1.1,1.1),expand=FALSE,clip="off") +
    labs(title=title,subtitle=subtitle,
      caption="Positive = SN-TSB response more positive in anc (stk/anc, stp/anc) or stk (stp/stk); negative = converse.\nSTRING edges are functional associations; their width reflects confidence, not binding strength or phosphorylation.") +
    theme_void(base_size=11) + theme(plot.title=element_text(face="bold",size=18),
      plot.subtitle=element_text(size=10,colour="#52677F"),plot.caption=element_text(hjust=0,size=8),
      plot.margin=margin(10,12,6,12))
  if(length(dots)) p <- p+geom_point(data=bind_rows(dots),aes(x,y),size=.65,colour="#1B2633")
  p
}
save_network <- function(tab,edges,stem,title,subtitle,width=16,height=13) {
  p <- network_plot(tab,edges,title,subtitle)
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1,heights=grid::unit(c(.83,.17),"null"))))
    print(p,newpage=FALSE,vp=grid::viewport(layout.pos.row=1))
    grid::pushViewport(grid::viewport(layout.pos.row=2)); draw_legend(); grid::popViewport(2)
  }
  grDevices::png(file.path(out,paste0(stem,".png")),width=width,height=height,units="in",res=220,bg="white")
  draw(); grDevices::dev.off()
  grDevices::pdf(file.path(out,paste0(stem,".pdf")),width=width,height=height,onefile=TRUE)
  draw(); grDevices::dev.off()
  grDevices::svg(file.path(out,paste0(stem,".svg")),width=width,height=height,bg="white")
  draw(); grDevices::dev.off()
}

seeds <- filter(nodes,node_role!="Added partner")
save_network(seeds,seed_edges,"01_selected_49_export_halos","Your selected proteins: exported STRING network",
  paste0(nrow(seeds)," nodes | ",nrow(seed_edges)," exported functional associations | original selected proteins retained, including isolates"))
save_network(seeds,high_edges,"02_selected_49_high_confidence_halos","Selected proteins: high-confidence sensitivity view",
  paste0(nrow(seeds)," nodes | ",nrow(high_edges)," associations at STRING confidence >= 0.7 | same coordinates as the expanded overview"))
save_network(nodes,expanded_edges,"03_expanded_network_halos","Selected proteins and one-round STRING expansion",
  paste0(nrow(nodes)," nodes | ",nrow(expanded_edges)," associations at confidence >= ",threshold,
    " | dashed centres identify added partners"),width=18,height=15)

# Anchor detail includes every retained direct neighbour; a blank bridge is not
# filled with a hypothesised biological link. Unconnected anchors remain visible.
anchors <- filter(nodes,node_role=="Biological anchor")$string_id
anchor_edges <- expanded_edges |> filter(node_a %in% anchors | node_b %in% anchors)
detail_ids <- unique(c(anchors,anchor_edges$node_a,anchor_edges$node_b))
detail_nodes <- nodes |> filter(string_id %in% detail_ids)
detail_edges <- expanded_edges |> filter(node_a %in% detail_ids,node_b %in% detail_ids)
if(nrow(detail_nodes)>0L) {
  local_positions <- layout_all(detail_nodes,detail_edges)
  detail_nodes <- detail_nodes |> select(-x,-y,-component_id) |>
    left_join(local_positions,by="string_id",relationship="one-to-one")
  save_network(detail_nodes,detail_edges,"04_Stk1_Stp1_neighbourhood_halos","Stk1 / Stp1: retained direct neighbourhood",
    paste0(nrow(detail_nodes)," nodes | ",nrow(detail_edges)," retained associations | includes both biological anchors"),width=13,height=10)
}

# Module zooms preserve the chosen module seeds (including isolates), both
# anchors, and only retained added partners directly linked to that module.
# The other module seeds are not shown, so these are explicitly scoped views.
module_stems <- c("Purine metabolism"="05_purine_module_halos",
  "Oxidative phosphorylation"="06_respiration_module_halos",
  "Nitrogen metabolism"="07_nitrogen_module_halos",
  "Aminoacyl-tRNA biosynthesis"="08_tRNA_module_halos")
module_manifest <- list()
for(module in names(module_stems)) {
  view_stem <- module_stems[[module]]
  seed_ids <- seeds$string_id[!is.na(seeds$modules) & grepl(module,seeds$modules,fixed=TRUE)]
  adjacent <- expanded_edges |> filter(node_a %in% seed_ids | node_b %in% seed_ids)
  added_ids <- intersect(filter(nodes,node_role=="Added partner")$string_id,unique(c(adjacent$node_a,adjacent$node_b)))
  ids <- unique(c(seed_ids,anchors,added_ids))
  tab <- nodes |> filter(string_id %in% ids)
  edges <- expanded_edges |> filter(node_a %in% ids,node_b %in% ids)
  coordinates <- layout_all(tab,edges)
  tab <- tab |> select(-x,-y,-component_id) |> left_join(coordinates,by="string_id",relationship="one-to-one")
  save_network(tab,edges,view_stem,paste0(module,": selected seeds and retained context"),
    paste0(length(seed_ids)," module seeds + 2 anchors + ",length(added_ids)," added neighbours | ",
      nrow(edges)," associations at confidence >= ",threshold),width=14,height=11)
  module_manifest[[module]] <- tibble(module=module,string_id=tab$string_id,view=view_stem)
}
write_csv(bind_rows(module_manifest),file.path(input,"module_plot_node_membership.csv"))
writeLines(c("# Three-halo network figures", "",
  "All original selected proteins and both anchors remain in the seed views, including disconnected nodes. Added partners are distinguished by dashed centre outlines. A dagger preserves a supplied pathway-membership review flag.","",
  "Ring order: inner = stk/anc; middle = stp/anc; outer = stp/stk. Ring width is constant. Colour is the signed protein interaction effect, with a common -2 to +2 scale. Effects outside that range are colour-clipped only; original values are preserved in nodes.csv. Grey means unavailable/unsuitable quantitative estimate, not zero. Black ring dots indicate protein FDR<0.1.","",
  "Seed and full overview figures use the same deterministic expanded-network coordinates, exported in reports/step06/string_network/network_plot_coordinates.csv. The anchor-detail and four module zooms are separately laid out for legibility. Module zooms contain that module's original seeds plus both anchors and retained added partners adjacent to a module seed; they omit other module seeds, so the complete network is always the overview. Membership is exported in module_plot_node_membership.csv. Coordinates, proximity, and components are drawing aids rather than evidence of biological distance or direction.","",
  "Node centre colours use the supplied pathway scope; they are not assertions that every membership is canonical. STRING confidence controls edges, whereas protein FDR and effect describe the experiment. An edge is not direct phosphorylation evidence. No NES was used to select or colour proteins.","",
  "The source script uses base-R force layout/component packing and ggplot2 annular polygons, avoiding extra graph-package dependencies."),file.path(out,"README.md"))
message("Eight three-halo network views exported as PNG, PDF and SVG.")
