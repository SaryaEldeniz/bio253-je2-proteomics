# Standalone poster labels; no analysis tables or existing plots are modified.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

export_volcano_labels <- function() {
  out <- file.path(REPO_ROOT,"steps/10_integration_poster/figures/volcano_labels")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  mm <- function(x) grid::unit(x,"mm")
  gp <- function(...) grid::gpar(fontfamily="Arial",...)
  # Measure on the same Cairo device used for the vector exports.
  measure_path <- tempfile(fileext=".svg")
  grDevices::svg(measure_path,width=20,height=2,family="Arial")
  text_width <- function(label,size=24,bold=FALSE) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(label,
      gp=gp(fontsize=size,fontface=if(bold)"bold" else "plain"))),"mm",valueOnly=TRUE)
  }
  labels <- c("Not selected","Negative","Positive","Observed","LOD-imputed")
  colours <- c("#B4BEC9","#0072B2","#D55E00","#526273","#526273")
  widths <- vapply(labels,text_width,numeric(1))
  padding <- 4; key_width <- 6; key_gap <- 3; item_gap <- 7
  legend_w <- 2*padding + sum(widths+key_width+key_gap) + 4*item_gap
  legend_h <- 28
  title_texts <- c(title_total_proteome_interactions="Total proteome: genotype \u00d7 medium interaction",
    title_phosphosite_interactions="Phosphosites: genotype \u00d7 medium interaction")
  interaction_formula <- expression(bold(Delta*Delta*log[2]~bold("abundance") ==
    log[2]*bgroup("(",frac(G[bold("1, TSB")],G[bold("1, PASN")]),")") -
    log[2]*bgroup("(",frac(G[bold("2, TSB")],G[bold("2, PASN")]),")")))
  title_widths <- vapply(title_texts,text_width,numeric(1),size=32,bold=TRUE)+4
  formula_w <- text_width(interaction_formula,24,bold=TRUE)+4
  formula_h <- grid::convertHeight(grid::grobHeight(grid::textGrob(interaction_formula,
    gp=gp(fontsize=24,fontface="bold"))),"mm",valueOnly=TRUE)+8
  grDevices::dev.off()
  unlink(measure_path)

  draw_legend <- function() {
    grid::grid.roundrect(x=mm(legend_w/2),y=mm(legend_h/2),
      width=mm(legend_w-1),height=mm(legend_h-1),r=mm(2),
      gp=grid::gpar(fill="white",col="#8A97A6",lwd=1))
    x <- padding
    for(i in seq_along(labels)) {
      if(i<5) grid::grid.circle(x=mm(x+key_width/2),y=mm(18.5),r=mm(2.75),
        gp=grid::gpar(fill=colours[i],col=NA)) else {
        grid::grid.segments(x0=mm(x+.5),y0=mm(16),x1=mm(x+5.5),y1=mm(21),
          gp=grid::gpar(col=colours[i],lwd=1.8))
        grid::grid.segments(x0=mm(x+.5),y0=mm(21),x1=mm(x+5.5),y1=mm(16),
          gp=grid::gpar(col=colours[i],lwd=1.8))
      }
      grid::grid.text(labels[i],x=mm(x+key_width+key_gap),y=mm(18.5),just="left",gp=gp(fontsize=24))
      x <- x+key_width+key_gap+widths[i]+item_gap
    }
    grid::grid.text("Selection: FDR < 0.10 and |log2 effect| > 1",
      x=mm(legend_w/2),y=mm(7.5),gp=gp(fontsize=24,col="#263A49"))
  }
  draw_title <- function(text,w) {
    grid::grid.text(text,x=mm(2),y=mm(7),just="left",gp=gp(fontsize=32,fontface="bold"))
  }
  draw_square_legend <- function(selection=TRUE) {
    side <- if(selection)100 else 80
    grid::grid.roundrect(x=mm(side/2),y=mm(side/2),width=mm(side-1),height=mm(side-1),
      r=mm(2),gp=grid::gpar(fill="white",col="#8A97A6",lwd=1))
    x <- (side-max(widths+key_width+key_gap))/2
    ys <- if(selection)c(90,78,66,52,40) else c(70,57,44,28,15)
    for(i in seq_along(labels)) {
      if(i<5) grid::grid.circle(x=mm(x+key_width/2),y=mm(ys[i]),r=mm(2.75),
        gp=grid::gpar(fill=colours[i],col=NA)) else {
        grid::grid.segments(x0=mm(x+.5),y0=mm(ys[i]-2.5),x1=mm(x+5.5),y1=mm(ys[i]+2.5),
          gp=grid::gpar(col=colours[i],lwd=1.8))
        grid::grid.segments(x0=mm(x+.5),y0=mm(ys[i]+2.5),x1=mm(x+5.5),y1=mm(ys[i]-2.5),
          gp=grid::gpar(col=colours[i],lwd=1.8))
      }
      grid::grid.text(labels[i],x=mm(x+key_width+key_gap),y=mm(ys[i]),just="left",gp=gp(fontsize=24))
    }
    if(selection) {
      grid::grid.lines(x=mm(c(8,side-8)),y=mm(c(31,31)),gp=grid::gpar(col="#CBD3DB",lwd=1))
      grid::grid.text("Selection",x=mm(side/2),y=mm(25),gp=gp(fontsize=24,fontface="bold",col="#263A49"))
      grid::grid.text("FDR < 0.10",x=mm(side/2),y=mm(16),gp=gp(fontsize=24,col="#263A49"))
      grid::grid.text("|log2 effect| > 1",x=mm(side/2),y=mm(7),gp=gp(fontsize=24,col="#263A49"))
    }
  }
  export <- function(name,w,h,draw) {
    svg_path <- file.path(out,paste0(name,".svg"))
    grDevices::svg(svg_path,width=w/25.4,height=h/25.4,family="Arial",bg="transparent")
    grid::grid.newpage();draw();grDevices::dev.off()
    svg <- readLines(svg_path,warn=FALSE)
    root <- grep("^<svg ",svg)[1]
    svg[root] <- sub('width="[^"]+"',sprintf('width="%.6fmm"',w),svg[root])
    svg[root] <- sub('height="[^"]+"',sprintf('height="%.6fmm"',h),svg[root])
    stopifnot(!any(grepl("<image",svg,fixed=TRUE)))
    writeLines(svg,svg_path,useBytes=TRUE)
    grDevices::png(file.path(out,paste0(name,".png")),width=w,height=h,units="mm",res=180,bg="white",type="cairo")
    grid::grid.newpage();draw();grDevices::dev.off()
    tibble(file=paste0(name,".svg"),width_mm=w,height_mm=h,font_pt=if(startsWith(name,"title_"))32 else 24)
  }
  manifest <- list(export("volcano_shared_legend_box",legend_w,legend_h,draw_legend))
  for(i in seq_along(title_texts)) {
    manifest[[i+1]] <- export(names(title_texts)[i],title_widths[i],14,
      function() draw_title(title_texts[i],title_widths[i]))
  }
  manifest[[length(manifest)+1]] <- export("interaction_formula_G",formula_w,formula_h,
    function() grid::grid.text(interaction_formula,x=mm(2),y=mm(formula_h/2),
      just="left",gp=gp(fontsize=24,fontface="bold")))
  manifest[[length(manifest)+1]] <- export("volcano_shared_legend_square",100,100,
    function() draw_square_legend(TRUE))
  manifest[[length(manifest)+1]] <- export("volcano_shared_legend_square_no_selection",80,80,
    function() draw_square_legend(FALSE))
  write_tsv(bind_rows(manifest),file.path(out,"asset_sizes.tsv"))
  zip::zipr(file.path(out,"volcano_legend_and_titles_SVGs.zip"),
    files=file.path(out,bind_rows(manifest)$file))
  print(bind_rows(manifest))
}
export_volcano_labels()
