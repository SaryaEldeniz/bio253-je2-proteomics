# Compact matched poster volcanoes. Source estimates and tests are never refit.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")

compact_poster_volcanoes <- function(sections=getOption("poster.volcano.sections",
    c("mutants","total","phospho_interactions"))) {
  stopifnot(length(sections)>0,all(sections %in% c("mutants","total","phospho_interactions")))
  total_dir <- file.path(REPO_ROOT,"steps/04_total_proteome/figures")
  phos_dir <- file.path(REPO_ROOT,"steps/08_mutant_venn/figures/volcanoes/compact_poster")
  interaction_dir <- file.path(REPO_ROOT,"steps/07_phosphosite_results/figures/volcanoes/compact_poster")
  compact_total_dir <- file.path(total_dir,"compact_poster")
  report_dir <- file.path(REPO_ROOT,"steps/10_integration_poster/reports/compact_volcanoes")
  cache_dir <- file.path(REPO_ROOT,"results/compact_poster_volcanoes")
  for (p in c(phos_dir,interaction_dir,compact_total_dir,report_dir,cache_dir)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
  mutants <- c("stk_vs_anc_at_TSB","stp_vs_anc_at_TSB","stk_vs_anc_at_SN","stp_vs_anc_at_SN")
  interactions <- filter(step04_contrast_plan(),contrast_type=="interaction")$contrast
  short_mutant <- c("TSB: stk vs ancestor","TSB: stp vs ancestor","PASN: stk vs ancestor","PASN: stp vs ancestor")
  short_interaction <- c("stk vs ancestor","stp vs ancestor","stp vs stk")
  phos_source <- read_tsv(file.path(REPO_ROOT,"results/step07_volcano_reproduction/all_annotated_coordinates.tsv.gz"),show_col_types=FALSE)
  phos <- phos_source |> filter(contrast %in% mutants,eligible_je2)
  # Keep the exact all-model population of the requested Step 7 source figure.
  phos_interactions <- phos_source |> filter(contrast %in% interactions)
  omitted <- phos_interactions |>
    filter(contrast==interactions[1],site_id=="lcl_CP020619.1_prot_ARG46254.1_1820_S538")
  stopifnot(nrow(omitted)==1L,omitted$estimate_type=="lod_imputed",-log10(omitted$FDR)>5,
    sum(-log10(phos_interactions$FDR)>5)==1L)
  total <- read_tsv(file.path(REPO_ROOT,"results/step04/total/annotated_all.tsv.gz"),show_col_types=FALSE) |>
    filter(eligible_cds,contrast %in% interactions)
  targets <- tibble::tribble(~site_id,~arrow_label,
    "lcl_CP020619.1_prot_ARG47142.1_2776_S194","IsaA S194",
    "lcl_CP020619.1_prot_ARG46254.1_1820_T537","Pyk T537",
    "lcl_CP020619.1_prot_ARG45539.1_1070_T552","QoxB T552")
  stopifnot(!anyDuplicated(phos[c("site_id","contrast")]),
            all(phos$estimate_type %in% c("observed","lod_imputed")),
            all(phos$selected_hit == (phos$FDR<.1 & abs(phos$diff)>1)))

  # Reuse the proven collision/arrow algorithm without executing the original
  # rendering script or rewriting any analysis tables. Only label bounds vary.
  parsed <- parse(file.path(REPO_ROOT,"steps/07_phosphosite_results/scripts/07_recreate_phospho_volcanoes.R"))
  definition <- Filter(function(e) is.call(e) && identical(e[[1]],as.name("<-")) &&
    identical(e[[2]],as.name("recreate_phospho_volcanoes")),as.list(parsed))[[1]]
  source_function <- eval(definition[[3]])
  helpers <- new.env(parent=environment())
  for (statement in as.list(body(source_function))[-1]) {
    if (is.call(statement) && identical(statement[[1]],as.name("<-")) &&
        as.character(statement[[2]]) %in% c("segment_crosses_label","place_site_labels")) eval(statement,helpers)
  }
  label_body <- as.list(body(helpers$place_site_labels))
  for (i in seq_along(label_body)) {
    e <- label_body[[i]]
    if (is.call(e) && identical(e[[1]],as.name("<-"))) {
      if (identical(e[[2]],as.name("half_width"))) label_body[[i]] <- quote(half_width <- poster_half_width)
      if (identical(e[[2]],as.name("half_height"))) label_body[[i]] <- quote(half_height <- poster_half_height)
    }
  }
  body(helpers$place_site_labels) <- as.call(label_body)
  # Large poster labels use the empty upper portions of the volcanoes first.
  # This keeps their boxes off the dense near-threshold points and crosses.
  prefer_upper_space <- function(e) {
    if(!is.call(e)) return(e)
    if(identical(e[[1]],as.name("<-")) && identical(e[[2]],as.name("distance"))) {
      e[[3]] <- substitute(BASE + 3 * (1-options$label_y/y_limit)^2,list(BASE=e[[3]]))
      return(e)
    }
    as.call(lapply(as.list(e),prefer_upper_space))
  }
  body(helpers$place_site_labels) <- prefer_upper_space(body(helpers$place_site_labels))

  style <- function(p,title,subtitle,individual=FALSE) {
    p + labs(title=if(individual) NULL else title,subtitle=if(individual) NULL else subtitle,
             caption=NULL,colour=NULL,shape=NULL,y="-log10(FDR)") +
      theme_classic(base_size=24,base_family="Arial") +
      theme(plot.title=element_text(size=32,face="bold",margin=margin(b=3)),
        plot.subtitle=element_text(size=24,margin=margin(b=5)),
        strip.background=element_blank(),strip.text=element_text(size=26,face="bold",hjust=0,margin=margin(b=4)),
        axis.title=element_text(size=24),axis.text=element_text(size=24),
        axis.title.x=element_text(margin=margin(t=4)),axis.title.y=element_text(margin=margin(r=4)),
        legend.position=if(individual) "none" else "bottom",legend.box="horizontal",
        legend.text=element_text(size=24),legend.title=element_blank(),
        legend.key.width=grid::unit(8,"mm"),legend.key.height=grid::unit(8,"mm"),
        legend.spacing.x=grid::unit(9,"pt"),legend.box.spacing=grid::unit(2,"pt"),
        legend.margin=margin(0,0,0,0),panel.spacing=grid::unit(7,"mm"),
        plot.margin=margin(3,5,2,3)) +
      guides(colour=guide_legend(order=1,nrow=1,override.aes=list(shape=16,alpha=1,size=5.5)),
        shape=guide_legend(order=2,nrow=1,override.aes=list(colour="#526273",alpha=1,size=4.8,stroke=.8)))
  }
  axis_signature <- function(p) lapply(ggplot_build(p)$layout$panel_params,function(a)
    list(x=a$x.range,y=a$y.range,xb=a$x$breaks,yb=a$y$breaks,xl=a$x$get_labels(),yl=a$y$get_labels()))
  draw_metrics <- function(p,w,h) {
    path <- tempfile(tmpdir=cache_dir,fileext=".svg")
    grDevices::svg(path,width=w/25.4,height=h/25.4,family="Arial")
    on.exit(grDevices::dev.off(),add=TRUE)
    grid::grid.newpage();grid::grid.draw(ggplotGrob(p));grid::grid.force()
    vp <- unique(grep("^panel[.-]",grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name,value=TRUE))[1]
    grid::seekViewport(vp)
    pw <- grid::convertWidth(grid::unit(1,"npc"),"mm",valueOnly=TRUE)
    ph <- grid::convertHeight(grid::unit(1,"npc"),"mm",valueOnly=TRUE)
    grid::upViewport(0)
    label <- grid::textGrob("QoxB\nT552",gp=grid::gpar(fontfamily="Arial",fontsize=24,lineheight=.9))
    list(panel_w=pw,panel_h=ph,label_w=grid::convertWidth(grid::grobWidth(label),"mm",valueOnly=TRUE)+5,
         label_h=grid::convertHeight(grid::grobHeight(label),"mm",valueOnly=TRUE)+5)
  }
  audits <- list();label_reports <- list()
  save_plot <- function(p,path,w,h) {
    for(ext in c("svg","png","pdf")) {
      args <- list(filename=paste0(path,".",ext),plot=p,width=w,height=h,units="mm",bg="white",limitsize=FALSE)
      if(ext=="svg") args$device <- grDevices::svg
      if(ext=="pdf") args$device <- grDevices::cairo_pdf
      if(ext=="png") args$dpi <- 180
      do.call(ggsave,args)
    }
    svg <- readLines(paste0(path,".svg"),warn=FALSE)
    stopifnot(!any(grepl("<image",svg,fixed=TRUE)))
    # Cairo writes unitless SVG dimensions. Explicit mm retain the intended
    # physical size (and poster font size) when inserted into PowerPoint.
    root_line <- grep("^<svg ",svg)[1]
    stopifnot(length(root_line)==1L,!is.na(root_line))
    svg[root_line] <- sub('width="[^"]+"',sprintf('width="%.6fmm"',w),svg[root_line])
    svg[root_line] <- sub('height="[^"]+"',sprintf('height="%.6fmm"',h),svg[root_line])
    writeLines(svg,paste0(path,".svg"),useBytes=TRUE)
  }
  make_phos <- function(which=mutants,w=529.333,h=190,individual=FALSE,
                        interaction=FALSE,full_range=FALSE) {
    cohort <- if(interaction) phos_interactions else phos
    full_cohort <- cohort
    if(interaction && !full_range) cohort <- anti_join(cohort,select(omitted,site_id,contrast),by=c("site_id","contrast"))
    panel_labels <- if(interaction) short_interaction[match(which,interactions)] else short_mutant[match(which,mutants)]
    d <- cohort |> filter(contrast %in% which) |>
      mutate(y=-log10(FDR),panel=factor(contrast,levels=which,labels=panel_labels),
        selection=factor(case_when(selected_hit & diff<0 ~ "Negative",selected_hit & diff>0 ~ "Positive",TRUE~"Not selected"),
          levels=c("Not selected","Negative","Positive")),
        model=factor(if_else(estimate_type=="observed","Observed","LOD-imputed"),levels=c("Observed","LOD-imputed"))) |>
      arrange(selected_hit,source_observed)
    # Shared axes in the combined and individual panels. Only the explicitly
    # identified Pyk S538 estimate is omitted from the interaction poster view.
    xmax <- max(2,ceiling(max(abs(full_cohort$diff)))); ymax <- max(2,ceiling(max(-log10(cohort$FDR))))
    p <- ggplot(d,aes(diff,y,colour=selection,shape=model)) +
      geom_vline(xintercept=c(-1,1),linetype="dashed",colour="#8A97A6",linewidth=.35) +
      geom_hline(yintercept=1,linetype="dashed",colour="#8A97A6",linewidth=.35) +
      geom_point(size=.95,alpha=.75,stroke=.45) +
      scale_colour_manual(values=c("Not selected"="#B4BEC9","Negative"="#0072B2","Positive"="#D55E00"),drop=FALSE) +
      scale_shape_manual(values=c("Observed"=16,"LOD-imputed"=4),drop=FALSE) +
      scale_x_continuous(limits=c(-xmax,xmax),breaks=pretty(c(-xmax,xmax),5)) +
      scale_y_continuous(limits=c(0,ymax),expand=expansion(mult=c(0,.03))) +
      facet_wrap(~panel,ncol=if(individual)1 else if(interaction)3 else 2) +
      labs(x=if(interaction)"Interaction effect (log2)" else "log2(mutant / ancestor)")
    title <- if(interaction) "Phosphosites: PASN response interactions" else "Phosphosites: mutant-ancestor comparisons"
    subtitle <- if(interaction && !full_range) {
      "FDR < 0.10, |effect| > 1 | One high LOD estimate omitted: Pyk S538"
    } else if(interaction) "FDR < 0.10 and |log2 effect| > 1 | Full range: all estimates shown" else {
      "Within-medium tests | FDR < 0.10 and |log2 effect| > 1"
    }
    p <- style(p,title,subtitle,individual)
    # Separate panels may be used without their combined subtitle. Attach the
    # omission notice to the affected panel itself, at the same 24-pt text size.
    if(interaction && individual && which==interactions[1] && !full_range) {
      p <- p + labs(subtitle="Pyk S538 (LOD) omitted above range")
    }
    metric <- draw_metrics(p,w,h)
    helpers$poster_half_width <- max(.065*2*xmax,metric$label_w/(2*metric$panel_w)*2*xmax)
    helpers$poster_half_height <- max(.065*ymax,metric$label_h/(2*metric$panel_h)*ymax)
    message("Panel mm: ",round(metric$panel_w,1)," x ",round(metric$panel_h,1),
      "; padded label mm: ",round(metric$label_w,1)," x ",round(metric$label_h,1))
    named <- inner_join(d,targets,by="site_id",relationship="many-to-one")
    named <- bind_rows(lapply(which,function(k) {
      message("Place labels: ",k)
      helpers$place_site_labels(filter(d,contrast==k),filter(named,contrast==k),xmax,ymax)
    }))
    before <- ggplot_build(p)
    labelled <- p +
      geom_segment(data=named,aes(x=arrow_x,y=arrow_y,xend=diff,yend=y),inherit.aes=FALSE,
        colour="#34495B",linewidth=.3,arrow=grid::arrow(length=grid::unit(1.35,"mm"),angle=22,type="open"),show.legend=FALSE) +
      geom_label(data=named,aes(x=label_x,y=label_y,label=arrow_label_plot),inherit.aes=FALSE,
        colour="#263A49",fill="white",size=24/.pt,linewidth=0,lineheight=.9,
        label.padding=grid::unit(.1,"lines"),label.r=grid::unit(0,"lines"),show.legend=FALSE)
    after <- ggplot_build(labelled)
    stopifnot(identical(before$data,after$data[1:3]),identical(axis_signature(p),axis_signature(labelled)),
              all(named$n_points_under_label==0),nrow(named)==length(which)*3)
    tag <- if(interaction) {
      paste0("phospho_",if(individual) which else "three_interactions",if(full_range)"_full_range" else "_poster")
    } else if(individual) which else "four_mutant_comparisons"
    n_omitted <- nrow(filter(full_cohort,contrast %in% which))-nrow(d)
    # Check exact source coordinates separately from the disclosed omission.
    stopifnot(identical(after$data[[3]]$x,d$diff),identical(after$data[[3]]$y,d$y))
    audits[[tag]] <<- tibble(figure=tag,original_point_layers_unchanged=n_omitted==0L,
      numeric_axes_unchanged=!(interaction && !full_range),source_estimates_unchanged=TRUE,display_omissions=n_omitted,
      observed_points=sum(d$source_observed),lod_points=sum(!d$source_observed),width_mm=w,height_mm=h,min_font_pt=24)
    label_reports[[tag]] <<- named |> transmute(figure=tag,contrast,arrow_label,diff,FDR,y,label_x,label_y,
      label_half_width,label_half_height,n_points_under_label)
    labelled
  }
  make_total <- function(which=interactions,individual=FALSE) {
    labels <- short_interaction[match(which,interactions)]
    base <- plot_step04_volcano(total,which,labels)
    # Separate panel files retain the limits and ticks of the full Figure 3.
    # This lets them be assembled on the poster without changing visual scales.
    reference_axes <- axis_signature(plot_step04_volcano(total,interactions,short_interaction))[[1]]
    if(individual) {
      full_panel <- ggplot_build(plot_step04_volcano(total,interactions,short_interaction))$layout$panel_params[[1]]
      base <- base +
        scale_x_continuous(limits=full_panel$x$limits,breaks=full_panel$x$breaks[is.finite(full_panel$x$breaks)]) +
        scale_y_continuous(limits=full_panel$y$limits,breaks=full_panel$y$breaks[is.finite(full_panel$y$breaks)])
      a <- axis_signature(base)[[1]]
      stopifnot(isTRUE(all.equal(a$x,reference_axes$x)),isTRUE(all.equal(a$y,reference_axes$y)),
        identical(a$xb[is.finite(a$xb)],reference_axes$xb[is.finite(reference_axes$xb)]),
        identical(a$yb[is.finite(a$yb)],reference_axes$yb[is.finite(reference_axes$yb)]))
    }
    p <- base + scale_colour_manual(values=c("Not selected"="grey75","Selected, negative effect"="#0072B2","Selected, positive effect"="#D55E00"),
      labels=c("Not selected","Negative","Positive"),drop=FALSE) +
      scale_shape_manual(values=c("Source observed"=16,"LOD-flagged model"=4,"Other / unknown"=1),
        labels=c("Source observed"="Observed","LOD-flagged model"="LOD-imputed","Other / unknown"="Other"),drop=TRUE) +
      labs(x="Interaction effect (log2)")
    p <- style(p,"Total proteome: PASN response interactions",
      "(A_TSB - A_SN) - (B_TSB - B_SN) | FDR < 0.10, |effect| > 1",individual)
    stopifnot(isTRUE(all.equal(ggplot_build(base)$data,ggplot_build(p)$data)),identical(axis_signature(base),axis_signature(p)))
    tag <- if(individual) paste0("total_",which) else "three_total_interactions"
    audits[[tag]] <<- tibble(figure=tag,original_point_layers_unchanged=TRUE,numeric_axes_unchanged=TRUE,
      observed_points=sum(p$data$estimate_type=="observed"),lod_points=sum(p$data$estimate_type=="lod_imputed"),
      width_mm=if(individual)171.111 else 529.333,height_mm=115,min_font_pt=24)
    p
  }
  cache_path <- file.path(cache_dir,"compact_plots.rds")
  cached <- if(file.exists(cache_path)) readRDS(cache_path) else list()
  if("mutants" %in% sections) {
    cached$phospho <- make_phos()
    save_plot(cached$phospho,file.path(phos_dir,"phospho_mutants_TSB_SN_JE2_all_models_compact"),529.333,190)
    for(k in mutants) save_plot(make_phos(k,171.111,115,TRUE),file.path(phos_dir,paste0("JE2_all_models_",k,"_compact")),171.111,115)
    # Measure the shared legend on a Cairo device with the actual Arial font.
    grDevices::svg(file.path(cache_dir,"legend_measure.svg"),width=529.333/25.4,height=30/25.4,family="Arial")
    gt <- ggplotGrob(cached$phospho); legend <- gt$grobs[[which(gt$layout$name=="guide-box-bottom")]]
    grid::grid.newpage();grid::grid.draw(legend)
    legend_w <- grid::convertWidth(grid::grobWidth(legend),"mm",valueOnly=TRUE)+2
    legend_h <- grid::convertHeight(grid::grobHeight(legend),"mm",valueOnly=TRUE)+2
    grDevices::dev.off()
    save_plot(legend,file.path(phos_dir,"shared_volcano_legend"),legend_w,legend_h)
  }
  replace_canonical <- function(source_dir,target_dir,name,archive_name) {
    archive <- file.path(REPO_ROOT,"90_archive/figures",archive_name)
    dir.create(archive,recursive=TRUE,showWarnings=FALSE)
    for(ext in c("svg","png","pdf")) {
      target <- file.path(target_dir,paste0(name,".",ext))
      backup <- file.path(archive,basename(target))
      if(file.exists(target) && !file.exists(backup)) stopifnot(file.copy(target,backup,overwrite=FALSE))
      stopifnot(file.copy(file.path(source_dir,basename(target)),target,overwrite=TRUE))
    }
  }
  if("total" %in% sections) {
    cached$total <- make_total()
    save_plot(cached$total,file.path(compact_total_dir,"figure03_total_interaction_volcanoes"),529.333,115)
    for(k in interactions) save_plot(make_total(k,TRUE),file.path(compact_total_dir,paste0("figure03_",k,"_compact")),171.111,115)
    replace_canonical(compact_total_dir,total_dir,"figure03_total_interaction_volcanoes","step04_before_compact_volcanoes")
    svg_paths <- list.files(total_dir,pattern="^figure0[2-6].*\\.svg$",full.names=TRUE)
    if(requireNamespace("zip",quietly=TRUE)) zip::zipr(file.path(total_dir,"step04_poster_SVGs.zip"),files=svg_paths,root=total_dir)
  }
  if("phospho_interactions" %in% sections) {
    cached$phospho_interactions <- make_phos(interactions,529.333,115,interaction=TRUE)
    cached$phospho_interactions_full_range <- make_phos(interactions,529.333,115,interaction=TRUE,full_range=TRUE)
    save_plot(cached$phospho_interactions,file.path(interaction_dir,"phospho_interactions_all_models"),529.333,115)
    save_plot(cached$phospho_interactions_full_range,file.path(interaction_dir,"phospho_interactions_all_models_full_range"),529.333,115)
    for(k in interactions) save_plot(make_phos(k,171.111,115,TRUE,interaction=TRUE),
      file.path(interaction_dir,paste0("all_models_",k,"_compact")),171.111,115)
    replace_canonical(interaction_dir,dirname(interaction_dir),"phospho_interactions_all_models","step07_before_compact_volcanoes")
    write_tsv(omitted |> select(contrast,site_id,site_label,uniprot_acc,diff,FDR,estimate_type,
      best_loc_prob,group_detection_pattern) |> mutate(reason="Poster display only: sole estimate above y=5; retained in full-range figure"),
      file.path(report_dir,"phospho_interaction_display_omission.tsv"))
    interaction_counts <- phos_interactions |> group_by(contrast,estimate_type) |>
      summarise(source_n=n(),selected=sum(selected_hit),.groups="drop") |>
      mutate(display_omissions=as.integer(contrast==interactions[1] & estimate_type=="lod_imputed"),
        poster_n=source_n-display_omissions)
    write_tsv(interaction_counts,file.path(report_dir,"phospho_interaction_population_and_hits.tsv"))
  }
  counts <- phos |> group_by(contrast,estimate_type) |> summarise(n=n(),significant=sum(selected_hit),.groups="drop")
  write_tsv(counts,file.path(report_dir,"phosphosite_population_and_hits.tsv"))
  update_report <- function(rows,name,quote="none") {
    path <- file.path(report_dir,name)
    old <- if(file.exists(path)) read_tsv(path,show_col_types=FALSE) else tibble()
    if(nrow(old)) old <- filter(old,!figure %in% rows$figure)
    write_tsv(bind_rows(old,rows),path,quote=quote)
  }
  update_report(bind_rows(audits),"preservation_checks.tsv")
  update_report(bind_rows(label_reports),"site_label_placement.tsv","needed")
  bundle_svgs <- c(list.files(phos_dir,pattern="\\.svg$",full.names=TRUE),
    list.files(compact_total_dir,pattern="\\.svg$",full.names=TRUE),
    list.files(interaction_dir,pattern="\\.svg$",full.names=TRUE))
  if(requireNamespace("zip",quietly=TRUE)) zip::zipr(file.path(report_dir,"compact_volcano_SVGs.zip"),files=bundle_svgs)
  saveRDS(cached,cache_path)
  print(counts);message("Compact matched SVG/PNG/PDF figures exported; source estimates unchanged; display omissions documented.")
  invisible(cached)
}
compact_volcano_results <- compact_poster_volcanoes()
