# Poster summary of existing observed-only enrichment; no statistical reruns.
source(if(file.exists("Sarya_George.Rproj"))"R/setup.R" else "Sarya_George/R/setup.R")

poster_functional_summary <- function() {
  out <- file.path(REPO_ROOT,"steps/05_enrichment/figures/poster_functional_summary")
  report <- file.path(REPO_ROOT,"steps/05_enrichment/reports/poster_functional_summary")
  for(p in c(out,report))dir.create(p,recursive=TRUE,showWarnings=FALSE)
  contrasts <- filter(step04_contrast_plan(),contrast_type=="interaction")$contrast
  gsea <- read_tsv(file.path(REPO_ROOT,"results/step05/observed_only/gsea_all.tsv.gz"),show_col_types=FALSE) |>
    filter(collection=="KEGG",contrast %in% contrasts) |>
    mutate(significant=is.finite(FDR_study)&FDR_study<.1,column=match(contrast,contrasts))
  ora <- read_tsv(file.path(REPO_ROOT,"steps/05_enrichment/reports/observed_only/ora_all_significant_GO_KEGG.tsv"),show_col_types=FALSE) |>
    filter(collection=="GO_BP",contrast %in% contrasts)
  stopifnot(all(gsea$policy=="observed_only"),all(gsea$background_size==1858),
    sum(gsea$significant)==13L,n_distinct(gsea$term_id[gsea$significant])==12L,
    identical(as.integer(table(factor(gsea$contrast[gsea$significant],levels=contrasts))),c(7L,3L,3L)),
    nrow(ora)==8L,all(ora$contrast==contrasts[1]),all(ora$background_size==1858))
  row_plan <- gsea |> filter(significant) |> arrange(column,desc(NES)) |>
    distinct(term_id,.keep_all=TRUE) |> transmute(term_id,term_description,row=row_number())
  matrix <- gsea |> semi_join(row_plan,by="term_id") |>
    left_join(select(row_plan,term_id,row),by="term_id",relationship="many-to-one") |>
    arrange(row,column)
  stopifnot(nrow(matrix)==36L,all(is.finite(matrix$NES)),!anyDuplicated(matrix[c("term_id","contrast")]))
  # Same contributing protein set identifies nested GO findings, not independent
  # mechanisms. Keep every term in the companion table; display lowest-q per set.
  ora <- ora |> mutate(overlap_group=vapply(strsplit(overlap_ids,";",fixed=TRUE),
    function(ids)paste(sort(unique(ids)),collapse=";"),character(1)))
  reps <- ora |> arrange(FDR_study,term_id) |> group_by(overlap_group) |>
    mutate(n_related_terms=n()) |> slice_head(n=1) |> ungroup() |> arrange(desc(fold_enrichment))
  stopifnot(nrow(reps)==3L,identical(reps$term_id,c("GO:0019835","GO:0009082","GO:0006189")),
    identical(as.integer(reps$n_related_terms),c(1L,2L,5L)))
  # Verify the overrepresentation ratios rather than treating the source labels
  # as aesthetic numbers. No p-values or correction families are changed here.
  stopifnot(all(abs(ora$fold_enrichment-(ora$overlap/ora$foreground_size)/
    (ora$term_size/ora$background_size))<1e-10))
  write_tsv(matrix,file.path(report,"KEGG_matrix_all_36_tests.tsv"))
  write_tsv(filter(matrix,significant),file.path(report,"KEGG_all_13_significant_pairs.tsv"))
  write_tsv(ora,file.path(report,"GO_BP_ORA_all_8_terms.tsv"))
  write_tsv(reps,file.path(report,"GO_BP_ORA_3_representatives.tsv"))
  count_table <- tibble(contrast=contrasts,
    KEGG_GSEA_significant=vapply(contrasts,function(k)sum(gsea$significant[gsea$contrast==k]),integer(1)),
    GO_BP_ORA_significant=vapply(contrasts,function(k)sum(ora$contrast==k),integer(1)))
  write_tsv(count_table,file.path(report,"enrichment_counts.tsv"))

  mm <- function(x)grid::unit(x,"mm")
  text <- function(label,x,y,size=24,bold=FALSE,just="left",col="#172E40",lineheight=1.05) {
    grid::grid.text(label,x=mm(x),y=mm(y),just=just,
      gp=grid::gpar(fontfamily="Arial",fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=lineheight))
  }
  # Green-weak-friendly blue/orange map, consistently keyed to the native sign.
  ramp <- grDevices::colorRamp(c("#70AED4","#FFFFFF","#EFAC7D"))
  nes_colour <- function(x) {
    rgb <- ramp((pmax(-2,pmin(2,x))+2)/4)
    grDevices::rgb(rgb[,1],rgb[,2],rgb[,3],maxColorValue=255)
  }
  w <- 529.333; h <- 250
  measure <- tempfile(fileext=".svg")
  grDevices::svg(measure,width=w/25.4,height=h/25.4,family="Arial")
  label_width <- max(vapply(row_plan$term_description,function(s)
    grid::convertWidth(grid::grobWidth(grid::textGrob(s,gp=grid::gpar(fontfamily="Arial",fontsize=24))),"mm",valueOnly=TRUE),numeric(1)))
  grDevices::dev.off();unlink(measure)
  cell_start <- label_width+9; cell_w <- (350-cell_start)/3
  stopifnot(cell_w>48)
  draw <- function() {
    grid::grid.newpage()
    text("Functional summaries of proteome interactions",4,240,32,TRUE)
    text("1,858 mapped proteins; observed models only; study-wide FDR < 0.10",4,226)
    text("A  KEGG GSEA: 12 pathways, 13 enrichments",4,211,26,TRUE)
    headers <- c("stk vs\nancestor","stp vs\nancestor","stp vs stk")
    for(j in 1:3)text(headers[j],cell_start+(j-.5)*cell_w,194,24,TRUE,"centre",lineheight=.9)
    for(i in seq_len(nrow(row_plan))) {
      y <- 174-(i-1)*12
      text(row_plan$term_description[i],cell_start-4,y,just="right")
      for(j in 1:3) {
        z <- filter(matrix,row==i,column==j)
        grid::grid.rect(x=mm(cell_start+(j-.5)*cell_w),y=mm(y),width=mm(cell_w-1),height=mm(10.8),
          gp=grid::gpar(fill=if(z$significant)nes_colour(z$NES) else "#F1F3F5",col=NA))
        text(if(z$significant)sprintf("%+.2f",z$NES) else "\u2014",
          cell_start+(j-.5)*cell_w,y,24,z$significant,"centre",if(z$significant)"#172E40" else "#667582")
      }
    }
    text("\u2014 FDR \u2265 0.10",4,22)
    text("NES",141,22)
    for(k in 1:100)grid::grid.rect(x=mm(175+(k-.5)*1.5),y=mm(25),width=mm(1.52),height=mm(4),
      gp=grid::gpar(fill=nes_colour(-2+4*(k-.5)/100),col=NA))
    for(k in 0:2)text(c("-2","0","+2")[k+1],175+k*75,14,just="centre")
    grid::grid.lines(x=mm(c(357,357)),y=mm(c(28,216)),gp=grid::gpar(col="#CBD3DB",lwd=1))
    bx <- 365
    text("B  GO-BP ORA",bx,211,26,TRUE)
    text("stk vs ancestor",bx,197,24,TRUE)
    card_labels <- c("Cytolysis","Branched-chain amino-acid\nbiosynthesis","De novo IMP biosynthesis")
    for(i in 1:3) {
      yy <- c(179,132,91)[i]; by <- c(161,120,79)[i]; qy <- c(150,109,68)[i]
      text(card_labels[i],bx,yy,lineheight=.95)
      length <- 122*reps$fold_enrichment[i]/20
      grid::grid.rect(x=mm(bx+length/2),y=mm(by),width=mm(length),height=mm(4),
        gp=grid::gpar(fill=if(reps$effect_sign[i]=="positive")"#D55E00" else "#0072B2",col=NA))
      text(sprintf("%.1f\u00d7",reps$fold_enrichment[i]),bx+length+3,by)
      q <- formatC(reps$FDR_study[i],format="g",digits=2)
      text(sprintf("%s; q=%s; %d proteins",if(reps$effect_sign[i]=="positive")"Positive" else "Negative",q,reps$overlap[i]),bx,qy)
    }
    text("Fold enrichment",bx,57)
    text("stp/ancestor and stp/stk:\nno significant GO-BP ORA terms",bx,41,lineheight=.95)
    text("3 representatives of 8 terms",bx,21)
    text("Positive NES: more positive SN\u2013TSB response in second genotype; negative NES: more positive in first.",4,5)
  }
  stem <- file.path(out,"functional_summary_interactions")
  grDevices::svg(paste0(stem,".svg"),width=w/25.4,height=h/25.4,family="Arial",bg="white")
  draw();grDevices::dev.off()
  svg <- readLines(paste0(stem,".svg"),warn=FALSE); root <- grep("^<svg ",svg)[1]
  svg[root] <- sub('width="[^"]+"',sprintf('width="%.6fmm"',w),svg[root])
  svg[root] <- sub('height="[^"]+"',sprintf('height="%.6fmm"',h),svg[root])
  stopifnot(!any(grepl("<image",svg,fixed=TRUE)))
  writeLines(svg,paste0(stem,".svg"),useBytes=TRUE)
  grDevices::png(paste0(stem,".png"),width=w,height=h,units="mm",res=180,type="cairo",bg="white")
  draw();grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=w/25.4,height=h/25.4,family="Arial",bg="white")
  draw();grDevices::dev.off()
  writeLines(c("Existing enrichment results rendered without refitting.",
    "PASS: all 12 KEGG pathways / 13 significant pairs; all 36 matrix tests retained in table.",
    "PASS: all 8 GO-BP ORA terms retained; 3 representatives chosen by exact hit overlap and smallest q.",
    "PASS: ORA fold enrichment independently checked; 1858-protein observed-only universe.",
    "PASS: SVG vector only, explicit physical size; Arial 24-pt minimum (headings 26/32)."),
    file.path(report,"render_checks.txt"))
  print(count_table)
  invisible(list(matrix=matrix,ora=ora,representatives=reps))
}
functional_summary <- poster_functional_summary()
