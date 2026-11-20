# Compact experimental-design schematic; each dot is a real sample-sheet row.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
samples<-read_samples()
design<-samples |> mutate(genotype=as.character(genotype),medium=as.character(medium)) |>
  count(genotype,medium,name="replicates")
stopifnot(nrow(samples)==27L,!anyDuplicated(samples$sample),nrow(design)==9L,
  all(design$replicates==3L),setequal(design$genotype,c("anc","stk","stp")),
  setequal(design$medium,c("TSB","SNpulse","SN")))
interactions<-step04_contrast_plan() |> filter(contrast_type=="interaction")
stopifnot(setequal(interactions$contrast,c("interaction_stk_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_anc_at_TSB_vs_SN","interaction_stp_vs_stk_at_TSB_vs_SN")))
out<-file.path(REPO_ROOT,"steps/01_design/figures/design")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

draw_design_summary<-function() {
  bg<-"#0B1727"; header<-"#1C3553"; cell<-"#14283F"
  ink<-"#FFFFFF"; secondary<-"#ADBDCF"; accent<-"#36C1B0"
  grid::grid.newpage()
  grid::grid.rect(gp=grid::gpar(fill=bg,col=NA))
  txt<-function(label,x,y,size=12,col=ink,bold=FALSE,just="left") {
    grid::grid.text(label,x=x,y=y,just=just,
      gp=grid::gpar(fontfamily="sans",fontsize=size,col=col,fontface=if(bold) "bold" else "plain"))
  }
  card<-function(x,y,w,h,fill=cell) grid::grid.roundrect(x=x,y=y,width=w,height=h,
    r=grid::unit(7,"pt"),gp=grid::gpar(fill=fill,col=NA))
  txt("JE2 experimental design",.04,.945,22,bold=TRUE)
  txt("3 genotypes  x  3 conditions  x  3 replicates  =  27 samples",.04,.889,12,col=secondary)
  # Three genotype columns and three exposure rows.
  cx<-c(.26,.42,.58); cy<-c(.62,.425,.23)
  names<-c("JE2 ancestor","JE2 stk mutant","JE2 stp mutant")
  roles<-c("wild type","kinase","phosphatase")
  for(j in 1:3) {
    card(cx[j],.786,.148,.117,header)
    txt(names[j],cx[j],.804,14,bold=TRUE,just="centre")
    txt(roles[j],cx[j],.764,11.5,col=secondary,just="centre")
  }
  media<-c("TSB","SNpulse","SN")
  descriptions<-c("control","15 min PASN","~15 h PASN")
  for(i in 1:3) {
    card(.106,cy[i],.132,.167,header)
    txt(media[i],.106,cy[i]+.02,16,bold=TRUE,just="centre")
    txt(descriptions[i],.106,cy[i]-.029,11.5,col=secondary,just="centre")
    for(j in 1:3) {
      card(cx[j],cy[i],.148,.167)
      n<-design$replicates[design$genotype==c("anc","stk","stp")[j] & design$medium==media[i]]
      for(k in seq_len(n)) grid::grid.circle(x=cx[j]+c(-.025,0,.025)[k],y=cy[i]+.018,
        r=grid::unit(.014,"snpc"),gp=grid::gpar(fill=ink,col=NA))
      txt("3 replicates",cx[j],cy[i]-.044,11.5,col=secondary,just="centre")
    }
  }
  grid::grid.lines(x=c(.687,.687),y=c(.145,.844),gp=grid::gpar(col="#2B4058",lwd=1))
  txt("Three focal interactions",.718,.814,17,bold=TRUE)
  txt("Genotype x medium: TSB vs SN",.718,.765,11.7,col=accent)
  pairs<-c("stk vs ancestor","stp vs ancestor","stp vs stk")
  descriptions<-c("Kinase mutant vs ancestor","Phosphatase mutant vs ancestor","Phosphatase vs kinase mutant")
  for(i in 1:3) {
    card(.844,cy[i],.252,.143)
    txt(sprintf("%02d",i),.737,cy[i]+.016,14,col=accent,bold=TRUE)
    txt(pairs[i],.777,cy[i]+.018,15.5,bold=TRUE)
    txt(descriptions[i],.777,cy[i]-.026,10.3,col=secondary)
  }
  txt("Each dot = one replicate.",.04,.092,11,col=secondary)
  txt("PASN = cell-free P. aeruginosa PAO1 supernatant.",.04,.046,10.8,col=secondary)
  txt("SNpulse is outside these focal interactions.",.718,.092,10.5,col=secondary)
}
png_path<-file.path(out,"JE2_27_samples_interactions.png")
grDevices::png(png_path,width=2400,height=1050,res=180,type="cairo")
draw_design_summary(); grDevices::dev.off()
grDevices::svg(file.path(out,"JE2_27_samples_interactions.svg"),width=13.333,height=5.833,bg="#0B1727")
draw_design_summary(); grDevices::dev.off()
grDevices::pdf(file.path(out,"JE2_27_samples_interactions.pdf"),width=13.333,height=5.833,
  bg="#0B1727",useDingbats=FALSE)
draw_design_summary(); grDevices::dev.off()
message("Design summary: nine groups, 27 dots and three verified TSB-vs-SN interactions. ",png_path)
