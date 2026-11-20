# Project workflow overview with the original step numbers; omit steps 2 and 11.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
draw_project_workflow<-function() {
  bg<-"#0B1727"; box<-"#162D46"; ink<-"#FFFFFF"; muted<-"#AEBDCD"; teal<-"#38C4B1"
  grid::grid.newpage(); grid::grid.rect(gp=grid::gpar(fill=bg,col=NA))
  txt<-function(value,x,y,size=17,col=ink,bold=FALSE) grid::grid.text(value,x=x,y=y,just="left",
    gp=grid::gpar(fontfamily="sans",fontsize=size,col=col,fontface=if(bold) "bold" else "plain"))
  txt("JE2 PASN project workflow",.045,.94,27,bold=TRUE)
  txt("Protein abundance, functional modules and phosphorylation",.045,.877,15,col=muted)
  numbers<-c(1,3,4,5,6,7,8,9,10)
  titles<-c("Design & input files","ID mapping & annotation","Three interaction contrasts",
    "Functional enrichment","Protein interaction networks","Phosphosite interpretation",
    "Mutant site overlaps","Residues & sequence context","Biological interpretation")
  details<-c("27 samples; total + phospho proteomes","JE2 IDs, UniProt and annotation terms","stk / anc, stp / anc, stp / stk",
    "GO / KEGG with ORA and GSEA","STRING associations; Stk1 / Stp1","Site localization and parent proteins",
    "Venn diagrams within TSB and SN","pS / pT / pY / pH; logos","Combine modules, sites and candidates")
  positions<-data.frame(x=c(.19,.5,.81,.81,.5,.19,.19,.5,.81),y=c(rep(.72,3),rep(.48,3),rep(.24,3)))
  # Orthogonal connectors make the reading order explicit without crossings.
  for(i in 1:8) {
    a<-positions[i,]; b<-positions[i+1,]
    if(a$y==b$y) {
      direction<-sign(b$x-a$x); ax<-a$x+direction*.137; bx<-b$x-direction*.137
      ay<-by<-a$y
    } else {ax<-bx<-a$x; ay<-a$y-.084; by<-b$y+.084}
    grid::grid.lines(x=c(ax,bx),y=c(ay,by),
      arrow=grid::arrow(type="closed",length=grid::unit(6,"pt")),gp=grid::gpar(col=teal,lwd=1.6,fill=teal))
  }
  for(i in seq_along(numbers)) {
    x<-positions$x[i]; y<-positions$y[i]
    grid::grid.roundrect(x=x,y=y,width=.266,height=.168,r=grid::unit(8,"pt"),gp=grid::gpar(fill=box,col=NA))
    txt(sprintf("%02d",numbers[i]),x-.116,y+.045,14,col=teal,bold=TRUE)
    txt(titles[i],x-.116,y+.003,16,bold=TRUE)
    txt(details[i],x-.116,y-.046,10.8,col=muted)
  }
  txt("Focal interactions compare the TSB-to-SN response between genotypes.",.045,.085,12,col=muted)
}
out<-file.path(REPO_ROOT,"steps/01_design/figures/design"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
grDevices::png(file.path(out,"JE2_project_workflow.png"),width=2400,height=1200,res=180,type="cairo")
draw_project_workflow(); grDevices::dev.off()
grDevices::svg(file.path(out,"JE2_project_workflow.svg"),width=13.333,height=6.667,bg="#0B1727")
draw_project_workflow(); grDevices::dev.off()
message("Workflow saved with steps 1,3,4,5,6,7,8,9,10.")
