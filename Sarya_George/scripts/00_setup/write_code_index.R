source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
folders<-file.path(REPO_ROOT,c("R","scripts","tests","notebook"))
paths<-unlist(lapply(folders,list.files,pattern="\\.R$",recursive=TRUE,full.names=TRUE),use.names=FALSE)
paths<-c(paths,file.path(REPO_ROOT,"run_step.R"))
index<-bind_rows(lapply(paths,function(path) {
  code_lines<-readLines(path,warn=FALSE)
  names_found<-sub("[[:space:]]*<-.*$","",grep("^[A-Za-z.][A-Za-z0-9._]*[[:space:]]*<-[[:space:]]*function",code_lines,value=TRUE))
  tibble(file=substring(path,nchar(REPO_ROOT)+2),lines=length(code_lines),
         functions=if(length(names_found)) paste(names_found,collapse="; ") else "(entry script / no top-level functions)")
}))
write_tsv(index,file.path(REPO_ROOT,"docs","R_code_index.tsv"))
message(nrow(index)," R files indexed.")
