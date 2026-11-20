source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
paths<-file.path(DATA_ROOT,"derived",c("samples.tsv","protein_annotation.tsv.gz",
  "dea_total.tsv.gz","dea_phospho_STYH.tsv.gz","abundances_total.tsv.gz",
  "abundances_phospho_STYH.tsv.gz","sites_STYH.tsv.gz"))
inventory<-input_provenance(paths)
out<-file.path(REPO_ROOT,"results","input_inventory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv(inventory,file.path(out,"required_inputs.tsv"))
folders<-file.path(DATA_ROOT,c("derived","raw/resources","raw/total","raw/phospho_STYH","raw/sites"))
available<-bind_rows(lapply(folders[dir.exists(folders)],function(folder) {
  search_course_files("\\.(xlsx|tsv|gz|fasta|gff|html)$",root=folder,recursive=FALSE)
}))
write_tsv(available,file.path(out,"available_source_files.tsv"))
print(inventory)
