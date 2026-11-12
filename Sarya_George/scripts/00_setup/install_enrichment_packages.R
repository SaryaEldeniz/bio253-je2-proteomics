# R 4.4 uses the matching Bioconductor 3.20 Windows binary repository.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
lib<-file.path(REPO_ROOT,"renv","library")
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
if(getRversion()>="4.4" && getRversion()<"4.5" && .Platform$OS.type=="windows") {
  repos<-c(BioCsoft="https://bioconductor.org/packages/3.20/bioc",CRAN="https://cloud.r-project.org")
  install.packages(c("fgsea","BiocParallel"),lib=lib,repos=repos,type="win.binary",dependencies=NA)
} else {
  if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",lib=lib,repos="https://cloud.r-project.org")
  BiocManager::install(c("fgsea","BiocParallel"),lib=lib,ask=FALSE,update=FALSE)
}
stopifnot(requireNamespace("fgsea",quietly=TRUE),requireNamespace("BiocParallel",quietly=TRUE))
cat("fgsea",as.character(packageVersion("fgsea")),"is available in the project-local library.\n")
