# Load the Sarya_George project from its own folder or the repository root.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2)
})
find_group_root <- function(start=getwd()) {
  path <- normalizePath(start,winslash="/",mustWork=TRUE)
  repeat {
    if(file.exists(file.path(path,"Sarya_George.Rproj"))) return(path)
    if(file.exists(file.path(path,"Sarya_George","Sarya_George.Rproj"))) return(file.path(path,"Sarya_George"))
    parent<-dirname(path)
    if(parent==path) stop("Open Sarya_George.Rproj or run from the repository root.")
    path<-parent
  }
}
# Backward-compatible name used by the existing analysis functions.
REPO_ROOT <- find_group_root()
DATA_ROOT <- Sys.getenv("JE2_DATA_DIR",unset="")
if(!nzchar(DATA_ROOT)) {
  config<-file.path(REPO_ROOT,"config","data_dir.txt")
  DATA_ROOT<-if(file.exists(config)) trimws(readLines(config,n=1,warn=FALSE)) else file.path(dirname(REPO_ROOT),"data")
}
DATA_ROOT<-normalizePath(DATA_ROOT,winslash="/",mustWork=FALSE)
# Recursive loading supports the purpose-specific R subfolders.
for(helper in sort(list.files(file.path(REPO_ROOT,"R"),pattern="\\.R$",full.names=TRUE,recursive=TRUE))) {
  if(basename(helper)!="setup.R") source(helper)
}
rm(helper)
theme_set(theme_bw(base_size=11))
