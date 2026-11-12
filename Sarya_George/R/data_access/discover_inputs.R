# File discovery and provenance without importing measurement contents into chat.
search_course_files <- function(pattern="", root=DATA_ROOT, recursive=FALSE) {
  if(!dir.exists(root)) stop("Source directory does not exist: ",root)
  paths<-list.files(root,pattern=pattern,recursive=recursive,full.names=TRUE,ignore.case=TRUE)
  paths<-paths[!file.info(paths)$isdir]
  info<-file.info(paths)
  tibble(path=normalizePath(paths,winslash="/",mustWork=TRUE),
         filename=basename(paths),bytes=info$size,modified=as.character(info$mtime))
}
input_provenance <- function(paths) {
  if(any(!file.exists(paths))) stop("One or more input files are missing.")
  tibble(file=basename(paths),bytes=file.info(paths)$size,md5=unname(tools::md5sum(paths)))
}
