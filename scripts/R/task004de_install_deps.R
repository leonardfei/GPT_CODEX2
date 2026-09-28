#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[[1]],mustWork=TRUE) else normalizePath(".",mustWork=TRUE)
lib <- file.path(root,".task004de_Rlib")
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,file.path(root,".task004b_Rlib"),.libPaths()))
if(!requireNamespace("BiocManager",quietly=TRUE)){
  install.packages("BiocManager",repos="https://cloud.r-project.org",lib=lib)
}
pkgs <- c("scDblFinder","SingleCellExperiment","SingleR","BiocParallel")
need <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(need)) BiocManager::install(need,ask=FALSE,update=FALSE,lib=lib)
bad <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(bad)) stop("Missing packages after install: ",paste(bad,collapse=", "))
cat("Task004d/e R dependencies ready in ",lib,"\n",sep="")
