#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name); hit <- which(args == key)
  if (length(hit) == 0L || hit == length(args)) return(default)
  args[[hit + 1L]]
}
project_root <- normalizePath(arg_value("project-root", "."), mustWork = TRUE)
lib_dir <- arg_value("lib", file.path(project_root, ".task004b_Rlib"))
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
ncpus <- max(1L, min(8L, parallel::detectCores(logical = TRUE)))
install_cran_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, lib = lib_dir, dependencies = NA, Ncpus = ncpus)
}
install_cran_if_missing(c("Rcpp", "RApiSerialize", "stringfish", "BH"))
if (!requireNamespace("qs", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.3.tar.gz",
    repos = NULL, type = "source", lib = lib_dir, dependencies = FALSE, Ncpus = ncpus
  )
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib_dir, dependencies = NA)
}
if (!requireNamespace("anndataR", quietly = TRUE) || !requireNamespace("rhdf5", quietly = TRUE)) {
  BiocManager::install(c("anndataR", "rhdf5"), lib = lib_dir, ask = FALSE, update = FALSE, Ncpus = ncpus)
}
required <- c("qs", "anndataR", "rhdf5")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("R export dependency installation incomplete: ", paste(missing, collapse = ", "))
cat("R export environment ready\n")
for (pkg in required) cat(pkg, as.character(utils::packageVersion(pkg)), "\n")
cat("R library:", normalizePath(lib_dir), "\n")
