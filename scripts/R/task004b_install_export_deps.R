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
options(timeout = 900)
ncpus <- max(1L, min(8L, parallel::detectCores(logical = TRUE)))
install_cran_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, lib = lib_dir, dependencies = NA, Ncpus = ncpus)
}
install_cran_if_missing(c("Rcpp", "RApiSerialize", "BH"))
stringfish_version <- if (requireNamespace("stringfish", quietly = TRUE)) {
  as.character(packageVersion("stringfish"))
} else {
  NA_character_
}
if (!identical(stringfish_version, "0.18.0")) {
  stringfish_src <- file.path(project_root, "tmp", "stringfish_0.18.0.tar.gz")
  stringfish_url <- "https://cran.r-project.org/src/contrib/Archive/stringfish/stringfish_0.18.0.tar.gz"
  if (file.exists(stringfish_src) && isTRUE(file.info(stringfish_src)$size > 500000)) {
    install.packages(stringfish_src, repos = NULL, type = "source", lib = lib_dir,
                     dependencies = FALSE, Ncpus = ncpus)
  } else {
    install.packages(stringfish_url, repos = NULL, type = "source", lib = lib_dir,
                     dependencies = FALSE, Ncpus = ncpus)
  }
}
if (!requireNamespace("qs", quietly = TRUE)) {
  qs_src <- file.path(project_root, "tmp", "qs_0.27.3.tar.gz")
  if (file.exists(qs_src) && isTRUE(file.info(qs_src)$size > 1000000)) {
    install.packages(qs_src, repos = NULL, type = "source", lib = lib_dir,
                     dependencies = FALSE, Ncpus = ncpus)
  } else {
    install.packages(
      "https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.3.tar.gz",
      repos = NULL, type = "source", lib = lib_dir, dependencies = FALSE, Ncpus = ncpus
    )
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib_dir, dependencies = NA)
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  BiocManager::install("rhdf5", lib = lib_dir, ask = FALSE, update = FALSE, Ncpus = ncpus)
}
if (!requireNamespace("anndataR", quietly = TRUE)) {
  # anndataR was added after the Bioconductor 3.18 repository used by the
  # server's R 4.3.3 runtime. Version 0.2.0 is the last official release
  # before the package raised its R minimum and uses the rhdf5 backend.
  install.packages(
    "https://github.com/scverse/anndataR/archive/refs/tags/v0.2.0.tar.gz",
    repos = NULL, type = "source", lib = lib_dir, dependencies = FALSE, Ncpus = ncpus
  )
}
required <- c("qs", "anndataR", "rhdf5")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("R export dependency installation incomplete: ", paste(missing, collapse = ", "))
cat("R export environment ready\n")
for (pkg in required) cat(pkg, as.character(utils::packageVersion(pkg)), "\n")
cat("R library:", normalizePath(lib_dir), "\n")
