#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else
  "/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas"
lib <- file.path(root, ".task004d_Rlib")
fallback_lib <- file.path(root, ".task004b_Rlib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, fallback_lib, .libPaths()))
env_prefix <- file.path(root, "tmp", "r_env")
Sys.setenv(PATH = paste(
  "/usr/bin", file.path(env_prefix, "bin"), Sys.getenv("PATH"),
  sep = .Platform$path.sep
))
Sys.setenv(
  PKG_CONFIG_PATH = paste(file.path(env_prefix, "lib", "pkgconfig"),
                          Sys.getenv("PKG_CONFIG_PATH"), sep = .Platform$path.sep),
  LD_LIBRARY_PATH = paste(file.path(env_prefix, "lib"),
                          Sys.getenv("LD_LIBRARY_PATH"), sep = .Platform$path.sep)
)
options(download.file.method = "curl", timeout = 1800)

# Fail before downloading R source packages if the project Conda prefix lacks
# the native build tools used by systemfonts, Cairo, XML, RCurl, and Rhtslib.
pkg_config <- Sys.which("pkg-config")
if (!nzchar(pkg_config)) {
  stop("Missing native dependency: pkg-config is required in ", env_prefix,
       ". Install pkg-config in the project environment before retrying.")
}
pkg_config_status <- system2(
  pkg_config,
  c("--exists", "cairo", "freetype2", "fontconfig", "libxml-2.0"),
  stdout = FALSE, stderr = FALSE
)
if (!identical(pkg_config_status, 0L)) {
  stop("pkg-config cannot locate cairo, freetype2, fontconfig, and libxml-2.0 ",
       "in the project environment: ", env_prefix)
}
native_files <- c(file.path(env_prefix, "include", "lzma.h"),
                  file.path(env_prefix, "lib", "liblzma.so"))
if (!all(file.exists(native_files))) {
  stop("Missing XZ/liblzma development files in the project environment: ",
       paste(native_files[!file.exists(native_files)], collapse = ", "))
}

if (getRversion() < "4.3.0" || getRversion() >= "4.4.0") {
  stop("This dependency lane requires R 4.3.x / Bioconductor 3.18; found ", R.version.string)
}

# Use the Bioconductor 3.18 repositories directly. BiocManager's config
# endpoint is not reachable from this server. Only software and annotation
# repositories are needed by these packages.
bioc_base <- "https://bioconductor.statistik.tu-dortmund.de/packages/3.18"
options(repos = c(
  BioCsoft = paste0(bioc_base, "/bioc"),
  BioCann = paste0(bioc_base, "/data/annotation"),
  CRAN = "https://cloud.r-project.org"
))

# Bioconductor 3.18's rtracklayer UCSC C sources are incompatible with the
# server's Conda GCC 16 prototype diagnostics. Build this dependency separately
# with Ubuntu's system GCC, then keep the normal project toolchain for others.
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  system_cc <- Sys.which("gcc")
  if (!nzchar(system_cc)) stop("rtracklayer is missing and system gcc is unavailable")
  makevars <- tempfile("task004d-system-gcc-", tmpdir = file.path(root, "tmp"))
  writeLines(paste0("CC = ", system_cc), makevars)
  old_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  Sys.setenv(R_MAKEVARS_USER = makevars)
  tryCatch(
    install.packages(
      "rtracklayer", lib = lib,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = 1L, quiet = FALSE
    ),
    finally = {
      if (is.na(old_makevars)) Sys.unsetenv("R_MAKEVARS_USER") else
        Sys.setenv(R_MAKEVARS_USER = old_makevars)
      unlink(makevars)
    }
  )
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer did not install with the system GCC compatibility lane")
  }
}

if (!requireNamespace("xgboost", quietly = TRUE) ||
    packageVersion("xgboost") < package_version("1.7.11.1")) {
  xgb_archive <- paste0(
    "https://cran.r-project.org/src/contrib/Archive/xgboost/",
    "xgboost_1.7.11.1.tar.gz"
  )
  install.packages(xgb_archive, repos = NULL, type = "source", lib = lib)
}

required <- c("scDblFinder", "SingleCellExperiment", "BiocParallel", "SingleR")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages from Bioconductor 3.18 / CRAN: ",
          paste(missing, collapse = ", "))
  install.packages(
    missing,
    lib = lib,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    Ncpus = 4L,
    quiet = FALSE
  )
}

.libPaths(c(lib, fallback_lib, .libPaths()))
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing dependencies after installation: ", paste(missing, collapse = ", "))
if (packageVersion("xgboost") >= package_version("3.1.0") &&
    packageVersion("scDblFinder") < package_version("1.24.8")) {
  stop("Incompatible pair: scDblFinder >=1.24.8 is required with xgboost >=3.1")
}

cat("R: ", R.version.string, "\n", sep = "")
cat("Bioconductor compatibility lane: 3.18\n")
cat("Task library: ", lib, "\n", sep = "")
for (pkg in c(required, "xgboost")) {
  cat(pkg, " ", as.character(packageVersion(pkg)), "\n", sep = "")
}

# Smoke-test both approved scDblFinder modes before loading the project object.
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(BiocParallel)
})
set.seed(40400)
smoke_counts <- Matrix::Matrix(
  matrix(stats::rpois(1500L * 600L, lambda = 0.12), nrow = 1500L),
  sparse = TRUE
)
colnames(smoke_counts) <- paste0("smoke_", seq_len(ncol(smoke_counts)))
rownames(smoke_counts) <- paste0("gene_", seq_len(nrow(smoke_counts)))
sce <- SingleCellExperiment(assays = list(counts = smoke_counts))
sce <- scDblFinder::scDblFinder(
  sce, clusters = TRUE, dbr = NULL, dbr.per1k = 0.008, dbr.sd = NULL,
  iter = 2, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE,
  nfeatures = 1352, dims = 20
)
stopifnot(all(c("scDblFinder.score", "scDblFinder.class") %in%
                colnames(SummarizedExperiment::colData(sce))))
cat("Cluster-mode smoke test: PASSED\n")

set.seed(40401)
sce_small <- SingleCellExperiment(assays = list(counts = smoke_counts[, seq_len(250L), drop = FALSE]))
sce_small <- scDblFinder::scDblFinder(
  sce_small, clusters = NULL, dbr = NULL, dbr.per1k = 0.008, dbr.sd = NULL,
  iter = 2, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE,
  nfeatures = 1352, dims = 20
)
stopifnot(all(c("scDblFinder.score", "scDblFinder.class") %in%
                colnames(SummarizedExperiment::colData(sce_small))))
cat("Random-mode smoke test: PASSED\n")
