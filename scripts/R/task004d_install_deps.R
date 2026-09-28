#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else "/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas"
lib <- file.path(root, ".task004d_Rlib")
fallback_lib <- file.path(root, ".task004b_Rlib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, fallback_lib, .libPaths()))
Sys.setenv(PATH = paste("/usr/bin", file.path(root, "tmp", "r_env", "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep))
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(download.file.method = "curl", timeout = 1800)

if (getRversion() < "4.3.0" || getRversion() >= "4.4.0") {
  stop("This installer is pinned to the server's R 4.3 / Bioconductor 3.18 compatibility lane; found ", R.version.string)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, quiet = FALSE)
}
bioc_version <- as.character(BiocManager::version())
if (!identical(bioc_version, "3.18")) {
  stop("R 4.3.3 is expected to use Bioconductor 3.18; BiocManager reports ", bioc_version)
}

# Bioconductor 3.18 is the release matched to the server's R 4.3 runtime.
# Its scDblFinder release predates the xgboost >=3.1 compatibility updates.
# Pin the last pre-3.x CRAN xgboost release in this isolated task library so
# dependency resolution cannot silently pair old scDblFinder with xgboost 3.x.
xgb_version <- if (requireNamespace("xgboost", quietly = TRUE)) {
  packageVersion("xgboost")
} else {
  package_version("0.0.0")
}
if (xgb_version < package_version("3.1.0")) {
  if (xgb_version < package_version("1.7.11.1")) {
    xgb_archive <- "https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.11.1.tar.gz"
    install.packages(xgb_archive, repos = NULL, type = "source", lib = lib, quiet = FALSE)
  }
}

bioc_base <- "https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.18"
options(repos = c(
  BioCsoft = paste0(bioc_base, "/bioc"),
  BioCann = paste0(bioc_base, "/data/annotation"),
  BioCexp = paste0(bioc_base, "/data/experiment"),
  BioCworkflows = paste0(bioc_base, "/workflows"),
  CRAN = "https://cloud.r-project.org"
))
install.packages(
  c("scDblFinder", "SingleCellExperiment", "BiocParallel"),
  lib = lib, dependencies = c("Depends", "Imports", "LinkingTo"),
  Ncpus = 4L, quiet = FALSE
)

.libPaths(c(lib, fallback_lib, .libPaths()))
required <- c("scDblFinder", "SingleCellExperiment", "BiocParallel", "xgboost")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Required dependencies missing after installation: ", paste(missing, collapse = ", "))

if (packageVersion("xgboost") >= package_version("3.1.0") &&
    packageVersion("scDblFinder") < package_version("1.24.8")) {
  stop("Incompatible package pair: scDblFinder >=1.24.8 is required with xgboost >=3.1")
}

cat("R: ", R.version.string, "\n", sep = "")
cat("Bioconductor: ", bioc_version, "\n", sep = "")
for (pkg in required) cat(pkg, " ", as.character(packageVersion(pkg)), "\n", sep = "")

# Smoke-test both the cluster-based and random artificial-doublet code paths
# with synthetic sparse counts before loading the project object for scoring.
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
if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(SummarizedExperiment::colData(sce)))) {
  stop("scDblFinder cluster-mode smoke test did not return score/class")
}
cat("Cluster-mode smoke test: PASSED\n")

set.seed(40401)
sce_small <- SingleCellExperiment(assays = list(counts = smoke_counts[, seq_len(250L), drop = FALSE]))
sce_small <- scDblFinder::scDblFinder(
  sce_small, clusters = NULL, dbr = NULL, dbr.per1k = 0.008, dbr.sd = NULL,
  iter = 2, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE,
  nfeatures = 1352, dims = 20
)
if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(SummarizedExperiment::colData(sce_small)))) {
  stop("scDblFinder random-mode smoke test did not return score/class")
}
cat("Random-mode smoke test: PASSED\n")
