#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name); hit <- which(args == key)
  if (length(hit) == 0L || hit == length(args)) return(default)
  args[[hit + 1L]]
}
project_root <- normalizePath(arg_value("project-root", "."), mustWork = TRUE)
lib_dir <- arg_value("lib", file.path(project_root, ".task004b_Rlib"))
.libPaths(c(lib_dir, .libPaths()))
checkpoint <- arg_value("checkpoint", file.path(project_root, "objects", "task004_merge_review", "checkpoint", "HCC_TA_8datasets_merged_review_v1_checkpoint.rds"))
out_path <- arg_value("out", file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.qs"))
validation_path <- arg_value("validation", file.path(project_root, "results", "task004b_qs_validation.csv"))
if (!file.exists(checkpoint)) stop("Missing final merge checkpoint: ", checkpoint)
if (!requireNamespace("qs", quietly = TRUE)) stop("R package 'qs' is not available")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
message("Reading validated final merge checkpoint...")
obj <- readRDS(checkpoint)
if (!inherits(obj, "Seurat")) stop("Checkpoint is not a Seurat object")
if (ncol(obj) != 1490852L || nrow(obj) != 68394L) stop("Checkpoint dimensions mismatch")
md <- obj[[]]
if (uniqueN(md$dataset) != 8L || uniqueN(md$project_sample_id) != 194L || uniqueN(md$project_patient_id) != 132L) stop("Checkpoint metadata cardinality mismatch")
if (anyDuplicated(colnames(obj))) stop("Checkpoint duplicate cell IDs")
if (!identical(SeuratObject::Layers(obj[["RNA"]]), "counts")) stop("Checkpoint RNA layer mismatch")
counts_layer_class <- class(obj[["RNA"]]@layers[["counts"]])[[1L]]
nthreads <- max(1L, min(8L, parallel::detectCores(logical = TRUE)))
message("Writing QS: ", out_path)
qs::qsave(obj, out_path, preset = "high", check_hash = TRUE, nthreads = nthreads)
rm(md, obj); gc(verbose = FALSE)
message("Reload-validating QS...")
obj2 <- qs::qread(out_path, use_alt_rep = FALSE, nthreads = nthreads)
if (!inherits(obj2, "Seurat")) stop("QS reload is not Seurat")
if (ncol(obj2) != 1490852L || nrow(obj2) != 68394L) stop("QS reload dimensions mismatch")
md2 <- obj2[[]]
if (uniqueN(md2$dataset) != 8L || uniqueN(md2$project_sample_id) != 194L || uniqueN(md2$project_patient_id) != 132L) stop("QS reload metadata mismatch")
if (anyDuplicated(colnames(obj2))) stop("QS reload duplicate cell IDs")
if (!identical(SeuratObject::Layers(obj2[["RNA"]]), "counts")) stop("QS reload RNA layer mismatch")
reload_counts_class <- class(obj2[["RNA"]]@layers[["counts"]])[[1L]]
fwrite(data.table(
  status="VALIDATED", qs_path=out_path, qs_size_bytes=file.info(out_path)$size,
  qs_package_version=as.character(utils::packageVersion("qs")),
  n_cells=ncol(obj2), n_features=nrow(obj2), n_datasets=uniqueN(md2$dataset),
  n_samples=uniqueN(md2$project_sample_id), n_patients=uniqueN(md2$project_patient_id),
  duplicated_cell_ids=anyDuplicated(colnames(obj2)),
  RNA_layers=paste(SeuratObject::Layers(obj2[["RNA"]]), collapse=";"),
  counts_layer_class_before=counts_layer_class,
  counts_layer_class_after_reload=reload_counts_class
), validation_path)
message("QS VALIDATED: ", out_path)
