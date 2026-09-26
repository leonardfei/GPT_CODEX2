#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name); hit <- which(args == key)
  if (length(hit) == 0L || hit == length(args)) return(default)
  args[[hit + 1L]]
}
project_root <- normalizePath(arg_value("project-root", "."), mustWork = TRUE)
lib_dir <- arg_value("lib", file.path(project_root, ".task004b_Rlib"))
.libPaths(c(lib_dir, .libPaths()))
cohort_summary_path <- arg_value("cohort-summary", file.path(project_root, "results", "task004b_merge_review_cohort_summary.csv"))
out_dir <- arg_value("out-dir", file.path(project_root, "tmp", "task004b_h5ad_parts"))
manifest_out <- arg_value("manifest", file.path(project_root, "results", "task004b_h5ad_parts_manifest.csv"))
if (!requireNamespace("anndataR", quietly = TRUE) || !requireNamespace("rhdf5", quietly = TRUE)) stop("Missing anndataR/rhdf5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cs <- fread(cohort_summary_path)
expected_datasets <- c("GSE282701","GSE242889","GSE326201","GSE149614","GSE299340","CRA002308","nature_xue","in_house")
if (!setequal(cs$dataset, expected_datasets)) stop("Cohort summary dataset mismatch")
add_project_ids <- function(md, d) {
  md$project_sample_id <- paste(d, as.character(md$sample_id), sep="::")
  md$project_patient_id <- paste(d, as.character(md$patient_id), sep="::")
  raw_pair <- if ("paired_id" %in% colnames(md)) as.character(md$paired_id) else rep(NA_character_, nrow(md))
  valid <- !is.na(raw_pair) & nzchar(raw_pair) & !(tolower(raw_pair) %in% c("unknown","na","none"))
  md$project_paired_id <- NA_character_
  md$project_paired_id[valid] <- paste(d, raw_pair[valid], sep="::")
  md
}
meta_kind <- function(x) {
  if (is.factor(x) || is.character(x)) return("character")
  if (is.double(x) || is.numeric(x)) return("numeric")
  if (is.integer(x)) return("integer")
  if (is.logical(x)) return("logical")
  "character"
}
field_kinds <- list(); all_fields <- character()
for (d in expected_datasets) {
  p <- cs[dataset == d, cohort_object_path][[1]]
  o <- readRDS(p); md <- add_project_ids(o[[]], d)
  all_fields <- union(all_fields, colnames(md))
  for (nm in colnames(md)) field_kinds[[nm]] <- unique(c(field_kinds[[nm]], meta_kind(md[[nm]])))
  rm(o, md); gc(verbose=FALSE)
}
target_kind <- vapply(all_fields, function(nm) {
  k <- field_kinds[[nm]]
  if ("character" %in% k) return("character")
  if ("numeric" %in% k) return("numeric")
  if ("integer" %in% k && "logical" %in% k) return("numeric")
  if ("integer" %in% k) return("integer")
  "logical"
}, character(1L)); names(target_kind) <- all_fields
coerce_field <- function(x, kind, n) {
  if (is.null(x)) return(switch(kind, character=rep(NA_character_,n), numeric=rep(NA_real_,n), integer=rep(NA_integer_,n), logical=rep(NA,n)))
  switch(kind, character=as.character(x), numeric=as.numeric(x), integer=as.integer(x), logical=as.logical(x))
}
rows <- list()
for (d in expected_datasets) {
  p <- cs[dataset == d, cohort_object_path][[1]]
  message("Preparing cohort H5AD: ", d)
  o <- readRDS(p)
  if (!inherits(o,"Seurat") || !"RNA" %in% Assays(o)) stop("Invalid cohort Seurat: ", d)
  if (!identical(SeuratObject::Layers(o[["RNA"]]), "counts")) stop("Cohort RNA layer mismatch: ", d)
  counts <- SeuratObject::LayerData(o[["RNA"]], layer="counts")
  if (!inherits(counts,"sparseMatrix")) counts <- as(counts,"dgCMatrix")
  if (anyDuplicated(rownames(counts))) stop("Duplicate features: ", d)
  old_ids <- colnames(counts); prefix <- paste0(d,"::")
  new_ids <- ifelse(startsWith(old_ids,prefix), old_ids, paste0(prefix,old_ids))
  if (anyDuplicated(new_ids)) stop("Duplicate global cell IDs: ", d)
  colnames(counts) <- new_ids
  md <- add_project_ids(o[[]], d); rownames(md) <- new_ids; n <- nrow(md)
  for (nm in all_fields) md[[nm]] <- coerce_field(if (nm %in% colnames(md)) md[[nm]] else NULL, target_kind[[nm]], n)
  md <- md[, all_fields, drop=FALSE]
  e <- CreateSeuratObject(counts=counts, assay="RNA", project=d, meta.data=md)
  out <- file.path(out_dir, paste0(d,".h5ad")); if (file.exists(out)) unlink(out)
  anndataR::write_h5ad(e, path=out, compression="gzip", chunk_size="auto", mode="w",
    assay_name="RNA", x_mapping="counts", layers_mapping=FALSE, obs_mapping=TRUE,
    var_mapping=FALSE, obsm_mapping=FALSE, varm_mapping=FALSE, obsp_mapping=FALSE,
    varp_mapping=FALSE, uns_mapping=FALSE)
  ad <- anndataR::read_h5ad(out, as="HDF5AnnData", mode="r", backed=TRUE)
  dims <- dim(ad)
  if (dims[[1L]] != ncol(counts) || dims[[2L]] != nrow(counts)) stop("Cohort H5AD dimension mismatch: ", d)
  rm(ad)
  nnz <- length(as(counts,"dgCMatrix")@x)
  rows[[length(rows)+1L]] <- data.table(dataset=d, source_rds=p, h5ad_path=out, n_obs=ncol(counts), n_vars=nrow(counts), nnz=nnz, file_size_bytes=file.info(out)$size)
  rm(o,counts,md,e); gc(verbose=FALSE)
}
m <- rbindlist(rows)
m[, dataset_order := match(dataset, expected_datasets)]
setorder(m, dataset_order); m[, dataset_order := NULL]
if (sum(m$n_obs) != 1490852L) stop("Cohort H5AD total cell mismatch")
fwrite(m, manifest_out)
message("Per-cohort H5AD export complete")
