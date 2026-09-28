#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas"
source_qs <- file.path(root, "objects", "merge", "HCC_TA_8datasets_merged_review_v1.qs")
results_dir <- file.path(root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(source_qs)) stop("Approved Task 004d input is absent: ", source_qs)

dim.task004b_chunked_counts <- function(x) c(length(x$feature_names), length(x$cell_names))
dimnames.task004b_chunked_counts <- function(x) list(x$feature_names, x$cell_names)

message("Reading approved input metadata from ", source_qs)
obj <- qs::qread(source_qs, use_alt_rep = FALSE, nthreads = 8L)
if (!inherits(obj, "Seurat") || !"RNA" %in% SeuratObject::Assays(obj)) {
  stop("Approved input is not the expected Seurat object")
}
if (!identical(SeuratObject::Layers(obj[["RNA"]]), "counts")) {
  stop("Expected exactly one RNA counts layer")
}
counts <- obj[["RNA"]]@layers[["counts"]]
if (!inherits(counts, "task004b_chunked_counts")) stop("Expected chunked raw counts")

md <- as.data.table(obj[[]])
required <- c(
  "dataset", "sample_id", "patient_id", "tissue", "project_sample_id",
  "project_patient_id", "nCount_RNA", "nFeature_RNA", "percent.mt",
  "project_broad_celltype"
)
missing_columns <- setdiff(required, names(md))
if (length(missing_columns)) stop("Required metadata columns missing: ", paste(missing_columns, collapse = ", "))

cell_ids <- colnames(obj)
feature_ids <- rownames(obj)
if (anyDuplicated(cell_ids)) stop("Cell identifiers are not globally unique")
if (anyDuplicated(feature_ids)) stop("HGNC feature symbols are not unique")
if (length(cell_ids) != 1490852L || length(feature_ids) != 38025L) {
  stop("Approved input dimensions differ from Task 004c validation: ",
       length(cell_ids), " cells x ", length(feature_ids), " features")
}
if (length(counts$cell_names) != length(cell_ids) || !identical(counts$cell_names, cell_ids)) {
  stop("Counts-layer cell names/order differ from the Seurat object")
}

missingness <- rbindlist(lapply(names(md), function(nm) {
  x <- md[[nm]]
  missing <- is.na(x)
  if (is.character(x) || is.factor(x)) missing <- missing | !nzchar(trimws(as.character(x)))
  data.table(field = nm, class = paste(class(x), collapse = ";"),
             n_missing = sum(missing), fraction_missing = mean(missing))
}))
fwrite(missingness, file.path(results_dir, "task004d_input_metadata_missingness.csv"))

required_missing <- missingness[field %in% required, n_missing]
if (any(required_missing > 0L)) {
  bad <- missingness[field %in% required & n_missing > 0L]
  stop("Missing values in required metadata: ", paste(bad$field, bad$n_missing, collapse = "; "))
}

md[, candidate_neutrophil := tolower(trimws(as.character(project_broad_celltype))) == "neutrophil"]
capture_qc <- md[, .(
  dataset = as.character(dataset[[1L]]),
  sample_id = as.character(sample_id[[1L]]),
  patient_id = as.character(project_patient_id[[1L]]),
  tissue = as.character(tissue[[1L]]),
  n_cells = .N,
  n_datasets = uniqueN(dataset),
  n_sample_ids = uniqueN(sample_id),
  n_patients = uniqueN(project_patient_id),
  n_tissues = uniqueN(tissue),
  n_candidate_neutrophils = sum(candidate_neutrophil)
), by = .(project_sample_id)]
if (any(capture_qc$n_datasets != 1L | capture_qc$n_sample_ids != 1L |
        capture_qc$n_patients != 1L | capture_qc$n_tissues != 1L)) {
  stop("At least one project_sample_id spans conflicting dataset/sample/patient/tissue metadata")
}
if (uniqueN(md$dataset) != 8L || uniqueN(md$project_sample_id) != 194L ||
    uniqueN(md$project_patient_id) != 132L) {
  stop("Dataset/sample/patient counts differ from the approved Task 004c input")
}
fwrite(capture_qc, file.path(results_dir, "task004d_preflight_by_sample.csv"))

dataset_qc <- md[, .(
  n_cells = .N,
  n_samples = uniqueN(project_sample_id),
  n_patients = uniqueN(project_patient_id),
  n_candidate_neutrophils = sum(candidate_neutrophil),
  n_zero_count_metadata = sum(nCount_RNA <= 0),
  n_zero_feature_metadata = sum(nFeature_RNA <= 0)
), by = .(dataset, tissue)]
fwrite(dataset_qc, file.path(results_dir, "task004d_preflight_by_dataset_tissue.csv"))

sha <- system2("sha256sum", source_qs, stdout = TRUE)
sha <- strsplit(sha, "\\s+")[[1L]][[1L]]
input_record <- data.table(
  source_qs = source_qs,
  source_size_bytes = as.numeric(file.info(source_qs)$size),
  source_sha256 = sha,
  n_cells = ncol(obj),
  n_features = nrow(obj),
  n_datasets = uniqueN(md$dataset),
  n_samples = uniqueN(md$project_sample_id),
  n_patients = uniqueN(md$project_patient_id),
  n_candidate_neutrophils = sum(md$candidate_neutrophil),
  counts_layer = "RNA/counts",
  counts_representation = class(counts)[[1L]],
  software_R = R.version.string,
  qs_version = as.character(packageVersion("qs")),
  Seurat_version = as.character(packageVersion("Seurat")),
  SeuratObject_version = as.character(packageVersion("SeuratObject"))
)
fwrite(input_record, file.path(results_dir, "task004d_input_preflight.csv"))
message("Preflight passed: ", ncol(obj), " cells; ", nrow(obj),
        " features; ", uniqueN(md$dataset), " datasets; ",
        uniqueN(md$project_sample_id), " samples; ",
        sum(md$candidate_neutrophil), " preliminary neutrophils.")
