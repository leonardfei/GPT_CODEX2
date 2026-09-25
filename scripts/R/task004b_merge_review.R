#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name)
  hit <- which(args == key)
  if (length(hit) == 0L || hit == length(args)) return(default)
  args[[hit + 1L]]
}

project_root <- normalizePath(arg_value("project-root", "."), mustWork = TRUE)
validation_csv <- arg_value(
  "validation-csv",
  file.path(project_root, "results", "task004_object_validation.csv")
)
out_path <- arg_value(
  "out",
  file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.rds")
)
cohort_root <- arg_value(
  "cohort-root",
  file.path(project_root, "objects", "task004_merge_review", "cohort_merged")
)
results_root <- arg_value("results-root", file.path(project_root, "results"))
report_path <- arg_value(
  "report",
  file.path(project_root, "reports", "task_004b_merge_review_report.md")
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
dir.create(cohort_root, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

expected_datasets <- c(
  "GSE282701", "GSE242889", "GSE326201", "GSE149614", "GSE299340",
  "CRA002308", "nature_xue", "in_house"
)

manifest <- fread(validation_csv)
required_cols <- c("dataset", "output_path", "annotated_object_cells")
missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols)) {
  stop("Missing required validation columns: ", paste(missing_cols, collapse = ", "))
}
if (!setequal(unique(manifest$dataset), expected_datasets)) {
  stop("Dataset set mismatch. Found: ", paste(sort(unique(manifest$dataset)), collapse = ", "))
}
if (nrow(manifest) != 103L) {
  stop("Expected 103 Task 004 annotated source objects, found ", nrow(manifest))
}
if (any(!file.exists(manifest$output_path))) {
  stop("Missing annotated object(s): ",
       paste(manifest$output_path[!file.exists(manifest$output_path)], collapse = "; "))
}

expected_total <- sum(as.numeric(manifest$annotated_object_cells))
if (expected_total != 1490852L) {
  stop("Unexpected Task 004 annotated cell total: ", expected_total, " (expected 1490852)")
}

required_meta <- c(
  "dataset", "sample_id", "patient_id", "tissue", "abundance_eligible",
  "project_broad_celltype", "neutrophil_confidence",
  "source_author_annotation", "qc_status"
)

get_counts <- function(obj) {
  if (!"RNA" %in% Assays(obj)) stop("Object lacks RNA assay")
  layers <- Layers(obj[["RNA"]])
  if (!"counts" %in% layers) stop("Object lacks RNA counts layer")
  x <- LayerData(obj[["RNA"]], layer = "counts")
  if (ncol(x) != ncol(obj)) stop("RNA counts/object cell mismatch")
  x
}

prepare_for_merge <- function(path, dataset) {
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Not a Seurat object: ", path)
  md <- obj[[]]
  miss <- setdiff(required_meta, colnames(md))
  if (length(miss)) {
    stop("Missing metadata in ", basename(path), ": ", paste(miss, collapse = ", "))
  }

  md$task004_qc_status_premerge <- as.character(md$qc_status)
  if ("rescued_by_corrected_qc" %in% colnames(md)) {
    md$task004_rescued_premerge <- as.logical(md$rescued_by_corrected_qc)
  } else {
    md$task004_rescued_premerge <- NA
  }

  if (dataset %in% c("CRA002308", "in_house")) {
    md$rescued_by_corrected_qc <- FALSE
    md$qc_status <- "corrected_qc_pass"
  } else if (dataset == "nature_xue") {
    md$rescued_by_corrected_qc <- FALSE
    md$qc_status <- "author_processed"
  } else {
    if (!"rescued_by_corrected_qc" %in% colnames(md)) {
      stop("Original five-dataset object lacks rescued_by_corrected_qc: ", path)
    }
  }

  md$annotation_status <- "preliminary_unvalidated_task004"
  md$merge_review_source_object <- path
  md$merge_review_dataset <- dataset

  counts <- get_counts(obj)
  if (any(counts < 0)) stop("Negative counts found: ", path)

  slim <- CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = dataset,
    meta.data = md,
    min.cells = 0,
    min.features = 0
  )
  if (!identical(colnames(slim), colnames(obj))) {
    stop("Cell IDs changed while slimming: ", path)
  }

  rm(obj, counts, md)
  gc(verbose = FALSE)
  slim
}

merge_two <- function(x, y) {
  tryCatch(
    merge(x = x, y = y, merge.data = FALSE, merge.dr = FALSE, collapse = TRUE),
    error = function(e) {
      message("collapse=TRUE merge failed; retrying without collapse: ", conditionMessage(e))
      merge(x = x, y = y, merge.data = FALSE, merge.dr = FALSE)
    }
  )
}

cohort_rows <- list()

for (dataset in expected_datasets) {
  files <- manifest[manifest$dataset == dataset, output_path]
  message("Preparing cohort ", dataset, " from ", length(files), " object(s)")

  acc <- NULL
  for (i in seq_along(files)) {
    message("  [", i, "/", length(files), "] ", basename(files[[i]]))
    slim <- prepare_for_merge(files[[i]], dataset)
    if (is.null(acc)) {
      acc <- slim
    } else {
      acc <- merge_two(acc, slim)
    }
    rm(slim)
    gc(verbose = FALSE)
  }

  cohort_path <- file.path(cohort_root, paste0(dataset, "_merged_counts_metadata.rds"))
  saveRDS(acc, cohort_path, compress = "gzip")
  cohort_rows[[length(cohort_rows) + 1L]] <- data.table(
    dataset = dataset,
    cohort_object_path = cohort_path,
    n_cells = ncol(acc),
    n_features = nrow(acc),
    n_samples = uniqueN(acc$sample_id),
    n_patients = uniqueN(acc$patient_id),
    n_tumor = sum(acc$tissue == "Tumor"),
    n_adjacent = sum(acc$tissue == "Adjacent"),
    n_preliminary_neutrophil = sum(acc$project_broad_celltype == "neutrophil", na.rm = TRUE)
  )
  rm(acc)
  gc(verbose = FALSE)
}

cohort_summary <- rbindlist(cohort_rows)
setorder(cohort_summary, dataset)
fwrite(cohort_summary, file.path(results_root, "task004b_merge_review_cohort_summary.csv"))

if (sum(cohort_summary$n_cells) != expected_total) {
  stop("Cohort-level merge changed total cell count")
}

merged <- NULL
for (dataset in expected_datasets) {
  path <- cohort_summary[cohort_summary$dataset == dataset, cohort_object_path][[1]]
  message("Final merge: ", dataset)
  obj <- readRDS(path)
  if (is.null(merged)) {
    merged <- obj
  } else {
    merged <- merge_two(merged, obj)
  }
  rm(obj)
  gc(verbose = FALSE)
}

if (ncol(merged) != expected_total) {
  stop("Final merged object cell count mismatch: ", ncol(merged), " vs ", expected_total)
}
if (anyDuplicated(colnames(merged))) {
  stop("Final merged object has duplicate cell IDs")
}
if (!setequal(unique(as.character(merged$dataset)), expected_datasets)) {
  stop("Final merged object dataset metadata mismatch")
}

merged@misc$merge_review <- list(
  task = "task_004b",
  purpose = "Eight-cohort unintegrated merge for manual annotation review",
  source = "Task 004 annotated objects",
  integration_performed = FALSE,
  normalization_performed = FALSE,
  reductions_preserved = FALSE,
  graphs_preserved = FALSE,
  task004_annotation_status = "preliminary_unvalidated",
  known_task004_rescued_bug_corrected_for = c("CRA002308", "in_house"),
  n_cells = ncol(merged),
  n_features = nrow(merged),
  datasets = expected_datasets
)

saveRDS(merged, out_path, compress = "gzip")

meta <- merged[[]]
dataset_summary <- as.data.table(meta)[, .(
  n_cells = .N,
  n_samples = uniqueN(sample_id),
  n_patients = uniqueN(patient_id),
  n_tumor = sum(tissue == "Tumor"),
  n_adjacent = sum(tissue == "Adjacent"),
  n_preliminary_neutrophil = sum(project_broad_celltype == "neutrophil", na.rm = TRUE),
  n_rescued_current = sum(rescued_by_corrected_qc %in% TRUE, na.rm = TRUE)
), by = dataset]
setorder(dataset_summary, dataset)
fwrite(dataset_summary, file.path(results_root, "task004b_merge_review_dataset_summary.csv"))

metadata_fields <- data.table(field = colnames(meta))
fwrite(metadata_fields, file.path(results_root, "task004b_merge_review_metadata_fields.csv"))

rna_layers <- Layers(merged[["RNA"]])
validation <- data.table(
  object_path = out_path,
  n_cells = ncol(merged),
  n_features = nrow(merged),
  n_datasets = uniqueN(meta$dataset),
  n_samples = uniqueN(meta$sample_id),
  n_patients = uniqueN(meta$patient_id),
  duplicated_cell_ids = anyDuplicated(colnames(merged)),
  RNA_layers = paste(rna_layers, collapse = ";"),
  has_counts_layer = any(grepl("^counts", rna_layers)),
  has_project_broad_celltype = "project_broad_celltype" %in% colnames(meta),
  has_source_author_annotation = "source_author_annotation" %in% colnames(meta),
  annotation_status = paste(unique(meta$annotation_status), collapse = ";"),
  integration_performed = FALSE
)
fwrite(validation, file.path(results_root, "task004b_merge_review_validation.csv"))

sha <- tryCatch(
  system2("sha256sum", out_path, stdout = TRUE, stderr = TRUE),
  error = function(e) paste("sha256sum unavailable:", conditionMessage(e))
)

report <- c(
  "# Task 004b report — eight-cohort unintegrated Seurat merge for annotation review",
  "",
  "## Purpose",
  "",
  "This object was created for manual inspection of the current Task 004 annotations. No RPCA/Harmony integration, batch correction, normalization, PCA, UMAP, or reclustering was performed.",
  "",
  "## Input",
  "",
  "- 103 Task 004 annotated Seurat objects",
  "- 8 cohorts",
  paste0("- Total cells: ", format(ncol(merged), big.mark = ",")),
  "",
  "## Merge content",
  "",
  "- RNA counts and cell-level metadata were retained.",
  "- Pre-existing reductions/graphs were intentionally discarded because they are not directly comparable across independently processed source objects.",
  "- Current project_broad_celltype and neutrophil_confidence are retained only for review and are marked annotation_status=preliminary_unvalidated_task004.",
  "- Xue author labels remain available through source_author_annotation.",
  "- The known Task 004 rescued-status bug was corrected for CRA002308 and in_house in this review object only.",
  "",
  "## Output",
  "",
  paste0("Merged review object: ", out_path),
  paste0("SHA256: ", paste(sha, collapse = " ")),
  "",
  "## Important limitation",
  "",
  "This is a pure merge object, not an integrated atlas. Dataset-driven structure is expected if the merged counts are normalized/PCA/UMAPed without batch correction.",
  ""
)
writeLines(report, report_path)

message("DONE: ", out_path)
message("Cells: ", ncol(merged), "; features: ", nrow(merged))
