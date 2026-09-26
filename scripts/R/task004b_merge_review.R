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
qs_path <- arg_value(
  "qs-out",
  file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.qs")
)
h5ad_path <- arg_value(
  "h5ad-out",
  file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.h5ad")
)
cohort_root <- arg_value(
  "cohort-root",
  file.path(project_root, "objects", "task004_merge_review", "cohort_merged")
)
checkpoint_path <- arg_value(
  "checkpoint",
  file.path(project_root, "objects", "task004_merge_review", "checkpoint",
            "HCC_TA_8datasets_merged_review_v1_checkpoint.rds")
)
results_root <- arg_value("results-root", file.path(project_root, "results"))
report_path <- arg_value(
  "report",
  file.path(project_root, "reports", "task_004b_merge_review_report.md")
)

dir.create(dirname(qs_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(h5ad_path), recursive = TRUE, showWarnings = FALSE)
dir.create(cohort_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(checkpoint_path), recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

expected_datasets <- c(
  "GSE282701", "GSE242889", "GSE326201", "GSE149614", "GSE299340",
  "CRA002308", "nature_xue", "in_house"
)

resume_from_checkpoint <- file.exists(checkpoint_path)

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
    meta.data = md
  )
  if (!identical(colnames(slim), colnames(obj))) {
    stop("Cell IDs changed while slimming: ", path)
  }

  rm(obj, counts, md)
  gc(verbose = FALSE)
  slim
}

merge_two <- function(x, y) {
  merge(x = x, y = y, merge.data = FALSE, merge.dr = FALSE)
}

if (resume_from_checkpoint) {
  message("Existing merge checkpoint detected; skipping 103-object merge: ", checkpoint_path)
  merged <- readRDS(checkpoint_path)
  if (!inherits(merged, "Seurat")) stop("Checkpoint is not a Seurat object")
} else {
cohort_rows <- list()

for (dataset_name in expected_datasets) {
  acc <- NULL
  cohort_path <- file.path(cohort_root, paste0(dataset_name, "_merged_counts_metadata.rds"))
  expected_dataset_cells <- sum(as.numeric(manifest[dataset == dataset_name, annotated_object_cells]))

  if (file.exists(cohort_path)) {
    message("Checking existing cohort checkpoint: ", cohort_path)
    acc <- tryCatch(readRDS(cohort_path), error = function(e) {
      warning("Existing cohort checkpoint is unreadable and will be rebuilt: ", conditionMessage(e))
      NULL
    })

    if (!is.null(acc)) {
      existing_layers <- if (inherits(acc, "Seurat") && "RNA" %in% Assays(acc)) {
        Layers(acc[["RNA"]])
      } else {
        character()
      }
      existing_count_layers <- existing_layers[grepl("^counts", existing_layers)]
      if (!inherits(acc, "Seurat") || !"RNA" %in% Assays(acc) ||
          !length(existing_count_layers) ||
          ncol(acc) != expected_dataset_cells ||
          anyDuplicated(colnames(acc)) ||
          !setequal(unique(as.character(acc$dataset)), dataset_name)) {
        warning("Existing cohort checkpoint failed structural validation and will be rebuilt: ", cohort_path)
        rm(acc)
        acc <- NULL
      } else if (!identical(existing_layers, "counts")) {
        message("Stripping non-count RNA layers from existing cohort checkpoint")
        if (length(existing_count_layers) > 1L) {
          acc[["RNA"]] <- JoinLayers(
            acc[["RNA"]],
            layers = existing_count_layers,
            new = "counts"
          )
        }
        md <- acc[[]]
        counts <- LayerData(acc[["RNA"]], layer = "counts")
        slim <- CreateSeuratObject(
          counts = counts,
          assay = "RNA",
          project = dataset_name,
          meta.data = md
        )
        if (!identical(colnames(slim), colnames(acc))) {
          stop("Cell IDs changed while sanitizing cohort checkpoint: ", cohort_path)
        }
        saveRDS(slim, cohort_path, compress = "gzip")
        rm(acc, counts, md)
        acc <- slim
        rm(slim)
        gc(verbose = FALSE)
      }
    }
  }

  if (is.null(acc)) {
    files <- manifest[dataset == dataset_name, output_path]
    message("Preparing cohort ", dataset_name, " from ", length(files), " source object(s)")
    for (i in seq_along(files)) {
      message("  [", i, "/", length(files), "] ", basename(files[[i]]))
      slim <- prepare_for_merge(files[[i]], dataset_name)
      if (is.null(acc)) {
        acc <- slim
      } else {
        acc <- merge_two(acc, slim)
      }
      rm(slim)
      gc(verbose = FALSE)
    }
    # Keep the resumable intermediate uncompressed to avoid the large peak RAM
    # incurred by gzip serialization of the nature_xue cohort.
    saveRDS(acc, cohort_path, compress = FALSE)
  }

  cohort_rows[[length(cohort_rows) + 1L]] <- data.table(
    dataset = dataset_name,
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
for (dataset_name in expected_datasets) {
  path <- cohort_summary[dataset == dataset_name, cohort_object_path][[1]]
  message("Final merge: ", dataset_name)
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
} # end fresh merge path

if (ncol(merged) != expected_total) {
  stop("Merged/checkpoint object cell count mismatch: ", ncol(merged), " vs ", expected_total)
}
if (anyDuplicated(colnames(merged))) {
  stop("Merged/checkpoint object has duplicate cell IDs")
}
if (!setequal(unique(as.character(merged$dataset)), expected_datasets)) {
  stop("Merged/checkpoint object dataset metadata mismatch")
}

meta_preexport <- merged[[]]
if (data.table::uniqueN(meta_preexport$dataset) != 8L) {
  stop("Final merge validation failed: expected 8 datasets")
}
if (data.table::uniqueN(meta_preexport$sample_id) != 194L) {
  stop(
    "Final merge validation failed: expected 194 sample IDs, found ",
    data.table::uniqueN(meta_preexport$sample_id)
  )
}
if (data.table::uniqueN(meta_preexport$patient_id) != 132L) {
  stop(
    "Final merge validation failed: expected 132 patient IDs, found ",
    data.table::uniqueN(meta_preexport$patient_id)
  )
}
if (sum(meta_preexport$tissue == "Tumor", na.rm = TRUE) != 1039293L) {
  stop("Final merge validation failed: Tumor cell count mismatch")
}
if (sum(meta_preexport$tissue == "Adjacent", na.rm = TRUE) != 451559L) {
  stop("Final merge validation failed: Adjacent cell count mismatch")
}

rna_layers_final <- Layers(merged[["RNA"]])
if (!identical(rna_layers_final, "counts")) {
  count_layers_final <- rna_layers_final[grepl("^counts", rna_layers_final)]
  if (!length(count_layers_final) || length(setdiff(rna_layers_final, count_layers_final))) {
    stop(
      "Merged review object has unexpected RNA layers; found: ",
      paste(rna_layers_final, collapse = ", ")
    )
  }
  message("Joining ", length(count_layers_final), " merged RNA counts layers")
  merged[["RNA"]] <- JoinLayers(
    merged[["RNA"]],
    layers = count_layers_final,
    new = "counts"
  )
  rna_layers_final <- Layers(merged[["RNA"]])
  if (!identical(rna_layers_final, "counts")) {
    stop(
      "Merged review object must contain exactly one RNA counts layer after JoinLayers; found: ",
      paste(rna_layers_final, collapse = ", ")
    )
  }
}

merge_validation_status <- "VALIDATED"
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

if (!resume_from_checkpoint) {
  message("Writing recoverable merge checkpoint: ", checkpoint_path)
  # The checkpoint is server-side and temporary; avoid compression peaks for
  # the full eight-cohort object.
  saveRDS(merged, checkpoint_path, compress = FALSE)
}

qs_status <- "NOT_ATTEMPTED"
h5ad_status <- "NOT_ATTEMPTED"

if (requireNamespace("qs", quietly = TRUE)) {
  message("Writing Seurat QS object: ", qs_path)
  qs::qsave(merged, qs_path, preset = "high")
  message("Reload-validating QS object...")
  qs_check <- qs::qread(qs_path)
  if (!inherits(qs_check, "Seurat")) stop("QS validation failed: object is not Seurat")
  if (ncol(qs_check) != expected_total) stop("QS validation failed: cell count mismatch")
  if (!identical(colnames(qs_check), colnames(merged))) {
    stop("QS validation failed: ordered cell IDs changed")
  }
  qs_status <- "VALIDATED"
  rm(qs_check)
  gc(verbose = FALSE)
} else {
  qs_status <- "BLOCKED_MISSING_qs"
  warning("QS export skipped: R package 'qs' is not installed. Merge checkpoint retained.")
}

if (requireNamespace("anndataR", quietly = TRUE) &&
    requireNamespace("rhdf5", quietly = TRUE)) {
  message("Writing AnnData H5AD natively from Seurat with anndataR: ", h5ad_path)
  anndataR::write_h5ad(
    merged,
    path = h5ad_path,
    compression = "gzip",
    assay_name = "RNA",
    x_mapping = "counts",
    layers_mapping = FALSE,
    obs_mapping = TRUE,
    var_mapping = FALSE,
    obsm_mapping = FALSE,
    varm_mapping = FALSE,
    obsp_mapping = FALSE,
    varp_mapping = FALSE,
    uns_mapping = FALSE
  )

  message("Backed-validating H5AD with anndataR...")
  ad_check <- anndataR::read_h5ad(h5ad_path, as = "HDF5AnnData", mode = "r")
  ad_dims <- dim(ad_check)
  if (length(ad_dims) != 2L) stop("H5AD validation failed: could not determine dimensions")
  if (as.integer(ad_dims[[1L]]) != expected_total) {
    stop("H5AD validation failed: observation count mismatch")
  }
  if (as.integer(ad_dims[[2L]]) != nrow(merged)) {
    stop("H5AD validation failed: feature count mismatch")
  }
  h5ad_status <- "VALIDATED"
  rm(ad_check)
  gc(verbose = FALSE)
} else {
  miss <- c()
  if (!requireNamespace("anndataR", quietly = TRUE)) miss <- c(miss, "anndataR")
  if (!requireNamespace("rhdf5", quietly = TRUE)) miss <- c(miss, "rhdf5")
  h5ad_status <- paste0("BLOCKED_MISSING_", paste(miss, collapse = "+"))
  warning(
    "H5AD export skipped; missing R package(s): ", paste(miss, collapse = ", "),
    ". Merge checkpoint retained."
  )
}

task004b_status <- if (
  identical(merge_validation_status, "VALIDATED") &&
  identical(qs_status, "VALIDATED") &&
  identical(h5ad_status, "VALIDATED")
) {
  "COMPLETED"
} else if (identical(merge_validation_status, "VALIDATED")) {
  "MERGE_COMPLETED_EXPORT_PARTIAL"
} else {
  "FAILED"
}

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
  task004b_status = task004b_status,
  merge_validation_status = merge_validation_status,
  qs_object_path = qs_path,
  h5ad_object_path = h5ad_path,
  checkpoint_path = checkpoint_path,
  checkpoint_size_bytes = if (file.exists(checkpoint_path)) file.info(checkpoint_path)$size else NA_real_,
  qs_status = qs_status,
  h5ad_status = h5ad_status,
  qs_size_bytes = if (file.exists(qs_path)) file.info(qs_path)$size else NA_real_,
  h5ad_size_bytes = if (file.exists(h5ad_path)) file.info(h5ad_path)$size else NA_real_,
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

sha_qs <- if (file.exists(qs_path)) tryCatch(
  system2("sha256sum", qs_path, stdout = TRUE, stderr = TRUE),
  error = function(e) paste("sha256sum unavailable:", conditionMessage(e))
) else "NOT_CREATED"
sha_h5ad <- if (file.exists(h5ad_path)) tryCatch(
  system2("sha256sum", h5ad_path, stdout = TRUE, stderr = TRUE),
  error = function(e) paste("sha256sum unavailable:", conditionMessage(e))
) else "NOT_CREATED"

report <- c(
  "# Task 004b report — eight-cohort unintegrated Seurat merge for annotation review",
  "",
  "## Status",
  "",
  paste0("Task 004b status: ", task004b_status),
  paste0("Final merge validation: ", merge_validation_status),
  paste0("QS export: ", qs_status),
  paste0("H5AD export: ", h5ad_status),
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
  "- Samples: 194",
  "- Patients: 132",
  "- Explicit paired Tumor-Adjacent patients: 59 (metadata design reference)",
  "- Tumor cells: 1,039,293",
  "- Adjacent cells: 451,559",
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
  paste0("Seurat QS object: ", qs_path),
  paste0("AnnData H5AD object: ", h5ad_path),
  paste0("QS status: ", qs_status),
  paste0("H5AD status: ", h5ad_status),
  paste0("QS SHA256: ", paste(sha_qs, collapse = " ")),
  paste0("H5AD SHA256: ", paste(sha_h5ad, collapse = " ")),
  paste0("Recoverable merge checkpoint: ", checkpoint_path),
  "",
  "The H5AD is written natively from the Seurat object using anndataR, with RNA counts in AnnData X and cell metadata in obs.",
  "",
  "## Important limitation",
  "",
  "This is a pure merge object, not an integrated atlas. Dataset-driven structure is expected if the merged counts are normalized/PCA/UMAPed without batch correction.",
  ""
)
writeLines(report, report_path)

if (identical(task004b_status, "COMPLETED")) {
  message("Task 004b COMPLETED: both requested exports validated; removing temporary checkpoint.")
  unlink(checkpoint_path)
  message("DONE QS: ", qs_path)
  message("DONE H5AD: ", h5ad_path)
} else {
  message("Task 004b PARTIAL: final merge validated but one or more format exports are blocked.")
  message("Checkpoint retained for export-only resume: ", checkpoint_path)
}
message("Cells: ", ncol(merged), "; features: ", nrow(merged))
