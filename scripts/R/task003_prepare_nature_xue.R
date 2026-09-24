#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) != 1L || hit == length(args)) stop("Missing argument: ", flag)
  args[[hit + 1L]]
}

source_path <- normalizePath(get_arg("--source"), mustWork = TRUE)
out_path <- get_arg("--out")
audit_path <- get_arg("--audit")
manifest_path <- get_arg("--manifest")

message("Loading author Seurat object read-only: ", source_path)
x <- readRDS(source_path)
if (!inherits(x, "Seurat")) stop("Author object is not a Seurat object: ", paste(class(x), collapse = ","))
original_meta <- x[[]]
original_cells <- ncol(x)
original_features <- nrow(x)
original_assays <- Assays(x)
original_layers <- setNames(lapply(original_assays, function(a) Layers(x[[a]])), original_assays)
original_reductions <- Reductions(x)
original_metadata_columns <- colnames(original_meta)

if (!"Cancer_type" %in% colnames(original_meta)) stop("Expected Cancer_type metadata is absent")
if (!"Sample" %in% colnames(original_meta)) stop("Expected Sample metadata is absent")
if (!"clusters" %in% colnames(original_meta)) stop("Expected clusters metadata is absent")

source_cancer_type <- as.character(original_meta$Cancer_type)
source_sample <- as.character(original_meta$Sample)
# AL is a tissue label and can also occur in adjacent samples from ICC/CHC or
# other non-HCC entities. Retain AL only when the source sample identifier
# explicitly identifies an HCC patient/sample; this avoids treating all AL as
# HCC adjacent liver.
keep <- source_cancer_type == "HCC" | (source_cancer_type == "AL" & grepl("(^|_)HCC(_|$)", source_sample))
if (!any(keep)) stop("No HCC/AL cells found")
keep_cells <- rownames(original_meta)[keep]
obj <- subset(x, cells = keep_cells)

# Preserve the original cell identity before applying globally unique project IDs.
source_cell_id <- colnames(obj)
source_meta <- obj[[]]
obj[["source_author_cell_id"]] <- source_cell_id
for (nm in original_metadata_columns) {
  obj[[paste0("source_author_", nm)]] <- source_meta[[nm]]
}
obj[["source_author_annotation"]] <- as.character(source_meta$clusters)

# Project metadata are separate from the author fields. Cancer_type HCC is used
# as tumor and AL as adjacent liver, exactly as encoded by the uploaded object.
obj[["dataset"]] <- "nature_xue"
project_sample_id <- as.character(source_meta$Sample)
project_tissue <- ifelse(as.character(source_meta$Cancer_type) == "AL", "Adjacent", "Tumor")
base_patient <- project_sample_id
base_patient <- sub("_HCC_N$", "", base_patient)
base_patient <- sub("_HCC_IM[0-9]*$", "", base_patient)
base_patient <- sub("_HCC$", "", base_patient)
if (any(!grepl("^A[0-9]+$", base_patient))) {
  stop("Unexpected HCC Sample naming after base-patient parsing: ",
       paste(sort(unique(project_sample_id[!grepl("^A[0-9]+$", base_patient)])), collapse = ","))
}
project_patient_id <- paste0("nature_xue_", base_patient)
sample_design <- unique(data.table(
  sample_id = project_sample_id,
  patient_id = project_patient_id,
  tissue = project_tissue
))
patient_design <- sample_design[, .(
  has_tumor = any(tissue == "Tumor"),
  has_adjacent = any(tissue == "Adjacent")
), by = patient_id]
patient_design[, paired_status := ifelse(has_tumor & has_adjacent, "paired", "unpaired")]
patient_design[, paired_id := ifelse(paired_status == "paired", patient_id, "unknown")]
pair_status_map <- setNames(patient_design$paired_status, patient_design$patient_id)
pair_id_map <- setNames(patient_design$paired_id, patient_design$patient_id)

obj[["sample_id"]] <- project_sample_id
obj[["patient_id"]] <- project_patient_id
obj[["tissue"]] <- project_tissue
obj[["paired_status"]] <- unname(pair_status_map[project_patient_id])
obj[["paired_id"]] <- unname(pair_id_map[project_patient_id])
obj[["etiology"]] <- "unknown"
obj[["MVI"]] <- "unknown"
obj[["platform"]] <- "author_processed_Seurat"
obj[["selection_strategy"]] <- "author_processed_object; uploaded metadata do not document whole-tissue versus enriched/sorted sampling"
obj[["matrix_type"]] <- "author_processed_Seurat_assay"
counts_available <- "counts" %in% Layers(obj[["RNA"]])
obj[["counts_available"]] <- counts_available
obj[["qc_provenance"]] <- "author_processed_object_no_task002_threshold_reapplication"
obj[["abundance_eligible"]] <- "CONDITIONAL"
obj[["analysis_inclusion"]] <- "YES"
obj[["qc_status"]] <- "AUTHOR_PROCESSED_REVIEW"

colnames(obj) <- make.unique(paste0("nature_xue__", source_cell_id))
obj@misc$task003 <- list(
  task = "task_003",
  source_file = source_path,
  source_object_preserved_unchanged = TRUE,
  selection_design = "not documented in uploaded object metadata; abundance eligibility is CONDITIONAL",
  subset_rule = "Cancer_type HCC plus Cancer_type AL only when Sample explicitly contains _HCC; AL harmonized to Adjacent; non-HCC adjacent samples excluded",
  author_annotation_field = "source_author_annotation (copied from source_author_clusters)",
  raw_counts_preserved = counts_available,
  counts_available = counts_available,
  original_dimensions = c(features = original_features, cells = original_cells),
  subset_dimensions = c(features = nrow(obj), cells = ncol(obj)),
  original_assays = original_assays,
  original_layers = original_layers,
  original_reductions = original_reductions,
  original_metadata_columns = original_metadata_columns
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, out_path, compress = TRUE)

md <- obj[[]]
sample_manifest <- as.data.table(md[, c("dataset", "sample_id", "patient_id", "tissue", "paired_status", "paired_id", "etiology", "MVI", "platform", "selection_strategy", "matrix_type", "counts_available", "qc_provenance", "qc_status", "abundance_eligible", "analysis_inclusion"), drop = FALSE])
sample_manifest[, n_cells := .N, by = .(dataset, sample_id, patient_id, tissue, paired_status, paired_id)]
sample_manifest <- unique(sample_manifest)
fwrite(sample_manifest, manifest_path)

audit <- data.table(
  dataset = "nature_xue",
  source_path = source_path,
  source_class = paste(class(x), collapse = ";"),
  seurat_version = as.character(packageVersion("Seurat")),
  original_n_features = original_features,
  original_n_cells = original_cells,
  subset_n_features = nrow(obj),
  subset_n_cells = ncol(obj),
  source_cancer_types_observed = paste(sort(unique(as.character(original_meta$Cancer_type))), collapse = ";"),
  retained_cancer_types = "HCC;AL (HCC-linked AL only)",
  excluded_cancer_types = paste(sort(setdiff(unique(as.character(original_meta$Cancer_type)), c("HCC", "AL"))), collapse = ";"),
  assays = paste(original_assays, collapse = ";"),
  RNA_layers = paste(original_layers[["RNA"]], collapse = ";"),
  reductions = paste(original_reductions, collapse = ";"),
  metadata_columns = paste(original_metadata_columns, collapse = ";"),
  source_author_annotation = "clusters copied to source_author_annotation and source_author_clusters",
  counts_available = counts_available,
  selection_strategy = "author-processed object; selection/sampling provenance not encoded in uploaded metadata",
  abundance_eligible = "CONDITIONAL",
  analysis_inclusion = "YES",
  original_object_unchanged = TRUE,
  global_cell_ids = !anyDuplicated(colnames(obj)),
  n_sample_rows = uniqueN(sample_manifest$sample_id),
  n_unique_patients = uniqueN(sample_manifest$patient_id),
  n_paired_patients = uniqueN(sample_manifest[paired_status == "paired", patient_id]),
  patient_pairing_rule = "strip _HCC, _HCC_N and _HCC_IM[0-9]* to base A-number; paired when patient has both Tumor and Adjacent"
)
fwrite(audit, audit_path)
message("Saved ", out_path, " with ", ncol(obj), " cells")
