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

project_root <- normalizePath(get_arg("--project-root"), mustWork = TRUE)
results_root <- normalizePath(get_arg("--results-root"), mustWork = TRUE)

atomic_save_rds <- function(obj, path) {
  tmp <- paste0(path, ".tmp_metadata_correction")
  saveRDS(obj, tmp, compress = TRUE)
  test <- readRDS(tmp)
  if (!inherits(test, "Seurat") || ncol(test) != ncol(obj)) stop("Validation failed for ", tmp)
  rm(test)
  if (!file.rename(tmp, path)) stop("Atomic replace failed for ", path)
}

nature_path <- file.path(project_root, "objects", "task003_extension", "nature_xue", "nature_xue_HCC_T_AL_prepared.rds")
if (!file.exists(nature_path)) stop("Missing prepared nature_xue object: ", nature_path)
message("Correcting nature_xue pairing metadata: ", nature_path)
obj <- readRDS(nature_path)
md <- obj[[]]
if (!all(c("sample_id", "tissue") %in% colnames(md))) stop("nature_xue prepared object lacks sample_id/tissue")

sample_id <- as.character(md$sample_id)
base_patient <- sample_id
base_patient <- sub("_HCC_N$", "", base_patient)
base_patient <- sub("_HCC_IM[0-9]*$", "", base_patient)
base_patient <- sub("_HCC$", "", base_patient)
if (any(!grepl("^A[0-9]+$", base_patient))) {
  bad <- sort(unique(sample_id[!grepl("^A[0-9]+$", base_patient)]))
  stop("Unexpected nature_xue HCC sample naming: ", paste(bad, collapse = ","))
}
patient_id <- paste0("nature_xue_", base_patient)

sample_dt <- unique(data.table(
  sample_id = sample_id,
  base_patient = base_patient,
  patient_id = patient_id,
  tissue = as.character(md$tissue)
))
patient_design <- sample_dt[, .(
  has_tumor = any(tissue == "Tumor"),
  has_adjacent = any(tissue == "Adjacent"),
  n_tumor_samples = sum(tissue == "Tumor"),
  n_adjacent_samples = sum(tissue == "Adjacent"),
  tumor_samples = paste(sort(sample_id[tissue == "Tumor"]), collapse = ";"),
  adjacent_samples = paste(sort(sample_id[tissue == "Adjacent"]), collapse = ";")
), by = .(base_patient, patient_id)]
patient_design[, paired_status := ifelse(has_tumor & has_adjacent, "paired", "unpaired")]
patient_design[, paired_id := ifelse(paired_status == "paired", patient_id, "unknown")]

pair_map <- setNames(patient_design$paired_status, patient_design$patient_id)
paired_id_map <- setNames(patient_design$paired_id, patient_design$patient_id)
obj[["patient_id"]] <- patient_id
obj[["paired_status"]] <- unname(pair_map[patient_id])
obj[["paired_id"]] <- unname(paired_id_map[patient_id])

obj@misc$task003_metadata_correction <- list(
  corrected = TRUE,
  correction_date = "2026-09-24",
  rule = "Sample base A-number defines patient: _HCC, _HCC_N and _HCC_IM[0-9]* map to the same patient",
  n_sample_rows = nrow(sample_dt),
  n_unique_patients = nrow(patient_design),
  n_paired_patients = sum(patient_design$paired_status == "paired"),
  multiregion_patients = patient_design[n_tumor_samples > 1, base_patient]
)
atomic_save_rds(obj, nature_path)

md2 <- obj[[]]
manifest <- as.data.table(md2[, c(
  "dataset", "sample_id", "patient_id", "tissue", "paired_status", "paired_id",
  "etiology", "MVI", "platform", "selection_strategy", "matrix_type",
  "counts_available", "qc_provenance", "qc_status", "abundance_eligible",
  "analysis_inclusion"
), drop = FALSE])
manifest[, n_cells := .N, by = .(
  dataset, sample_id, patient_id, tissue, paired_status, paired_id
)]
manifest <- unique(manifest)
setorder(manifest, sample_id)
fwrite(manifest, file.path(results_root, "task003_nature_xue_sample_manifest.csv"))

pairing_audit <- merge(
  patient_design,
  manifest[, .(
    tumor_cells = sum(n_cells[tissue == "Tumor"]),
    adjacent_cells = sum(n_cells[tissue == "Adjacent"])
  ), by = patient_id],
  by = "patient_id",
  all.x = TRUE
)
setorder(pairing_audit, base_patient)
fwrite(pairing_audit, file.path(results_root, "task003_nature_xue_pairing_audit.csv"))

audit_path <- file.path(results_root, "task003_nature_xue_subset_audit.csv")
audit <- fread(audit_path)
audit[, n_sample_rows := nrow(manifest)]
audit[, n_unique_patients := uniqueN(manifest$patient_id)]
audit[, n_paired_patients := uniqueN(manifest[paired_status == "paired", patient_id])]
audit[, patient_pairing_rule := "strip _HCC, _HCC_N, and _HCC_IM[0-9]* to base A-number; paired when patient has both Tumor and Adjacent"]
fwrite(audit, audit_path)

cra_dir <- file.path(project_root, "objects", "task003_extension", "CRA002308")
cra_files <- sort(list.files(cra_dir, pattern = "_prepared\\.rds$", full.names = TRUE))
if (length(cra_files) != 14L) stop("Expected 14 CRA002308 prepared RDS files, found ", length(cra_files))
for (path in cra_files) {
  x <- readRDS(path)
  x[["abundance_eligible"]] <- "CONDITIONAL"
  x[["selection_strategy"]] <- "flow-sorted live nucleated cells after doublet exclusion; no lineage-specific immune enrichment documented; paired within-cohort abundance sensitivity analysis allowed"
  x@misc$task003_metadata_correction <- list(
    corrected = TRUE,
    correction_date = "2026-09-24",
    abundance_rule = "CONDITIONAL: live nucleated-cell sorting can alter composition, but no lineage-specific immune enrichment is documented; use only paired within-cohort abundance sensitivity analysis"
  )
  atomic_save_rds(x, path)
}

message("Task 003 metadata correction complete: nature_xue patients/pairs corrected and CRA002308 abundance eligibility set to CONDITIONAL.")
