#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) != 1L || hit == length(args)) stop("Missing argument: ", flag)
  args[[hit + 1L]]
}

stage_manifest_path <- normalizePath(get_arg("--stage-manifest"), mustWork = TRUE)
out_root <- normalizePath(get_arg("--out-root"), mustWork = FALSE)
validation_path <- get_arg("--validation")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

stage <- fread(stage_manifest_path, sep = ",", na.strings = c("", "NA"), data.table = TRUE)
required <- c("dataset", "sample_id", "matrix_path", "features_path", "barcodes_path", "cell_qc_path")
if (!all(required %in% names(stage))) stop("Stage manifest is missing required columns")

read_features <- function(path) {
  x <- fread(path, sep = "\t", header = FALSE, colClasses = "character", data.table = FALSE)
  if (ncol(x) < 2L) stop("Feature file has fewer than two columns: ", path)
  list(id = x[[1L]], name = make.unique(x[[2L]]), type = if (ncol(x) >= 3L) x[[3L]] else rep("Gene Expression", nrow(x)))
}

validation <- list()
for (i in seq_len(nrow(stage))) {
  row <- stage[i]
  dataset <- row$dataset
  sample_id <- row$sample_id
  message("Creating Task 003 Seurat object: ", dataset, "/", sample_id)

  counts <- readMM(row$matrix_path)
  counts <- as(counts, "dgCMatrix")
  features <- read_features(row$features_path)
  barcodes <- readLines(row$barcodes_path, warn = FALSE)
  if (nrow(counts) != length(features$name) || ncol(counts) != length(barcodes)) stop("Staged dimensions do not match identifiers for ", sample_id)
  rownames(counts) <- features$name
  colnames(counts) <- make.unique(barcodes)

  cell_qc <- fread(cmd = paste("gzip -cd", shQuote(row$cell_qc_path)), sep = "\t", data.table = TRUE)
  required_qc <- c("barcode", "candidate_neutrophil", "passes_qc", "percent.mt", "nCount_RNA", "nFeature_RNA")
  if (!all(required_qc %in% names(cell_qc))) stop("Cell QC file is missing required columns for ", sample_id)
  qc_index <- match(barcodes, cell_qc$barcode)
  if (anyNA(qc_index)) stop("Cell QC identifiers do not match staged barcodes for ", sample_id)
  cell_qc <- cell_qc[qc_index]
  meta <- as.data.frame(cell_qc[, setdiff(names(cell_qc), c("barcode", "nCount_RNA", "nFeature_RNA")), with = FALSE])
  rownames(meta) <- colnames(counts)
  meta$task003_nCount_RNA <- as.numeric(cell_qc$nCount_RNA)
  meta$task003_nFeature_RNA <- as.numeric(cell_qc$nFeature_RNA)
  meta$percent.mt <- as.numeric(cell_qc$percent.mt)
  meta$passes_qc <- as.logical(cell_qc$passes_qc)
  meta$candidate_neutrophil <- as.logical(cell_qc$candidate_neutrophil)
  meta$doublet_risk_proxy <- as.logical(cell_qc$doublet_risk_proxy)

  obj <- CreateSeuratObject(counts = counts, project = paste(dataset, sample_id, sep = "_"), min.cells = 0, min.features = 0, meta.data = meta)
  obj <- subset(obj, subset = passes_qc)
  obj@misc$task003 <- list(
    task = "task_003",
    dataset = dataset,
    sample_id = sample_id,
    source_matrix = row$source_member,
    raw_counts_preserved = TRUE,
    source_representation = "uploaded 10x cell-called count matrix",
    filtering = "corrected Task 002 identity-independent all-source-cell q01/q01/q98 bounded thresholds; doublet proxy retained as review flag",
    abundance_eligible = row$abundance_eligible,
    n_cells_before_qc = nrow(cell_qc),
    n_cells_after_qc = ncol(obj)
  )

  out_dir <- file.path(out_root, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, paste0(sample_id, "_prepared.rds"))
  saveRDS(obj, out_path, compress = TRUE)
  count_layer <- GetAssayData(obj, assay = "RNA", layer = "counts")
  validation[[length(validation) + 1L]] <- data.table(
    dataset = dataset,
    sample_id = sample_id,
    rds_path = out_path,
    n_features = nrow(obj),
    n_cells_before_qc = nrow(cell_qc),
    n_cells_after_qc = ncol(obj),
    raw_count_layer_present = "counts" %in% Layers(obj[["RNA"]]),
    counts_nonnegative = all(count_layer@x >= 0),
    n_candidate_after_qc = sum(obj$candidate_neutrophil, na.rm = TRUE),
    globally_unique_cell_ids = !anyDuplicated(colnames(obj))
  )
  rm(obj, counts, cell_qc, meta)
  gc(verbose = FALSE)
}

fwrite(rbindlist(validation, fill = TRUE), validation_path)
message("Wrote ", length(validation), " Task 003 matrix-cohort Seurat objects")
