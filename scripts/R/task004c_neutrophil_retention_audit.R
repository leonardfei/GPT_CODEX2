#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas"
qs_path <- file.path(root, "objects", "HCC_TA_8datasets_merged_review_v1.qs")
out_path <- file.path(root, "results", "task004c_neutrophil_retention_audit.csv")

dim.task004b_chunked_counts <- function(x) c(length(x$feature_names), length(x$cell_names))

obj <- qs::qread(qs_path, use_alt_rep = FALSE, nthreads = 8L)
if (ncol(obj) != 1490852L || nrow(obj) != 38025L) stop("Final QS dimensions mismatch")

md <- as.data.table(obj[[]], keep.rownames = "cell_id")
required <- c(
  "dataset", "project_sample_id", "project_patient_id", "tissue",
  "project_broad_celltype", "source_nCount_RNA_before_feature_harmonisation",
  "source_nFeature_RNA_before_feature_harmonisation", "hgnc_nCount_RNA",
  "hgnc_nFeature_RNA", "percent.mt"
)
missing <- setdiff(required, names(md))
if (length(missing)) stop("Final QS lacks audit fields: ", paste(missing, collapse = ", "))
if (anyDuplicated(md$cell_id) || !identical(md$cell_id, colnames(obj))) {
  stop("Final QS cell IDs are duplicated or out of order")
}

# This is the preliminary Task 004 broad label, retained for review; it is not a
# new cell-identity call. Task 004c performs no cell filtering.
md[, qc_group := fifelse(
  !is.na(project_broad_celltype) & project_broad_celltype == "neutrophil",
  "preliminary_neutrophil",
  "other_cells"
)]
md[, candidate_neutrophil := qc_group == "preliminary_neutrophil"]
n_candidates <- sum(md$candidate_neutrophil)
if (n_candidates != 65546L) stop("Preliminary neutrophil count differs from Task 004 hold-point report")

sample_keys <- c("dataset", "project_sample_id", "project_patient_id", "tissue")
retention <- md[, .(
  candidate_neutrophils_before = sum(candidate_neutrophil),
  candidate_neutrophils_after = sum(candidate_neutrophil),
  cells_before = .N,
  cells_after = .N,
  removed_cells = 0L,
  removal_reason = "none; Task 004c did not filter cells"
), by = sample_keys]
retention[, retention_fraction := fifelse(
  candidate_neutrophils_before > 0L,
  candidate_neutrophils_after / candidate_neutrophils_before,
  NA_real_
)]

distribution <- md[, .(
  n_cells = .N,
  source_nCount_RNA_q25 = as.numeric(quantile(source_nCount_RNA_before_feature_harmonisation, 0.25)),
  source_nCount_RNA_median = as.numeric(median(source_nCount_RNA_before_feature_harmonisation)),
  source_nCount_RNA_q75 = as.numeric(quantile(source_nCount_RNA_before_feature_harmonisation, 0.75)),
  hgnc_nCount_RNA_q25 = as.numeric(quantile(hgnc_nCount_RNA, 0.25)),
  hgnc_nCount_RNA_median = as.numeric(median(hgnc_nCount_RNA)),
  hgnc_nCount_RNA_q75 = as.numeric(quantile(hgnc_nCount_RNA, 0.75)),
  source_nFeature_RNA_q25 = as.numeric(quantile(source_nFeature_RNA_before_feature_harmonisation, 0.25)),
  source_nFeature_RNA_median = as.numeric(median(source_nFeature_RNA_before_feature_harmonisation)),
  source_nFeature_RNA_q75 = as.numeric(quantile(source_nFeature_RNA_before_feature_harmonisation, 0.75)),
  hgnc_nFeature_RNA_q25 = as.numeric(quantile(hgnc_nFeature_RNA, 0.25)),
  hgnc_nFeature_RNA_median = as.numeric(median(hgnc_nFeature_RNA)),
  hgnc_nFeature_RNA_q75 = as.numeric(quantile(hgnc_nFeature_RNA, 0.75)),
  percent_mt_q25 = as.numeric(quantile(percent.mt, 0.25)),
  percent_mt_median = as.numeric(median(percent.mt)),
  percent_mt_q75 = as.numeric(quantile(percent.mt, 0.75))
), by = c(sample_keys, "qc_group")]

audit <- merge(distribution, retention, by = sample_keys, sort = FALSE)
setcolorder(audit, c(
  sample_keys, "qc_group", "candidate_neutrophils_before",
  "candidate_neutrophils_after", "retention_fraction", "cells_before",
  "cells_after", "removed_cells", "removal_reason", "n_cells",
  setdiff(names(audit), c(sample_keys, "qc_group", "candidate_neutrophils_before",
                          "candidate_neutrophils_after", "retention_fraction",
                          "cells_before", "cells_after", "removed_cells",
                          "removal_reason", "n_cells"))
))
fwrite(audit, out_path)

cat(
  "Task 004c neutrophil-retention audit written:", out_path,
  "\ncandidate neutrophils before/after:", n_candidates, "/", n_candidates,
  "\nretention fraction: 1.0; cells removed: 0\n"
)
