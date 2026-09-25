#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  key <- paste0("--", name)
  hit <- which(args == key)
  if (length(hit) == 0L || hit == length(args)) return(default)
  args[[hit + 1L]]
}

project_root <- normalizePath(arg_value("project-root", "."), mustWork = TRUE)
out_root <- arg_value("out-root", file.path(project_root, "objects", "task004_annotated"))
results_root <- arg_value("results-root", file.path(project_root, "results"))
figures_root <- arg_value("figures-root", file.path(project_root, "figures"))
report_path <- arg_value("report", file.path(project_root, "reports", "task_004_report.md"))

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

known_datasets <- c("GSE282701", "GSE242889", "GSE326201", "GSE149614", "GSE299340",
                    "CRA002308", "nature_xue", "in_house")
abundance_role <- c(
  GSE282701 = "YES", GSE242889 = "YES", GSE326201 = "YES", GSE149614 = "YES",
  GSE299340 = "YES", CRA002308 = "CONDITIONAL", nature_xue = "CONDITIONAL",
  in_house = "CONDITIONAL"
)

programs <- list(
  `hepatocyte/tumor epithelial` = c("ALB", "APOA1", "APOA2", "TTR", "KRT8", "KRT18", "KRT19", "EPCAM", "KRT7"),
  `T/NK` = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB", "NKG7", "GNLY", "KLRD1", "CCL5"),
  B = c("CD79A", "CD79B", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD22", "CD83"),
  plasma = c("MZB1", "JCHAIN", "XBP1", "SEC11C", "DERL3", "IGHG1", "IGKC"),
  `monocyte/macrophage` = c("LST1", "FCER1G", "TYROBP", "CTSS", "LILRB1", "CTSD", "C1QC", "APOC1", "SPP1", "SAT1"),
  neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "S100A8", "S100A9", "FCAR", "ELANE", "MPO", "NAMPT", "MCEMP1", "FFAR2", "CXCR4"),
  dendritic = c("FCER1A", "CD1C", "CLEC10A", "CST3", "GZMB", "IRF7", "TCF4", "SERPINF1"),
  mast = c("TPSAB1", "TPSB2", "KIT", "MS4A2", "HDC", "CPA3"),
  endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "RAMP2", "ENG", "ESM1"),
  `fibroblast/mesenchymal` = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1", "COL6A2")
)
program_names <- names(programs)
neut_core <- c("FCGR3B", "CSF3R", "FPR1", "ELANE", "MPO", "FCAR")

as_character <- function(x, fallback = "unknown") {
  y <- as.character(x)
  y[is.na(y) | !nzchar(y)] <- fallback
  y
}

as_logical_safe <- function(x, default = FALSE) {
  if (is.null(x)) return(rep(default, 0L))
  if (is.logical(x)) return(!is.na(x) & x)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y", "pass", "passed")
}

get_meta_column <- function(md, candidates, fallback) {
  for (nm in candidates) if (nm %in% colnames(md)) return(as_character(md[[nm]], fallback))
  rep(fallback, nrow(md))
}

get_dataset <- function(path) {
  hits <- known_datasets[vapply(known_datasets, function(x) grepl(x, path, fixed = TRUE), logical(1))]
  if (length(hits) == 0L) stop("Cannot identify dataset from source path: ", path)
  hits[[1L]]
}

source_files <- sort(unique(c(
  list.files(file.path(project_root, "objects", "task002_corrected_seurat"),
             pattern = "\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
  list.files(file.path(project_root, "objects", "task003_extension"),
             pattern = "\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
)))
if (length(source_files) != 103L) {
  stop("Expected 103 corrected Task 002/003 RDS objects, found ", length(source_files))
}
message("Task 004 source objects: ", length(source_files))

score_programs <- function(expr, library_size) {
  score_mat <- matrix(0, nrow = ncol(expr), ncol = length(programs),
                      dimnames = list(colnames(expr), program_names))
  hit_mat <- matrix(0, nrow = ncol(expr), ncol = length(programs),
                    dimnames = list(colnames(expr), program_names))
  feature_upper <- toupper(rownames(expr))
  for (j in seq_along(programs)) {
    marker_idx <- match(toupper(programs[[j]]), feature_upper)
    marker_idx <- unique(marker_idx[!is.na(marker_idx)])
    if (length(marker_idx) == 0L) next
    submat <- expr[marker_idx, , drop = FALSE]
    umi <- Matrix::colSums(submat)
    hits <- Matrix::colSums(submat > 0)
    denom <- max(length(marker_idx), 1L)
    score_mat[, j] <- log1p((umi / pmax(library_size, 1)) * 10000 / denom) + hits / denom
    hit_mat[, j] <- hits
  }
  list(score = score_mat, hits = hit_mat)
}

extract_expression <- function(obj) {
  assay_name <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
  assay <- obj[[assay_name]]
  layer_names <- tryCatch(Layers(assay), error = function(e) character())
  layer_name <- if ("counts" %in% layer_names) "counts" else if ("data" %in% layer_names) "data" else layer_names[[1L]]
  if (is.null(layer_name) || !nzchar(layer_name)) stop("No usable RNA layer found")
  expr <- tryCatch(
    GetAssayData(obj, assay = assay_name, layer = layer_name),
    error = function(e) GetAssayData(obj, assay = assay_name, slot = layer_name)
  )
  list(expr = expr, assay = assay_name, layer = layer_name)
}

atomic_save_rds <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(obj, tmp, compress = "xz")
  if (!file.rename(tmp, path)) stop("Failed to move annotated object into place: ", path)
}

celltype_rows <- list()
sample_rows <- list()
cluster_rows <- list()
crosswalk_rows <- list()
dot_rows <- list()
embedding_rows <- list()
object_rows <- list()

for (file_index in seq_along(source_files)) {
  source_path <- source_files[[file_index]]
  dataset <- get_dataset(source_path)
  message(sprintf("[%d/%d] Reading %s", file_index, length(source_files), source_path))
  obj <- readRDS(source_path)
  md <- obj[[]]
  n_cells <- nrow(md)
  if (n_cells != ncol(obj)) stop("Metadata/object cell mismatch: ", source_path)

  sample_id <- get_meta_column(md, c("sample_id", "sample", "orig.ident"),
                               tools::file_path_sans_ext(basename(source_path)))
  patient_id <- get_meta_column(md, c("patient_id", "patient", "donor_id"), sample_id)
  tissue <- get_meta_column(md, c("tissue", "tissue_group", "source_tissue"), "unknown")
  role <- get_meta_column(md, c("abundance_eligible"), abundance_role[[dataset]])
  role[!role %in% c("YES", "CONDITIONAL", "NO")] <- abundance_role[[dataset]]
  author_annotation <- if ("source_author_annotation" %in% colnames(md)) {
    as_character(md$source_author_annotation, "(missing)")
  } else {
    rep("(missing)", n_cells)
  }

  passes_qc <- if ("passes_qc" %in% colnames(md)) as_logical_safe(md$passes_qc) else rep(TRUE, n_cells)
  passes_initial <- if ("passes_initial_qc" %in% colnames(md)) as_logical_safe(md$passes_initial_qc) else rep(TRUE, n_cells)
  rescued <- passes_qc & !passes_initial
  if (dataset == "nature_xue") {
    qc_status <- rep("author_processed", n_cells)
    rescued[] <- FALSE
  } else {
    qc_status <- ifelse(rescued, "corrected_qc_rescued", "corrected_qc_pass")
  }

  expr_info <- extract_expression(obj)
  expr <- expr_info$expr
  library_size <- Matrix::colSums(expr)
  scored <- score_programs(expr, library_size)
  scores <- scored$score
  hits <- scored$hits
  neut_i <- match("neutrophil", program_names)
  mono_i <- match("monocyte/macrophage", program_names)
  neut_hits <- hits[, neut_i]
  mono_hits <- hits[, mono_i]
  core_idx <- match(toupper(neut_core), toupper(rownames(expr)))
  core_idx <- unique(core_idx[!is.na(core_idx)])
  neut_core_hits <- if (length(core_idx)) Matrix::colSums(expr[core_idx, , drop = FALSE] > 0) else rep(0, n_cells)
  margin <- scores[, neut_i] - scores[, mono_i]

  confidence <- rep("none", n_cells)
  confidence[neut_hits >= 1] <- "low"
  confidence[neut_hits >= 3 & margin >= 0] <- "medium"
  confidence[neut_core_hits >= 1 & neut_hits >= 2 & margin > 0.05] <- "medium"
  confidence[neut_core_hits >= 1 & neut_hits >= 4 & margin > 0.10] <- "high"
  ambiguous <- neut_hits >= 2 & mono_hits >= 2 & abs(margin) < 0.10
  confidence[ambiguous] <- "ambiguous_myeloid"

  non_neut <- setdiff(seq_along(program_names), neut_i)
  non_neut_score <- scores[, non_neut, drop = FALSE]
  non_neut_hits <- hits[, non_neut, drop = FALSE]
  best_non_neut <- apply(non_neut_score, 1L, which.max)
  best_non_neut_score <- non_neut_score[cbind(seq_len(n_cells), best_non_neut)]
  best_non_neut_hits <- non_neut_hits[cbind(seq_len(n_cells), best_non_neut)]
  broad <- rep("other/uncertain", n_cells)
  strong_non_neut <- best_non_neut_hits >= 2 & best_non_neut_score > 0.15
  broad[strong_non_neut] <- program_names[non_neut][best_non_neut[strong_non_neut]]
  confident_neut <- confidence %in% c("high", "medium")
  broad[confident_neut] <- "neutrophil"
  broad[ambiguous] <- ifelse(scores[ambiguous, mono_i] >= scores[ambiguous, neut_i],
                             "monocyte/macrophage", "other/uncertain")
  broad[broad == "neutrophil" & !confident_neut] <- "other/uncertain"

  assigned_program <- match(broad, program_names)
  project_score <- rep(0, n_cells)
  project_hits <- rep(0, n_cells)
  valid_assigned <- !is.na(assigned_program)
  project_score[valid_assigned] <- scores[cbind(which(valid_assigned), assigned_program[valid_assigned])]
  project_hits[valid_assigned] <- hits[cbind(which(valid_assigned), assigned_program[valid_assigned])]

  cluster_id <- if ("seurat_clusters" %in% colnames(md)) {
    paste0("seurat_cluster_", as_character(md$seurat_clusters, "unknown"))
  } else if (dataset == "nature_xue" && any(author_annotation != "(missing)")) {
    paste0("author_cluster_", author_annotation)
  } else {
    paste0("sample_unit_", dataset, "__", sample_id)
  }
  cluster_type <- if (grepl("^seurat_cluster_", cluster_id[[1L]])) "seurat_cluster" else if (grepl("^author_cluster_", cluster_id[[1L]])) "author_annotation_unit" else "sample_review_unit"

  obj_meta <- data.frame(
    dataset = dataset,
    sample_id = sample_id,
    patient_id = patient_id,
    tissue = tissue,
    abundance_eligible = role,
    source_author_annotation = author_annotation,
    qc_status = qc_status,
    rescued_by_corrected_qc = rescued,
    project_broad_celltype = broad,
    neutrophil_confidence = confidence,
    neutrophil_marker_hits = as.integer(neut_hits),
    neutrophil_core_hits = as.integer(neut_core_hits),
    monocyte_marker_hits = as.integer(mono_hits),
    neutrophil_score_margin = as.numeric(margin),
    project_marker_score = as.numeric(project_score),
    project_marker_hits = as.integer(project_hits),
    task004_cluster_id = cluster_id,
    task004_cluster_type = cluster_type,
    row.names = rownames(md)
  )
  obj <- AddMetaData(obj, metadata = obj_meta)

  output_name <- sub("\\.rds$", "_annotated.rds", basename(source_path), ignore.case = TRUE)
  output_path <- file.path(out_root, dataset, output_name)
  atomic_save_rds(obj, output_path)

  cell_dt <- data.table(
    dataset = dataset, sample_id = sample_id, patient_id = patient_id, tissue = tissue,
    abundance_eligible = role, source_author_annotation = author_annotation,
    qc_status = qc_status, rescued_by_corrected_qc = rescued,
    project_broad_celltype = broad, neutrophil_confidence = confidence,
    task004_cluster_id = cluster_id, task004_cluster_type = cluster_type
  )

  ct <- cell_dt[, .(n_cells = .N, n_rescued = sum(rescued_by_corrected_qc),
                    n_samples = uniqueN(sample_id), n_patients = uniqueN(patient_id)),
               by = .(dataset, tissue, abundance_eligible, project_broad_celltype)]
  ct[, dataset_tissue_total := sum(n_cells), by = .(dataset, tissue)]
  ct[, fraction_within_dataset_tissue := n_cells / pmax(dataset_tissue_total, 1)]
  celltype_rows[[length(celltype_rows) + 1L]] <- ct

  sample_stat <- cell_dt[, .(
    n_cells = .N,
    n_neutrophil = sum(project_broad_celltype == "neutrophil"),
    n_neutrophil_high = sum(neutrophil_confidence == "high"),
    n_neutrophil_medium = sum(neutrophil_confidence == "medium"),
    n_ambiguous_myeloid = sum(neutrophil_confidence == "ambiguous_myeloid"),
    n_rescued = sum(rescued_by_corrected_qc),
    n_rescued_neutrophil = sum(rescued_by_corrected_qc & project_broad_celltype == "neutrophil"),
    n_patients = uniqueN(patient_id)
  ), by = .(dataset, sample_id, patient_id, tissue, abundance_eligible)]
  sample_stat[, summary_level := "sample"]
  sample_stat[, neutrophil_fraction := n_neutrophil / pmax(n_cells, 1)]
  sample_stat[, rescued_fraction := n_rescued / pmax(n_cells, 1)]
  patient_stat <- cell_dt[, .(
    sample_id = paste(sort(unique(sample_id)), collapse = ";"),
    n_cells = .N,
    n_neutrophil = sum(project_broad_celltype == "neutrophil"),
    n_neutrophil_high = sum(neutrophil_confidence == "high"),
    n_neutrophil_medium = sum(neutrophil_confidence == "medium"),
    n_ambiguous_myeloid = sum(neutrophil_confidence == "ambiguous_myeloid"),
    n_rescued = sum(rescued_by_corrected_qc),
    n_rescued_neutrophil = sum(rescued_by_corrected_qc & project_broad_celltype == "neutrophil"),
    n_patients = 1L
  ), by = .(dataset, patient_id, tissue, abundance_eligible)]
  patient_stat[, summary_level := "patient_tissue"]
  patient_stat[, neutrophil_fraction := n_neutrophil / pmax(n_cells, 1)]
  patient_stat[, rescued_fraction := n_rescued / pmax(n_cells, 1)]
  sample_rows[[length(sample_rows) + 1L]] <- rbindlist(list(sample_stat, patient_stat), fill = TRUE, use.names = TRUE)

  cluster_stat <- cell_dt[, .(
    n_cells = .N,
    n_rescued = sum(rescued_by_corrected_qc),
    rescued_fraction = sum(rescued_by_corrected_qc) / .N,
    n_neutrophil = sum(project_broad_celltype == "neutrophil"),
    n_ambiguous_myeloid = sum(neutrophil_confidence == "ambiguous_myeloid"),
    top_project_broad_celltype = names(sort(table(project_broad_celltype), decreasing = TRUE))[1L]
  ), by = .(dataset, sample_id, patient_id, tissue, abundance_eligible, task004_cluster_type, task004_cluster_id)]
  cluster_stat[, cluster_review_flag := ifelse(rescued_fraction > 0.5, "REVIEW_RESCUED_DOMINATED", "none")]
  cluster_stat[, review_note := ifelse(cluster_review_flag == "REVIEW_RESCUED_DOMINATED",
                                       "Retain; inspect marker coherence and QC rather than deleting rescued cells.",
                                       "No rescued-cell dominance flag.")]
  cluster_rows[[length(cluster_rows) + 1L]] <- cluster_stat

  cw <- cell_dt[, .(n_cells = .N), by = .(dataset, source_author_annotation, project_broad_celltype)]
  cw[, fraction_within_source_author_annotation := n_cells / sum(n_cells), by = .(dataset, source_author_annotation)]
  crosswalk_rows[[length(crosswalk_rows) + 1L]] <- cw

  take <- if (n_cells > 500L) sample.int(n_cells, 500L) else seq_len(n_cells)
  dot_dt <- as.data.table(scores[take, , drop = FALSE])
  dot_dt[, project_broad_celltype := broad[take]]
  dot_dt[, dataset := dataset]
  dot_long <- melt(dot_dt, id.vars = c("dataset", "project_broad_celltype"),
                   variable.name = "program", value.name = "score")
  dot_rows[[length(dot_rows) + 1L]] <- dot_long[, .(
    mean_score = mean(score), detection_fraction = mean(score > 0)
  ), by = .(dataset, project_broad_celltype, program)]
  embedding_rows[[length(embedding_rows) + 1L]] <- cbind(
    data.table(dataset = dataset, tissue = tissue[take], project_broad_celltype = broad[take]),
    as.data.table(scores[take, , drop = FALSE])
  )
  object_rows[[length(object_rows) + 1L]] <- data.table(
    dataset = dataset, source_path = source_path, output_path = output_path,
    source_object_cells = n_cells, annotated_object_cells = ncol(obj),
    expression_assay = expr_info$assay, expression_layer = expr_info$layer,
    n_rescued = sum(rescued), n_neutrophil = sum(broad == "neutrophil"),
    n_source_author_annotation = sum(author_annotation != "(missing)"),
    status = "saved"
  )
  rm(obj, expr, scores, hits, cell_dt, dot_dt, dot_long)
  gc(verbose = FALSE)
}

celltype_all <- rbindlist(celltype_rows, fill = TRUE, use.names = TRUE)
sample_all <- rbindlist(sample_rows, fill = TRUE, use.names = TRUE)
cluster_all <- rbindlist(cluster_rows, fill = TRUE, use.names = TRUE)
crosswalk_all <- rbindlist(crosswalk_rows, fill = TRUE, use.names = TRUE)
dot_all <- rbindlist(dot_rows, fill = TRUE, use.names = TRUE)
embedding_all <- rbindlist(embedding_rows, fill = TRUE, use.names = TRUE)
object_all <- rbindlist(object_rows, fill = TRUE, use.names = TRUE)

fwrite(celltype_all, file.path(results_root, "task004_celltype_counts.csv"))
fwrite(sample_all, file.path(results_root, "task004_neutrophil_by_sample.csv"))
fwrite(cluster_all, file.path(results_root, "task004_cluster_rescued_cell_audit.csv"))
fwrite(crosswalk_all, file.path(results_root, "task004_author_vs_project_annotation_crosswalk.csv"))
fwrite(dot_all, file.path(results_root, "task004_marker_program_summary.csv"))
fwrite(object_all, file.path(results_root, "task004_object_validation.csv"))

set.seed(4004)
plot_n <- min(nrow(embedding_all), 20000L)
plot_rows <- sample(seq_len(nrow(embedding_all)), plot_n)
plot_dt <- embedding_all[plot_rows]
plot_matrix <- as.matrix(plot_dt[, ..program_names])
plot_matrix[!is.finite(plot_matrix)] <- 0
embedding_method <- "PCA fallback"
if (nrow(plot_matrix) >= 3L && requireNamespace("uwot", quietly = TRUE)) {
  coords <- uwot::umap(plot_matrix, n_neighbors = min(30L, nrow(plot_matrix) - 1L),
                       min_dist = 0.3, n_threads = 1, verbose = FALSE)
  embedding_method <- "uwot UMAP in marker-program space"
} else {
  scaled <- scale(plot_matrix)
  scaled[!is.finite(scaled)] <- 0
  coords <- stats::prcomp(scaled, rank. = 2, center = FALSE, scale. = FALSE)$x[, 1:2, drop = FALSE]
}
plot_dt[, embedding_1 := coords[, 1L]]
plot_dt[, embedding_2 := coords[, 2L]]
pdf(file.path(figures_root, "task004_marker_program_umap.pdf"), width = 11, height = 8)
print(ggplot(plot_dt, aes(embedding_1, embedding_2, color = project_broad_celltype)) +
        geom_point(size = 0.25, alpha = 0.45) +
        labs(title = paste("Task 004 broad annotation review -", embedding_method),
             subtitle = "Source-wise marker-program embedding; no cross-cohort integration") +
        theme_classic() + theme(legend.position = "right"))
dev.off()

pdf(file.path(figures_root, "task004_marker_program_dotplot.pdf"), width = 15, height = 8)
print(ggplot(dot_all, aes(program, project_broad_celltype)) +
        geom_point(aes(size = detection_fraction, color = mean_score)) +
        facet_wrap(~ dataset, scales = "free_y", ncol = 2) +
        scale_x_discrete(drop = FALSE) +
        labs(title = "Task 004 marker-program dot plot",
             subtitle = "Sampled cells per source object; scores use coherent multi-marker programs",
             x = "Marker program", y = "Project broad cell type") +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

dataset_qc <- celltype_all[, .(
  n_cells = sum(n_cells), n_rescued = sum(n_rescued),
  rescued_fraction = sum(n_rescued) / sum(n_cells),
  n_neutrophil = sum(n_cells[project_broad_celltype == "neutrophil"])
), by = dataset]
dataset_qc[, neutrophil_fraction := n_neutrophil / pmax(n_cells, 1)]
cluster_qc <- cluster_all[, .(n_rescued_dominated = sum(cluster_review_flag == "REVIEW_RESCUED_DOMINATED")), by = dataset]
dataset_qc <- merge(dataset_qc, cluster_qc, by = "dataset", all.x = TRUE)
dataset_qc[is.na(n_rescued_dominated), n_rescued_dominated := 0L]
qc_long <- rbind(
  dataset_qc[, .(dataset, metric = "neutrophil_fraction", value = neutrophil_fraction)],
  dataset_qc[, .(dataset, metric = "rescued_fraction", value = rescued_fraction)]
)
pdf(file.path(figures_root, "task004_qc_annotation_review.pdf"), width = 12, height = 8)
print(ggplot(qc_long, aes(dataset, value, fill = metric)) +
        geom_col(position = "dodge") +
        labs(title = "Task 004 QC and annotation review", subtitle = "Rescued cells are retained; flagged review units are not deleted") +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

total_cells <- sum(object_all$annotated_object_cells)
total_neut <- sum(object_all$n_neutrophil)
total_rescued <- sum(object_all$n_rescued)
flagged <- cluster_all[cluster_review_flag == "REVIEW_RESCUED_DOMINATED"]
dataset_lines <- paste(sprintf("- `%s`: %s cells; %s neutrophil-classified; abundance role `%s`.",
                              dataset_qc$dataset, format(dataset_qc$n_cells, big.mark = ","),
                              format(dataset_qc$n_neutrophil, big.mark = ","),
                              abundance_role[dataset_qc$dataset]), collapse = "\n")
report_lines <- c(
  "# Task 004 report - harmonized broad annotation and neutrophil confirmation",
  "",
  "## Execution status",
  "",
  "Task 004 completed on the corrected Task 002 objects and corrected Task 003 extension objects. The source objects were not overwritten. Annotated copies were written server-side under `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task004_annotated/`. Task 005 integration was not executed.",
  "",
  sprintf("Processed %d source objects and %s annotated cells. The project broad layer classified %s cells as neutrophil; %s cells were marked `rescued_by_corrected_qc`.",
          nrow(object_all), format(total_cells, big.mark = ","), format(total_neut, big.mark = ","), format(total_rescued, big.mark = ",")),
  "",
  "## Method and provenance",
  "",
  "Each object was read and annotated independently using its RNA counts layer when available, otherwise its RNA data layer. Coherent multi-marker programs were scored for the ten requested lineage classes. Neutrophil confidence required coherent neutrophil marker hits and a score margin against the monocyte/macrophage program; ambiguous myeloid cells were not silently promoted to neutrophils.",
  "",
  "`source_author_annotation` was preserved. For `nature_xue`, `qc_status` is `author_processed`; source author annotations remain in the annotated objects and are summarized in the crosswalk output. For corrected count cohorts, cells rescued by corrected QC are explicitly labeled and retained.",
  "",
  "Objects without a pre-existing `seurat_clusters` field were audited at the sample review-unit level; `nature_xue` objects with source author labels were audited at the author-annotation-unit level. This is explicitly recorded in `task004_cluster_type` and does not imply a new clustering analysis.",
  "",
  "## Dataset summary",
  "",
  dataset_lines,
  "",
  "## Rescued-cell review",
  "",
  sprintf("There are %d rescued-dominated review units (rescued fraction > 0.5). None were deleted solely because of rescue status. The cluster audit CSV records the review flag, rescued fraction, neutrophil count, and top broad type.", nrow(flagged)),
  "",
  "## Abundance rules",
  "",
  "The five original cohorts retain abundance role `YES` for primary analyses subject to their study design. `CRA002308`, `nature_xue`, and `in_house` retain `CONDITIONAL` roles and must not be pooled as unbiased whole-tissue absolute fractions. Patient-level summaries use `patient_id`; multiple `nature_xue` regions remain sample-level records but are not treated as independent patient replicates.",
  "",
  "## Required outputs",
  "",
  "- `results/task004_celltype_counts.csv`",
  "- `results/task004_neutrophil_by_sample.csv` (sample and patient-tissue summaries)",
  "- `results/task004_cluster_rescued_cell_audit.csv`",
  "- `results/task004_author_vs_project_annotation_crosswalk.csv`",
  "- `figures/task004_marker_program_umap.pdf`, `figures/task004_marker_program_dotplot.pdf`, and `figures/task004_qc_annotation_review.pdf`",
  "- server-side annotated objects under `objects/task004_annotated/`",
  "",
  "## Hold point",
  "",
  "The corrected cell set and harmonized broad annotation layer are frozen at this Task 004 hold point. Do not execute Task 005 integration until Web GPT/user review."
)
writeLines(report_lines, report_path)
message("Task 004 annotation complete. Results: ", results_root)
