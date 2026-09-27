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
r_lib <- arg_value("r-lib", file.path(project_root, ".task004b_Rlib"))
.libPaths(c(r_lib, .libPaths()))
if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required")

source_qs <- arg_value(
  "source-qs",
  file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.qs")
)
out_qs <- arg_value(
  "out-qs",
  file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.qs.task004c_tmp")
)
final_qs <- arg_value(
  "final-qs",
  file.path(project_root, "objects", "HCC_TA_8datasets_merged_review_v1.qs")
)
cohort_summary_path <- arg_value(
  "cohort-summary",
  file.path(project_root, "results", "task004b_merge_review_cohort_summary.csv")
)
reference_dir <- file.path(project_root, "references", "task004c")
results_dir <- file.path(project_root, "results")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_qs), recursive = TRUE, showWarnings = FALSE)

hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
hgnc_path <- file.path(reference_dir, "hgnc_complete_set.txt")
required_hgnc <- c(
  "hgnc_id", "symbol", "status", "locus_group", "locus_type",
  "alias_symbol", "prev_symbol", "ensembl_gene_id", "entrez_id"
)
hgnc_tmp <- paste0(hgnc_path, ".partial")
message("Downloading current HGNC complete set to a temporary file...")
unlink(hgnc_tmp)
options(timeout = 1800)
curl_status <- suppressWarnings(system2(
  "curl",
  args = c(
    "--fail", "--location", "--retry", "3", "--connect-timeout", "30",
    "--max-time", "1800", "--silent", "--show-error", "--output",
    shQuote(hgnc_tmp), shQuote(hgnc_url)
  ),
  stdout = TRUE,
  stderr = TRUE
))
if (!is.null(attr(curl_status, "status")) && attr(curl_status, "status") != 0L) {
  stop("HGNC reference download failed: ", paste(curl_status, collapse = "\n"))
}
if (!file.exists(hgnc_tmp) || file.info(hgnc_tmp)$size < 1e7) {
  stop("HGNC reference download is missing or unexpectedly small")
}
hgnc_header <- names(fread(hgnc_tmp, nrows = 0L))
if (length(setdiff(required_hgnc, hgnc_header))) {
  stop("Downloaded HGNC reference has an unexpected header")
}
if (!file.rename(hgnc_tmp, hgnc_path)) stop("Could not atomically install downloaded HGNC reference")
if (!file.exists(hgnc_path) || file.info(hgnc_path)$size < 1e7) stop("HGNC reference download failed")

sha256 <- tryCatch(
  strsplit(system2("sha256sum", hgnc_path, stdout = TRUE), "\\s+")[[1L]][[1L]],
  error = function(e) NA_character_
)
writeLines(
  c(
    paste0("retrieved_date=", format(Sys.Date(), "%Y-%m-%d")),
    paste0("source_url=", hgnc_url),
    paste0("sha256=", sha256)
  ),
  file.path(reference_dir, "hgnc_reference_provenance.txt")
)

hgnc <- fread(hgnc_path, na.strings = c("", "NA"))
miss <- setdiff(required_hgnc, names(hgnc))
if (length(miss)) stop("HGNC reference lacks: ", paste(miss, collapse = ", "))
hgnc <- hgnc[status == "Approved" & !is.na(symbol) & nzchar(symbol)]
if (anyDuplicated(hgnc$symbol)) stop("HGNC approved symbols are duplicated")

approved <- hgnc[, .(
  canonical_symbol = symbol,
  hgnc_id,
  locus_group,
  locus_type,
  ensembl_gene_id,
  entrez_id
)]

make_scalar_map <- function(dt, key_col) {
  z <- dt[!is.na(get(key_col)) & nzchar(get(key_col)),
          .(lookup_key = as.character(get(key_col)), canonical_symbol = symbol,
            hgnc_id, locus_group, locus_type)]
  if (!nrow(z)) return(list(unique = z, ambiguous = character()))
  amb <- z[, .(n = uniqueN(canonical_symbol)), by = lookup_key][n > 1L, lookup_key]
  u <- unique(z[!lookup_key %in% amb], by = "lookup_key")
  setnames(u, "lookup_key", "key")
  list(unique = u, ambiguous = amb)
}

make_token_map <- function(dt, key_col) {
  x <- dt[!is.na(get(key_col)) & nzchar(get(key_col)),
          .(canonical_symbol = symbol, hgnc_id, locus_group, locus_type,
            tokens = as.character(get(key_col)))]
  if (!nrow(x)) return(list(unique = data.table(), ambiguous = character()))
  sp <- strsplit(x$tokens, "\\|")
  z <- data.table(
    lookup_key = trimws(unlist(sp, use.names = FALSE)),
    canonical_symbol = rep(x$canonical_symbol, lengths(sp)),
    hgnc_id = rep(x$hgnc_id, lengths(sp)),
    locus_group = rep(x$locus_group, lengths(sp)),
    locus_type = rep(x$locus_type, lengths(sp))
  )
  z <- z[!is.na(lookup_key) & nzchar(lookup_key)]
  amb <- z[, .(n = uniqueN(canonical_symbol)), by = lookup_key][n > 1L, lookup_key]
  u <- unique(z[!lookup_key %in% amb], by = "lookup_key")
  setnames(u, "lookup_key", "key")
  list(unique = u, ambiguous = amb)
}

ens_map <- make_scalar_map(hgnc, "ensembl_gene_id")
entrez_map <- make_scalar_map(hgnc, "entrez_id")
prev_map <- make_token_map(hgnc, "prev_symbol")
alias_map <- make_token_map(hgnc, "alias_symbol")

fill_from_map <- function(out, rows, keys, mp, status, ambiguous_status) {
  if (!length(rows)) return(out)
  hit <- match(keys, mp$unique$key)
  good <- !is.na(hit)
  if (any(good)) {
    rr <- rows[good]
    mm <- hit[good]
    out[rr, canonical_symbol := mp$unique$canonical_symbol[mm]]
    out[rr, hgnc_id := mp$unique$hgnc_id[mm]]
    out[rr, locus_group := mp$unique$locus_group[mm]]
    out[rr, locus_type := mp$unique$locus_type[mm]]
    out[rr, mapping_status := status]
  }
  amb <- keys %in% mp$ambiguous
  if (any(amb)) out[rows[amb], mapping_status := ambiguous_status]
  out
}

map_features <- function(features, dataset = "merged_union") {
  out <- data.table(
    dataset = dataset,
    raw_feature = as.character(features),
    feature_index = seq_along(features)
  )
  out[, feature_class := fifelse(
    grepl("^ENSG[0-9]+(\\.[0-9]+)?$", raw_feature), "ensembl",
    fifelse(grepl("^[0-9]+$", raw_feature), "numeric_entrez_candidate", "symbol_or_other")
  )]
  out[, lookup_key := raw_feature]
  out[feature_class == "ensembl", lookup_key := sub("\\.[0-9]+$", "", raw_feature)]
  out[, canonical_symbol := NA_character_]
  out[, hgnc_id := NA_character_]
  out[, locus_group := NA_character_]
  out[, locus_type := NA_character_]
  out[, mapping_status := NA_character_]

  hit <- match(out$raw_feature, approved$canonical_symbol)
  ok <- !is.na(hit)
  if (any(ok)) {
    out[ok, canonical_symbol := approved$canonical_symbol[hit[ok]]]
    out[ok, hgnc_id := approved$hgnc_id[hit[ok]]]
    out[ok, locus_group := approved$locus_group[hit[ok]]]
    out[ok, locus_type := approved$locus_type[hit[ok]]]
    out[ok, mapping_status := "approved_symbol_exact"]
  }

  rows <- which(is.na(out$canonical_symbol) & is.na(out$mapping_status) &
                  out$feature_class == "ensembl")
  out <- fill_from_map(out, rows, out$lookup_key[rows], ens_map,
                       "ensembl_to_hgnc", "ambiguous_ensembl")
  rows <- which(is.na(out$canonical_symbol) & is.na(out$mapping_status) &
                  out$feature_class == "numeric_entrez_candidate")
  out <- fill_from_map(out, rows, out$lookup_key[rows], entrez_map,
                       "entrez_to_hgnc", "ambiguous_entrez")
  rows <- which(is.na(out$canonical_symbol) & is.na(out$mapping_status) &
                  out$feature_class == "symbol_or_other")
  out <- fill_from_map(out, rows, out$raw_feature[rows], prev_map,
                       "previous_symbol_to_current", "ambiguous_previous_symbol")
  rows <- which(is.na(out$canonical_symbol) & is.na(out$mapping_status) &
                  out$feature_class == "symbol_or_other")
  out <- fill_from_map(out, rows, out$raw_feature[rows], alias_map,
                       "alias_symbol_to_current", "ambiguous_alias_symbol")

  out[is.na(mapping_status), mapping_status := "unmapped_preserved_in_v1_only"]
  out[, mapped_hgnc := !is.na(canonical_symbol)]
  out[, ensembl_version_removed := feature_class == "ensembl" & raw_feature != lookup_key]
  out
}

# Methods needed for the custom Task 004b counts layer after qread.
dim.task004b_chunked_counts <- function(x) c(length(x$feature_names), length(x$cell_names))
dimnames.task004b_chunked_counts <- function(x) list(x$feature_names, x$cell_names)

message("Reading merged QS: ", source_qs)
obj <- qs::qread(source_qs, use_alt_rep = FALSE, nthreads = 8L)
if (!inherits(obj, "Seurat")) stop("Source QS is not Seurat")
if (ncol(obj) != 1490852L || nrow(obj) != 68394L) stop("Source QS dimensions mismatch")
if (!identical(SeuratObject::Layers(obj[["RNA"]]), "counts")) stop("Source QS must have one counts layer")

counts_layer <- obj[["RNA"]]@layers[["counts"]]
if (!inherits(counts_layer, "task004b_chunked_counts")) {
  stop("Expected task004b_chunked_counts in merged QS")
}
old_features <- counts_layer$feature_names
if (length(old_features) != 68394L || anyDuplicated(old_features)) stop("Source feature names invalid")

input_meta <- obj[[]]
required_metadata <- c(
  "dataset", "sample_id", "patient_id", "tissue", "project_sample_id",
  "project_patient_id", "nCount_RNA", "nFeature_RNA"
)
missing_metadata <- setdiff(required_metadata, colnames(input_meta))
if (length(missing_metadata)) {
  stop("Source QS lacks required metadata: ", paste(missing_metadata, collapse = ", "))
}
if (!identical(rownames(input_meta), colnames(obj))) stop("Source metadata row order differs from cell order")
if (anyDuplicated(colnames(obj))) stop("Source QS contains duplicated cell IDs")
required_missing <- vapply(required_metadata, function(nm) {
  x <- input_meta[[nm]]
  miss <- is.na(x)
  if (is.character(x) || is.factor(x)) miss <- miss | !nzchar(trimws(as.character(x)))
  sum(miss)
}, integer(1L))
if (any(required_missing > 0L)) {
  stop("Missing values in required metadata: ",
       paste(names(required_missing)[required_missing > 0L], required_missing[required_missing > 0L], collapse = ", "))
}
if (data.table::uniqueN(input_meta$dataset) != 8L ||
    data.table::uniqueN(input_meta$project_sample_id) != 194L ||
    data.table::uniqueN(input_meta$project_patient_id) != 132L) {
  stop("Source QS dataset/sample/patient identifiers differ from the validated Task 004b counts")
}
if (sum(input_meta$tissue == "Tumor") != 1039293L ||
    sum(input_meta$tissue == "Adjacent") != 451559L) {
  stop("Source QS tissue cell counts differ from the validated Task 004b counts")
}
metadata_missingness <- data.table(
  field = colnames(input_meta),
  class = vapply(input_meta, function(x) paste(class(x), collapse = ";"), character(1L)),
  n_missing = vapply(input_meta, function(x) {
    miss <- is.na(x)
    if (is.character(x) || is.factor(x)) miss <- miss | !nzchar(trimws(as.character(x)))
    sum(miss)
  }, integer(1L))
)
metadata_missingness[, fraction_missing := n_missing / nrow(input_meta)]
fwrite(metadata_missingness, file.path(results_dir, "task004c_input_metadata_missingness.csv"))
rm(input_meta, metadata_missingness)

mapping <- map_features(old_features)
mapping[, canonical_index := NA_integer_]
new_symbols <- unique(mapping$canonical_symbol[mapping$mapped_hgnc])
mapping[mapped_hgnc == TRUE, canonical_index := match(canonical_symbol, new_symbols)]

canonical_meta <- unique(mapping[mapped_hgnc == TRUE, .(
  canonical_symbol, hgnc_id, locus_group, locus_type, canonical_index
)], by = "canonical_symbol")
setorder(canonical_meta, canonical_index)
if (!identical(canonical_meta$canonical_symbol, new_symbols)) stop("Canonical feature ordering mismatch")
if (anyDuplicated(new_symbols)) stop("Canonical symbols duplicated")

fwrite(mapping, file.path(results_dir, "task004c_merged_feature_mapping.csv"))
fwrite(canonical_meta, file.path(results_dir, "task004c_hgnc_var.csv"))
fwrite(mapping[, .N, by = mapping_status][order(-N)],
       file.path(results_dir, "task004c_mapping_status_counts.csv"))

mapped_rows <- which(mapping$mapped_hgnc)
aggregator <- sparseMatrix(
  i = mapping$canonical_index[mapped_rows],
  j = mapped_rows,
  x = 1,
  dims = c(length(new_symbols), length(old_features))
)
rownames(aggregator) <- new_symbols
colnames(aggregator) <- old_features

old_total_counts <- 0
new_total_counts <- 0
new_chunks <- vector("list", length(counts_layer$chunks))
source_chunk_count <- length(counts_layer$chunks)
hgnc_nCount <- numeric(ncol(obj))
hgnc_nFeature <- integer(ncol(obj))
cell_cursor <- 0L

message("Harmonising merged QS counts in ", length(counts_layer$chunks), " sparse chunks...")
for (k in seq_len(source_chunk_count)) {
  message("  chunk ", k, "/", source_chunk_count)
  x <- counts_layer$chunks[[k]]
  if (!inherits(x, "sparseMatrix")) x <- as(x, "dgCMatrix")
  if (nrow(x) != length(old_features)) stop("Chunk feature dimension mismatch")
  old_total_counts <- old_total_counts + sum(x@x)

  y <- aggregator %*% x
  y <- as(y, "dgCMatrix")
  rownames(y) <- new_symbols
  colnames(y) <- colnames(x)
  y <- drop0(y)

  new_total_counts <- new_total_counts + sum(y@x)
  cols <- seq_len(ncol(y)) + cell_cursor
  hgnc_nCount[cols] <- Matrix::colSums(y)
  hgnc_nFeature[cols] <- diff(y@p)
  cell_cursor <- cell_cursor + ncol(y)

  new_chunks[[k]] <- y
  # Preserve list positions so later chunks retain their original indices.
  counts_layer$chunks[k] <- list(NULL)
  # Replace the source layer in the Seurat object as well, releasing the
  # processed sparse block instead of retaining a second full source copy.
  assay <- obj[["RNA"]]
  assay@layers$counts <- counts_layer
  obj@assays$RNA <- assay
  rm(assay)
  rm(x, y)
  gc(verbose = FALSE)
}
if (cell_cursor != ncol(obj)) stop("Processed cell count mismatch")

offsets <- c(0L, cumsum(vapply(new_chunks, ncol, integer(1L))))
new_counts <- structure(
  list(
    chunks = new_chunks,
    feature_names = new_symbols,
    cell_names = counts_layer$cell_names,
    offsets = offsets
  ),
  class = "task004b_chunked_counts"
)

assay <- obj[["RNA"]]
assay@layers$counts <- new_counts
assay@features <- SeuratObject:::LogMap(new_symbols)
assay@meta.data <- data.frame(row.names = new_symbols)
obj@assays$RNA <- assay
obj$source_nCount_RNA_before_feature_harmonisation <- obj$nCount_RNA
obj$source_nFeature_RNA_before_feature_harmonisation <- obj$nFeature_RNA
obj$hgnc_nCount_RNA <- hgnc_nCount
obj$hgnc_nFeature_RNA <- hgnc_nFeature
obj$feature_harmonisation <- "HGNC_approved_symbol"
obj@misc$feature_harmonisation <- list(
  task = "task_004c",
  source_object = source_qs,
  source_feature_count = length(old_features),
  harmonised_feature_count = length(new_symbols),
  reference = "HGNC complete set",
  reference_url = hgnc_url,
  reference_sha256 = sha256,
  reference_retrieved_date = format(Sys.Date(), "%Y-%m-%d"),
  mapping_rule = "approved symbol > version-stripped Ensembl > Entrez > unique previous symbol > unique alias",
  duplicate_rule = "raw counts summed for source features mapping to the same HGNC-approved symbol",
  unmapped_rule = "excluded from the overwritten merged analysis object; recoverable only from upstream cohort/source objects or retained Task 004b checkpoint",
  old_total_counts = old_total_counts,
  hgnc_total_counts = new_total_counts,
  count_retention_fraction = new_total_counts / old_total_counts
)

# Cheap source-cohort presence audit only; no cohort harmonised objects are written.
cohort_summary <- fread(cohort_summary_path)
expected_datasets <- c(
  "GSE282701", "GSE242889", "GSE326201", "GSE149614",
  "GSE299340", "CRA002308", "nature_xue", "in_house"
)
presence_list <- list()
cohort_audit <- list()
for (d in expected_datasets) {
  p <- cohort_summary[dataset == d, cohort_object_path][[1L]]
  o <- readRDS(p)
  f <- rownames(SeuratObject::LayerData(o[["RNA"]], layer = "counts"))
  mp <- map_features(f, d)
  syms <- unique(mp$canonical_symbol[mp$mapped_hgnc])
  presence_list[[d]] <- data.table(canonical_symbol = syms, dataset = d)
  cohort_audit[[d]] <- data.table(
    dataset = d,
    source_features = length(f),
    mapped_feature_rows = sum(mp$mapped_hgnc),
    mapped_unique_hgnc_symbols = length(syms),
    mapping_rate = sum(mp$mapped_hgnc) / length(f),
    ambiguous = sum(grepl("^ambiguous_", mp$mapping_status)),
    unmapped = sum(mp$mapping_status == "unmapped_preserved_in_v1_only")
  )
  rm(o, mp)
  gc(verbose = FALSE)
}
presence_long <- rbindlist(presence_list)
presence <- dcast(
  unique(presence_long)[, present := TRUE],
  canonical_symbol ~ dataset,
  value.var = "present",
  fill = FALSE
)
presence[, n_datasets := rowSums(.SD), .SDcols = expected_datasets]
setcolorder(presence, c("canonical_symbol", "n_datasets", expected_datasets))
setorder(presence, -n_datasets, canonical_symbol)
shared8 <- presence[n_datasets == 8L, canonical_symbol]
shared7 <- presence[n_datasets >= 7L, canonical_symbol]
if (length(shared8) < 10000L) stop("Unexpectedly small strict 8/8 HGNC set: ", length(shared8))

fwrite(rbindlist(cohort_audit), file.path(results_dir, "task004c_feature_audit_by_cohort.csv"))
fwrite(presence, file.path(results_dir, "task004c_hgnc_feature_presence.csv"))
writeLines(shared8, file.path(results_dir, "task004c_shared_hgnc_features_8of8.txt"))
writeLines(shared7, file.path(results_dir, "task004c_shared_hgnc_features_7plus.txt"))

message("Writing harmonised merged QS: ", out_qs)
qs::qsave(obj, out_qs, preset = "high", check_hash = TRUE, nthreads = 8L)

expected_cells <- ncol(obj)
expected_features <- nrow(obj)
expected_cell_ids <- colnames(obj)
rm(obj, assay, new_counts, new_chunks, aggregator, counts_layer)
gc(verbose = FALSE)

message("Reload-validating harmonised QS...")
chk <- qs::qread(out_qs, use_alt_rep = FALSE, nthreads = 8L)
if (!inherits(chk, "Seurat")) stop("Harmonised QS reload is not Seurat")
if (ncol(chk) != expected_cells || nrow(chk) != expected_features) stop("Harmonised QS dimensions mismatch")
if (!identical(colnames(chk), expected_cell_ids)) stop("Harmonised QS cell order changed")
if (!identical(rownames(chk), new_symbols)) stop("Harmonised QS feature order changed")
if (anyDuplicated(rownames(chk))) stop("Harmonised QS has duplicate gene symbols")
if (!all(rownames(chk) %in% hgnc$symbol)) stop("Harmonised QS contains non-HGNC row names")
if (!identical(SeuratObject::Layers(chk[["RNA"]]), "counts")) stop("Harmonised QS RNA layer mismatch")
if (!inherits(chk[["RNA"]]@layers[["counts"]], "task004b_chunked_counts")) stop("Harmonised QS counts class mismatch")
reloaded_counts <- chk[["RNA"]]@layers[["counts"]]
reloaded_total_counts <- 0
if (length(reloaded_counts$chunks) != source_chunk_count) stop("Harmonised QS chunk count changed")
for (k in seq_len(source_chunk_count)) {
  chunk <- reloaded_counts$chunks[[k]]
  if (!inherits(chunk, "sparseMatrix")) stop("Reloaded QS chunk is not sparse: ", k)
  reloaded_total_counts <- reloaded_total_counts + sum(chunk@x)
}
if (!identical(as.numeric(reloaded_total_counts), as.numeric(new_total_counts))) {
  stop("Reloaded QS count total differs from the transformation total")
}
if (!identical(as.numeric(chk$hgnc_nCount_RNA), as.numeric(hgnc_nCount))) {
  stop("Reloaded QS hgnc_nCount_RNA metadata mismatch")
}
if (!identical(as.integer(chk$hgnc_nFeature_RNA), as.integer(hgnc_nFeature))) {
  stop("Reloaded QS hgnc_nFeature_RNA metadata mismatch")
}

validation <- data.table(
  status = "VALIDATED",
  source_qs = source_qs,
  harmonised_qs_temp = out_qs,
  harmonised_qs_final = final_qs,
  qs_size_bytes = file.info(out_qs)$size,
  n_cells = ncol(chk),
  source_features = length(old_features),
  hgnc_features = nrow(chk),
  strict_shared_hgnc_8of8 = length(shared8),
  shared_hgnc_7plus = length(shared7),
  duplicated_gene_symbols = anyDuplicated(rownames(chk)),
  old_total_counts = old_total_counts,
  hgnc_total_counts = new_total_counts,
  count_retention_fraction = new_total_counts / old_total_counts,
  reference_sha256 = sha256
)
fwrite(validation, file.path(results_dir, "task004c_qs_validation.csv"))
message("Task 004c QS harmonisation VALIDATED")
