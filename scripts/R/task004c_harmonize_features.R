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
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("Package 'qs' is required. Reuse the isolated Task 004b R library or install qs 0.27.3.")
}

cohort_summary_path <- arg_value(
  "cohort-summary",
  file.path(project_root, "results", "task004b_merge_review_cohort_summary.csv")
)
reference_dir <- arg_value(
  "reference-dir",
  file.path(project_root, "references", "task004c")
)
object_dir <- arg_value(
  "object-dir",
  file.path(project_root, "objects", "task004c_feature_harmonized")
)
results_dir <- arg_value("results-dir", file.path(project_root, "results"))
report_path <- arg_value(
  "report",
  file.path(project_root, "reports", "task_004c_report.md")
)

dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
hgnc_path <- file.path(reference_dir, "hgnc_complete_set.txt")
if (!file.exists(hgnc_path) || file.info(hgnc_path)$size < 1e6) {
  message("Downloading current HGNC complete set...")
  download.file(hgnc_url, hgnc_path, mode = "wb", quiet = FALSE)
}
if (!file.exists(hgnc_path) || file.info(hgnc_path)$size < 1e6) {
  stop("HGNC reference download failed or is unexpectedly small")
}
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
required_hgnc <- c(
  "hgnc_id", "symbol", "status", "locus_group", "locus_type",
  "alias_symbol", "prev_symbol", "ensembl_gene_id", "entrez_id"
)
missing_hgnc <- setdiff(required_hgnc, names(hgnc))
if (length(missing_hgnc)) {
  stop("HGNC reference lacks required fields: ", paste(missing_hgnc, collapse = ", "))
}
hgnc <- hgnc[status == "Approved" & !is.na(symbol) & nzchar(symbol)]
if (anyDuplicated(hgnc$symbol)) stop("HGNC approved symbol table contains duplicates")

approved_map <- hgnc[, .(
  canonical_symbol = symbol,
  hgnc_id,
  locus_group,
  locus_type,
  ensembl_gene_id,
  entrez_id
)]

make_unique_scalar_map <- function(dt, key_col, source_type) {
  z <- dt[!is.na(get(key_col)) & nzchar(get(key_col)),
          .(canonical_symbol = symbol, hgnc_id, locus_group, locus_type,
            key = as.character(get(key_col)))]
  if (!nrow(z)) return(list(unique = z, ambiguous = character()))
  counts <- z[, .(n_symbol = uniqueN(canonical_symbol)), by = key]
  ambiguous <- counts[n_symbol > 1L, key]
  u <- z[!key %in% ambiguous]
  if (nrow(u)) u <- unique(u, by = "key")
  u[, mapping_source := source_type]
  list(unique = u, ambiguous = ambiguous)
}

make_token_map <- function(dt, key_col, source_type) {
  x <- dt[!is.na(get(key_col)) & nzchar(get(key_col)),
          .(canonical_symbol = symbol, hgnc_id, locus_group, locus_type,
            tokens = as.character(get(key_col)))]
  if (!nrow(x)) return(list(unique = data.table(), ambiguous = character()))
  pieces <- strsplit(x$tokens, "\\|")
  z <- data.table(
    key = trimws(unlist(pieces, use.names = FALSE)),
    canonical_symbol = rep(x$canonical_symbol, lengths(pieces)),
    hgnc_id = rep(x$hgnc_id, lengths(pieces)),
    locus_group = rep(x$locus_group, lengths(pieces)),
    locus_type = rep(x$locus_type, lengths(pieces))
  )
  z <- z[!is.na(key) & nzchar(key)]
  counts <- z[, .(n_symbol = uniqueN(canonical_symbol)), by = key]
  ambiguous <- counts[n_symbol > 1L, key]
  u <- z[!key %in% ambiguous]
  if (nrow(u)) u <- unique(u, by = "key")
  u[, mapping_source := source_type]
  list(unique = u, ambiguous = ambiguous)
}

ens_map <- make_unique_scalar_map(hgnc, "ensembl_gene_id", "ensembl")
entrez_map <- make_unique_scalar_map(hgnc, "entrez_id", "entrez")
prev_map <- make_token_map(hgnc, "prev_symbol", "previous_symbol")
alias_map <- make_token_map(hgnc, "alias_symbol", "alias_symbol")

lookup_from <- function(keys, map) {
  if (!nrow(map)) {
    return(data.table(
      key = keys, canonical_symbol = NA_character_, hgnc_id = NA_character_,
      locus_group = NA_character_, locus_type = NA_character_,
      mapping_source = NA_character_
    ))
  }
  m <- map[match(keys, map$key)]
  data.table(
    key = keys,
    canonical_symbol = m$canonical_symbol,
    hgnc_id = m$hgnc_id,
    locus_group = m$locus_group,
    locus_type = m$locus_type,
    mapping_source = m$mapping_source
  )
}

map_features <- function(features, dataset_name) {
  out <- data.table(
    dataset = dataset_name,
    raw_feature = features,
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

  hit <- match(out$raw_feature, approved_map$canonical_symbol)
  ok <- !is.na(hit)
  if (any(ok)) {
    out[ok, canonical_symbol := approved_map$canonical_symbol[hit[ok]]]
    out[ok, hgnc_id := approved_map$hgnc_id[hit[ok]]]
    out[ok, locus_group := approved_map$locus_group[hit[ok]]]
    out[ok, locus_type := approved_map$locus_type[hit[ok]]]
    out[ok, mapping_status := "approved_symbol_exact"]
  }

  todo <- is.na(out$canonical_symbol) & out$feature_class == "ensembl"
  if (any(todo)) {
    m <- lookup_from(out$lookup_key[todo], ens_map$unique)
    rows <- which(todo)
    good <- !is.na(m$canonical_symbol)
    if (any(good)) {
      rr <- rows[good]
      out[rr, canonical_symbol := m$canonical_symbol[good]]
      out[rr, hgnc_id := m$hgnc_id[good]]
      out[rr, locus_group := m$locus_group[good]]
      out[rr, locus_type := m$locus_type[good]]
      out[rr, mapping_status := "ensembl_to_hgnc"]
    }
    amb <- out$lookup_key[rows] %in% ens_map$ambiguous
    if (any(amb)) out[rows[amb], mapping_status := "ambiguous_ensembl"]
  }

  todo <- is.na(out$canonical_symbol) & out$feature_class == "numeric_entrez_candidate"
  if (any(todo)) {
    m <- lookup_from(out$lookup_key[todo], entrez_map$unique)
    rows <- which(todo)
    good <- !is.na(m$canonical_symbol)
    if (any(good)) {
      rr <- rows[good]
      out[rr, canonical_symbol := m$canonical_symbol[good]]
      out[rr, hgnc_id := m$hgnc_id[good]]
      out[rr, locus_group := m$locus_group[good]]
      out[rr, locus_type := m$locus_type[good]]
      out[rr, mapping_status := "entrez_to_hgnc"]
    }
    amb <- out$lookup_key[rows] %in% entrez_map$ambiguous
    if (any(amb)) out[rows[amb], mapping_status := "ambiguous_entrez"]
  }

  todo <- is.na(out$canonical_symbol) & out$feature_class == "symbol_or_other"
  if (any(todo)) {
    m <- lookup_from(out$raw_feature[todo], prev_map$unique)
    rows <- which(todo)
    good <- !is.na(m$canonical_symbol)
    if (any(good)) {
      rr <- rows[good]
      out[rr, canonical_symbol := m$canonical_symbol[good]]
      out[rr, hgnc_id := m$hgnc_id[good]]
      out[rr, locus_group := m$locus_group[good]]
      out[rr, locus_type := m$locus_type[good]]
      out[rr, mapping_status := "previous_symbol_to_current"]
    }
    amb <- out$raw_feature[rows] %in% prev_map$ambiguous
    if (any(amb)) out[rows[amb], mapping_status := "ambiguous_previous_symbol"]
  }

  todo <- is.na(out$canonical_symbol) & out$feature_class == "symbol_or_other"
  if (any(todo)) {
    m <- lookup_from(out$raw_feature[todo], alias_map$unique)
    rows <- which(todo)
    good <- !is.na(m$canonical_symbol)
    if (any(good)) {
      rr <- rows[good]
      out[rr, canonical_symbol := m$canonical_symbol[good]]
      out[rr, hgnc_id := m$hgnc_id[good]]
      out[rr, locus_group := m$locus_group[good]]
      out[rr, locus_type := m$locus_type[good]]
      out[rr, mapping_status := "alias_symbol_to_current"]
    }
    amb <- out$raw_feature[rows] %in% alias_map$ambiguous
    if (any(amb)) out[rows[amb], mapping_status := "ambiguous_alias_symbol"]
  }

  out[is.na(mapping_status), mapping_status := "unmapped_preserved_in_source_only"]
  out[, mapped_hgnc := !is.na(canonical_symbol)]
  out[, ensembl_version_removed := feature_class == "ensembl" & raw_feature != lookup_key]
  out
}

collapse_to_hgnc <- function(counts, mapping, dataset_name) {
  keep <- which(mapping$mapped_hgnc)
  if (!length(keep)) stop("No HGNC-mapped features for ", dataset_name)
  x <- counts[keep, , drop = FALSE]
  canonical <- mapping$canonical_symbol[keep]

  if (!anyDuplicated(canonical)) {
    rownames(x) <- canonical
    return(x)
  }

  new_symbols <- unique(canonical)
  group_index <- match(canonical, new_symbols)
  aggregator <- sparseMatrix(
    i = group_index,
    j = seq_along(group_index),
    x = 1,
    dims = c(length(new_symbols), length(group_index))
  )
  y <- aggregator %*% x
  rownames(y) <- new_symbols
  colnames(y) <- colnames(x)
  if (anyDuplicated(rownames(y))) stop("Duplicate HGNC symbols remain after collapse: ", dataset_name)
  y
}

cohort_summary <- fread(cohort_summary_path)
expected_datasets <- c(
  "GSE282701", "GSE242889", "GSE326201", "GSE149614",
  "GSE299340", "CRA002308", "nature_xue", "in_house"
)
if (!setequal(cohort_summary$dataset, expected_datasets)) {
  stop("Dataset mismatch in Task 004b cohort summary")
}

mapping_list <- list()
audit_rows <- list()
presence_list <- list()
object_rows <- list()

for (dataset_name in expected_datasets) {
  source_path <- cohort_summary[dataset == dataset_name, cohort_object_path][[1L]]
  message("Harmonising features: ", dataset_name)
  obj <- readRDS(source_path)
  if (!inherits(obj, "Seurat") || !"RNA" %in% Assays(obj)) stop("Invalid cohort object: ", source_path)
  if (!identical(SeuratObject::Layers(obj[["RNA"]]), "counts")) stop("Expected exactly one counts layer in ", dataset_name)

  counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
  original_features <- rownames(counts)
  if (anyDuplicated(original_features)) stop("Source object has duplicate feature names: ", dataset_name)

  mp <- map_features(original_features, dataset_name)
  mapping_list[[dataset_name]] <- mp
  mapped_unique <- unique(mp$canonical_symbol[mp$mapped_hgnc])

  audit_rows[[dataset_name]] <- data.table(
    dataset = dataset_name,
    source_features = length(original_features),
    mapped_feature_rows = sum(mp$mapped_hgnc),
    mapped_unique_hgnc_symbols = length(mapped_unique),
    mapping_rate_feature_rows = sum(mp$mapped_hgnc) / length(original_features),
    approved_symbol_exact = sum(mp$mapping_status == "approved_symbol_exact"),
    ensembl_to_hgnc = sum(mp$mapping_status == "ensembl_to_hgnc"),
    entrez_to_hgnc = sum(mp$mapping_status == "entrez_to_hgnc"),
    previous_symbol_to_current = sum(mp$mapping_status == "previous_symbol_to_current"),
    alias_symbol_to_current = sum(mp$mapping_status == "alias_symbol_to_current"),
    ambiguous = sum(grepl("^ambiguous_", mp$mapping_status)),
    unmapped = sum(mp$mapping_status == "unmapped_preserved_in_source_only"),
    ensembl_version_removed = sum(mp$ensembl_version_removed),
    duplicate_rows_collapsed = sum(mp$mapped_hgnc) - length(mapped_unique)
  )

  presence_list[[dataset_name]] <- data.table(
    canonical_symbol = mapped_unique,
    dataset = dataset_name
  )

  md <- obj[[]]
  if ("nCount_RNA" %in% colnames(md)) {
    md$source_nCount_RNA_before_feature_harmonisation <- md$nCount_RNA
  }
  if ("nFeature_RNA" %in% colnames(md)) {
    md$source_nFeature_RNA_before_feature_harmonisation <- md$nFeature_RNA
  }

  hcounts <- collapse_to_hgnc(counts, mp, dataset_name)
  hobj <- CreateSeuratObject(
    counts = hcounts,
    assay = "RNA",
    project = paste0(dataset_name, "_HGNC"),
    meta.data = md
  )
  hobj$feature_harmonisation <- "HGNC_approved_symbol"
  hobj@misc$feature_harmonisation <- list(
    task = "task_004c",
    reference = "HGNC complete set",
    reference_url = hgnc_url,
    reference_sha256 = sha256,
    reference_retrieved_date = format(Sys.Date(), "%Y-%m-%d"),
    original_feature_count = length(original_features),
    mapped_feature_rows = sum(mp$mapped_hgnc),
    harmonised_feature_count = nrow(hobj),
    mapping_rule = "approved symbol > Ensembl > Entrez > unique previous symbol > unique alias symbol",
    duplicate_rule = "sum raw counts for source features mapping to the same approved HGNC symbol",
    unmapped_rule = "retain in source object only; exclude from cross-cohort human integration"
  )

  out_path <- file.path(object_dir, paste0(dataset_name, "_HGNC_harmonized.qs"))
  detected_cores <- parallel::detectCores()
  if (is.na(detected_cores) || detected_cores < 1L) detected_cores <- 1L
  nthreads <- max(1L, min(8L, detected_cores))
  expected_cells <- ncol(hobj)
  expected_features <- nrow(hobj)
  qs::qsave(hobj, out_path, preset = "high", check_hash = TRUE, nthreads = nthreads)

  # Free the source and newly materialised matrices before reload validation.
  # This is important for nature_xue and avoids holding two full Seurat objects
  # plus the validation copy in RAM at the same time.
  rm(obj, counts, md, hcounts, hobj)
  gc(verbose = FALSE)

  chk <- qs::qread(out_path, use_alt_rep = FALSE, nthreads = nthreads)
  if (!inherits(chk, "Seurat")) stop("Reload failed for ", dataset_name)
  if (ncol(chk) != expected_cells) stop("Cell count changed after harmonisation: ", dataset_name)
  if (nrow(chk) != expected_features) stop("Feature count changed after harmonisation: ", dataset_name)
  if (anyDuplicated(rownames(chk))) stop("Duplicate harmonised features after reload: ", dataset_name)
  if (!all(rownames(chk) %in% hgnc$symbol)) stop("Non-HGNC feature found in harmonised object: ", dataset_name)

  object_rows[[dataset_name]] <- data.table(
    dataset = dataset_name,
    source_object = source_path,
    harmonised_object = out_path,
    n_cells = ncol(chk),
    n_features_hgnc = nrow(chk),
    file_size_bytes = file.info(out_path)$size,
    validation = "VALIDATED"
  )

  rm(chk, mp)
  gc(verbose = FALSE)
}

mapping_dt <- rbindlist(mapping_list, use.names = TRUE, fill = TRUE)
audit_dt <- rbindlist(audit_rows, use.names = TRUE, fill = TRUE)
object_dt <- rbindlist(object_rows, use.names = TRUE, fill = TRUE)
presence_long <- rbindlist(presence_list, use.names = TRUE)

presence <- dcast(
  unique(presence_long)[, present := TRUE],
  canonical_symbol ~ dataset,
  value.var = "present",
  fill = FALSE
)
dataset_cols <- intersect(expected_datasets, names(presence))
presence[, n_datasets := rowSums(.SD), .SDcols = dataset_cols]
setcolorder(presence, c("canonical_symbol", "n_datasets", expected_datasets))
setorder(presence, -n_datasets, canonical_symbol)

shared8 <- presence[n_datasets == 8L, canonical_symbol]
shared7plus <- presence[n_datasets >= 7L, canonical_symbol]

if (length(shared8) < 10000L) {
  stop("Strict 8/8 HGNC shared feature set unexpectedly small: ", length(shared8))
}

pairwise <- rbindlist(lapply(seq_along(expected_datasets), function(i) {
  rbindlist(lapply(seq_along(expected_datasets), function(j) {
    a <- expected_datasets[[i]]
    b <- expected_datasets[[j]]
    data.table(
      dataset_1 = a,
      dataset_2 = b,
      shared_hgnc_symbols = sum(presence[[a]] & presence[[b]])
    )
  }))
}))

fwrite(mapping_dt, file.path(results_dir, "task004c_feature_mapping.csv"))
fwrite(audit_dt, file.path(results_dir, "task004c_feature_audit_by_cohort.csv"))
fwrite(presence, file.path(results_dir, "task004c_hgnc_feature_presence.csv"))
fwrite(pairwise, file.path(results_dir, "task004c_pairwise_feature_overlap.csv"))
fwrite(object_dt, file.path(results_dir, "task004c_harmonized_objects.csv"))
writeLines(shared8, file.path(results_dir, "task004c_shared_hgnc_features_8of8.txt"))
writeLines(shared7plus, file.path(results_dir, "task004c_shared_hgnc_features_7plus.txt"))
status_counts <- mapping_dt[, .N, by = .(dataset, mapping_status)]
fwrite(status_counts, file.path(results_dir, "task004c_mapping_status_counts.csv"))

report <- c(
  "# Task 004c report - HGNC feature harmonisation",
  "",
  "## Status",
  "",
  "COMPLETED",
  "",
  "## Strategy",
  "",
  "- Source Task 004/004b objects were not modified.",
  "- Current HGNC complete-set data were downloaded and checksum-recorded.",
  "- Mapping priority: approved HGNC symbol > Ensembl gene ID (version removed) > Entrez ID > unique previous symbol > unique alias symbol.",
  "- Ambiguous aliases/previous symbols were not force-mapped.",
  "- Source features mapping to the same approved HGNC symbol were collapsed by summing raw counts.",
  "- Unmapped/custom/non-HGNC features remain available in source objects but are excluded from the cross-cohort human integration feature universe.",
  "",
  "## Summary",
  "",
  paste0("Strict HGNC shared genes present in all 8 cohorts: ", format(length(shared8), big.mark = ",")),
  paste0("HGNC genes present in at least 7/8 cohorts: ", format(length(shared7plus), big.mark = ",")),
  "",
  "Task 005 remains paused pending review of the harmonisation audit."
)
writeLines(report, report_path)
message("Task 004c COMPLETED")
message("Strict 8/8 shared HGNC symbols: ", length(shared8))
