#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[[1]],mustWork=TRUE) else normalizePath(".",mustWork=TRUE)
res <- file.path(root,"results")
report <- file.path(root,"reports","task_004d_report.md")
source_qs <- file.path(root,"objects","merge","HCC_TA_8datasets_merged_review_v1.qs")
singlet_qs <- file.path(root,"objects","merge","HCC_TA_8datasets_singlets_v1.qs")
final_qs <- file.path(root,"objects","merge","HCC_TA_8datasets_singlets_broad_v1.qs")
needed <- c(
  file.path(res,"task004d_scdblfinder_validation.csv"),
  file.path(res,"task004d_scdblfinder_by_sample.csv"),
  file.path(res,"task004d_scdblfinder_by_dataset.csv"),
  file.path(res,"task004d_neutrophil_doublet_retention.csv"),
  file.path(res,"task004d_broad_annotation_validation.csv"),
  file.path(res,"task004d_broad_celltype_counts.csv"),
  file.path(res,"task004d_cluster_annotation.csv"),
  file.path(res,"task004d_sketch_usage.csv")
)
miss <- needed[!file.exists(needed)]
if(length(miss)) stop("Missing Task004d outputs: ",paste(miss,collapse=", "))
if(!all(file.exists(c(source_qs,singlet_qs,final_qs)))) stop("Missing one or more Task004d QS objects")

dval <- fread(file.path(res,"task004d_scdblfinder_validation.csv"))
samp <- fread(file.path(res,"task004d_scdblfinder_by_sample.csv"))
dset <- fread(file.path(res,"task004d_scdblfinder_by_dataset.csv"))
neut <- fread(file.path(res,"task004d_neutrophil_doublet_retention.csv"))
aval <- fread(file.path(res,"task004d_broad_annotation_validation.csv"))
broad <- fread(file.path(res,"task004d_broad_celltype_counts.csv"))
sketch <- fread(file.path(res,"task004d_sketch_usage.csv"))

if(dval$status[[1]]!="VALIDATED"||aval$status[[1]]!="VALIDATED") stop("Task004d validation not complete")
if(dval$n_cells_after[[1]]!=aval$n_cells[[1]]) stop("Singlet/broad cell counts disagree")

sha256 <- function(path){
  z <- system2("sha256sum",path,stdout=TRUE,stderr=TRUE)
  if(!length(z)) return(NA_character_)
  strsplit(z[[1]],"\\s+")[[1]][[1]]
}
fmt <- function(x) format(as.numeric(x),big.mark=",",scientific=FALSE,trim=TRUE)
pct <- function(x) sprintf("%.2f%%",100*as.numeric(x))

file_manifest <- data.table(
  role=c("Task004c input","Task004d singlets","Task004d singlets+broad"),
  path=c(source_qs,singlet_qs,final_qs),
  size_bytes=file.info(c(source_qs,singlet_qs,final_qs))$size,
  sha256=vapply(c(source_qs,singlet_qs,final_qs),sha256,character(1))
)
fwrite(file_manifest,file.path(res,"task004d_object_checksums.csv"))

broad_total <- broad[,.(n_cells=sum(n_cells)),by=broad_celltype]
setorder(broad_total,-n_cells)
dataset_lines <- paste(
  sprintf("- %s: %s/%s doublets (%s); %s singlets retained.",
          dset$dataset,fmt(dset$n_doublet),fmt(dset$n_cells_before),
          pct(dset$doublet_fraction),fmt(dset$n_retained)),
  collapse="\n"
)
broad_lines <- paste(
  sprintf("- %s: %s cells (%s).",
          broad_total$broad_celltype,fmt(broad_total$n_cells),
          pct(broad_total$n_cells/sum(broad_total$n_cells))),
  collapse="\n"
)
sketch_lines <- paste(
  sprintf("- %s: %s cells used for clustering of %s retained cells; sketch=%s; method=%s.",
          sketch$dataset,fmt(sketch$n_analysis_cells),fmt(sketch$n_cells),
          ifelse(sketch$used_sketch,"YES","NO"),sketch$sampling_method),
  collapse="\n"
)

sample_high <- samp[review_flag!="none"]
neut_before <- sum(neut$n_candidate_neutrophils_before,na.rm=TRUE)
neut_doublet <- sum(neut$n_candidate_neutrophils_removed_doublet,na.rm=TRUE)
neut_ret <- if(neut_before>0) 1-neut_doublet/neut_before else NA_real_

lines <- c(
  "# Task 004d report — scDblFinder filtering and corrected broad annotation",
  "",
  "## Status",
  if ("n_neutrophils" %in% names(aval) && aval$n_neutrophils[[1]] == 0)
    "PARTIAL — technically validated; no final Neutrophil labels, scientific review required"
  else "COMPLETED — VALIDATED",
  "",
  "## Reproducibility",
  paste0("- Runtime: ",R.version.string,
         "; Seurat ",as.character(packageVersion("Seurat")),
         "; scDblFinder ",dval$scDblFinder_version[[1]],
         "; SingleCellExperiment ",as.character(packageVersion("SingleCellExperiment")),
         "; BiocParallel ",as.character(packageVersion("BiocParallel")),
         "; SingleR ",as.character(packageVersion("SingleR")),
         "; xgboost ",as.character(packageVersion("xgboost")),"."),
  "- Phase A/B: `TASK004D_REUSE_CALL_TABLE=TRUE bash scripts/bash/task004de_doublet_and_broad_annotation.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas` (validated all-cell calls reused; no models rerun).",
  "- Phase C: `Rscript scripts/R/task004e_broad_annotation.R --project-root /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas --source-qs objects/merge/HCC_TA_8datasets_singlets_v1.qs --out-qs objects/merge/HCC_TA_8datasets_singlets_broad_v1.qs --shared-features results/task004c_shared_hgnc_features_8of8.txt --xue-map config/task004d_xue_author_to_broad.tsv --sketch-cells 50000 --projection-block 2000 --seed 40500`.",
  "- Final QC reconciliation: `Rscript scripts/R/task004d_reconcile_cluster_annotations.R /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`; reused the saved marker table and per-dataset Seurat caches, without rerunning doublet calls, clustering, or marker tests.",
  "- Broad QC figure-only re-render: `Rscript scripts/R/task004d_render_annotation_qc.R /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`; reused saved Seurat reductions and corrected cluster labels, with a compact bottom legend for the many nature_xue author labels.",
  "- Annotation parameters: LogNormalize scale factor 10,000; up to 3,000 HVGs after mitochondrial/ribosomal exclusion; 30-PC PCA; resolutions 0.4/0.8 with 0.8 used; leverage sketch target 50,000; projection blocks 2,000; seed 40500.",
  "- H5AD interoperability export was not generated: the existing Task004c exporter targets an in-memory standard RNA counts matrix, while this final QS preserves counts in the project's custom chunked representation; safe direct reuse was not verified. The validated deliverable is the QS object with raw counts and v2 metadata.",
  "",
  "## Input and output objects",
  paste0("- Input HGNC merged QS: ",source_qs," (",
         fmt(file_manifest[role=="Task004c input",size_bytes])," bytes; SHA256 ",
         file_manifest[role=="Task004c input",sha256],")."),
  paste0("- Singlet-only QS: ",singlet_qs," (",
         fmt(file_manifest[role=="Task004d singlets",size_bytes])," bytes; SHA256 ",
         file_manifest[role=="Task004d singlets",sha256],")."),
  paste0("- Final singlet+broad QS: ",final_qs," (",
         fmt(file_manifest[role=="Task004d singlets+broad",size_bytes])," bytes; SHA256 ",
         file_manifest[role=="Task004d singlets+broad",sha256],")."),
  "",
  "## scDblFinder",
  paste0("- scDblFinder version: ",dval$scDblFinder_version[[1]],"."),
  "- Detection unit: each project_sample_id independently.",
  if ("scDblFinder_rate_api" %in% names(dval) &&
      dval$scDblFinder_rate_api[[1]]=="legacy_dbr_explicit_per_1000_equivalent")
    "- Rate API compatibility: this scDblFinder release lacks dbr.per1k; per sample it received the explicit equivalent dbr=min(1, 0.008*n_cells/1000), preserving the approved 0.8%-per-1,000-cells expectation (not the package default). dbr.sd=NULL, iter=2, SerialParam()."
  else "- Parameters: dbr=NULL, dbr.per1k=0.008, dbr.sd=NULL, iter=2, SerialParam().",
  if("scDblFinder_calls_reused"%in%names(dval) && isTRUE(dval$scDblFinder_calls_reused[[1]]))
    "- Resume record: all-cell calls from the preceding 194-capture scoring pass were reused only after validating cell order, capture/sample identifiers, score completeness, package-rate API, and per-sample expected rates; models were not recomputed."
  else "- Resume record: all 194 capture-level models were scored during this execution.",
  "- Samples with >=500 cells used cluster-based mode; smaller samples used random artificial-doublet mode.",
  "- Seeds: 44000 + sample_index.",
  paste0("- Cells before: ",fmt(dval$n_cells_before[[1]]),"."),
  paste0("- Predicted doublets removed: ",fmt(dval$n_doublets_removed[[1]]),
         " (",pct(dval$n_doublets_removed[[1]]/dval$n_cells_before[[1]]),")."),
  paste0("- Singlets retained: ",fmt(dval$n_cells_after[[1]]),"."),
  paste0("- Preliminary Task004 neutrophil retention: ",
         ifelse(is.na(neut_ret),"NA",pct(neut_ret)),
         " (",fmt(neut_before-neut_doublet),"/",fmt(neut_before),")."),
  paste0("- Samples with doublet-rate review flags: ",nrow(sample_high),"."),
  "- CRA002308 source preprocessing reported doublet exclusion; the same project scDblFinder audit was still run because no per-cell known-doublet truth labels are available.",
  "- nature_xue is an author-processed object: all author labels were preserved, and scDblFinder was applied as an additional project-level computational filter.",
  "",
  "### Per-dataset doublet summary",
  dataset_lines,
  "",
  "## Broad annotation v2",
  "- Annotation was performed within each dataset without cross-dataset batch integration.",
  "- Feature universe started from the Task004c 8/8 shared HGNC set; mitochondrial and ribosomal HVGs were excluded from PCA features.",
  "- LogNormalize (scale factor 10,000), ~3,000 HVGs, ScaleData, 30-PC PCA, neighbors and clustering at resolutions 0.4 and 0.8 were used; resolution 0.8 is the working partition. After a documented full-cell OOM in nature_xue, Seurat leverage-score sketching (up to 50,000 cells) was used for that dataset; the same method is available to any >250,000-cell dataset exceeding the conservative memory guard. Sketch cluster labels were projected to all cells by block-wise PCA-centroid distance.",
  "- Cluster markers and canonical lineage programs were combined with Xue author-label pseudobulk reference transfer.",
  "- Neutrophil calls required coherent granulocyte evidence including at least one core marker; S100A8/S100A9 alone were insufficient.",
  "- No malignant-cell call and no cross-dataset integration were performed.",
  paste0("- Uncertain/Mixed: ",fmt(aval$n_uncertain[[1]]),
         " (",pct(aval$uncertain_fraction[[1]]),")."),
  paste0("- Xue broad author-label concordance (excluding unmapped author labels): ",
         ifelse(is.na(aval$xue_author_broad_concordance[[1]]),"NA",
                pct(aval$xue_author_broad_concordance[[1]])),"."),
  "",
  "### Final broad cell counts",
  broad_lines,
  "",
  "### Scalability / sketch use",
  sketch_lines,
  "",
  "## Key QC outputs",
  "- results/task004d_scdblfinder_cell_calls.csv.gz",
  "- results/task004d_scdblfinder_by_sample.csv",
  "- results/task004d_scdblfinder_by_dataset.csv",
  "- results/task004d_scdblfinder_by_dataset_tissue.csv",
  "- results/task004d_scdblfinder_overall.csv",
  "- results/task004d_scdblfinder_cell_qc_distributions.csv",
  "- results/task004d_neutrophil_doublet_retention.csv",
  "- results/task004d_neutrophil_retention_distributions.csv",
  "- results/task004d_broad_celltype_counts.csv",
  "- results/task004d_cluster_annotation.csv",
  "- results/task004d_cluster_markers.csv.gz",
  "- results/task004d_xue_author_vs_v2.csv",
  "- results/task004d_task004_vs_v2.csv",
  "- results/task004d_uncertain_clusters.csv",
  "- results/task004d_neutrophil_marker_coherence.csv",
  "- figures/task004d_scdblfinder_qc.pdf",
  "- figures/task004d_broad_annotation_qc.pdf (required review copy; generated source: figures/task004e/task004d_broad_annotation_qc.pdf)",
  "",
  "## Deviations / errors",
  if(nrow(sample_high)) paste0(
    "Doublet-rate review flags were retained for inspection and did not trigger threshold retuning: ",
    paste(sample_high$project_sample_id,collapse=", "),"."
  ) else "No sample-level doublet-rate review flag was triggered.",
  "Resolved runtime issues: the server's scDblFinder 1.16 API required the explicit legacy expected-rate equivalent; Seurat accessor validation was namespaced; saved gzip calls were streamed without optional R.utils; data.table singlet selection and grouped median output types were corrected. The validated first all-cell scoring pass was preserved and reused. A full-cell nature_xue attempt reached >135 GB RSS and was OOM-killed; the approved Seurat leverage-score sketch (50,000 cells) and block-wise PCA-centroid projection were then used. Sketch cells/counts/metadata were aligned explicitly by cell ID, and mapped cluster-label vectors were made unnamed before Seurat metadata assignment to prevent its names-as-cell-IDs overlap error. Seurat's optional `presto` acceleration was not installed; the standard Wilcoxon marker test completed.",
  "During final QC, a scoping error was found in the original core-marker audit: neutrophil-core and CD3 hits had been copied from the last cluster into every cluster row. Per-cluster hit counts and all downstream reference decisions were recomputed from the preserved marker table and Seurat analysis caches using the original thresholds; doublet calls, cluster definitions, and raw counts were not changed.",
  "No secondary UMI/nFeature doublet cutoff was used. No predicted doublet was rescued because of a preliminary cell-type label.",
  "",
  "## Hold point",
  "Task 005 remains paused. Review Task 004d QC and broad annotation before integration."
)
writeLines(lines,report)
cat("Final Task004d report written: ",report,"\n",sep="")
