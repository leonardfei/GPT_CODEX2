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
  sprintf("- %s: %s cells used for clustering of %s retained cells; sketch=%s.",
          sketch$dataset,fmt(sketch$n_analysis_cells),fmt(sketch$n_cells),
          ifelse(sketch$used_sketch,"YES","NO")),
  collapse="\n"
)

sample_high <- samp[review_flag!="none"]
neut_before <- sum(neut$n_neutrophil_before,na.rm=TRUE)
neut_doublet <- sum(neut$n_neutrophil_doublet,na.rm=TRUE)
neut_ret <- if(neut_before>0) 1-neut_doublet/neut_before else NA_real_

lines <- c(
  "# Task 004d report — scDblFinder filtering and corrected broad annotation",
  "",
  "## Status",
  "COMPLETED — VALIDATED",
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
  "- Parameters: dbr=NULL, dbr.per1k=0.008, dbr.sd=NULL, iter=2, SerialParam().",
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
  "",
  "### Per-dataset doublet summary",
  dataset_lines,
  "",
  "## Broad annotation v2",
  "- Annotation was performed within each dataset without cross-dataset batch integration.",
  "- Feature universe started from the Task004c 8/8 shared HGNC set; mitochondrial and ribosomal HVGs were excluded from PCA features.",
  "- LogNormalize (scale factor 10,000), ~3,000 HVGs, ScaleData, 30-PC PCA, neighbors and clustering at resolutions 0.4 and 0.8 were used; resolution 0.8 is the working partition.",
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
  "- results/task004d_neutrophil_doublet_retention.csv",
  "- results/task004d_broad_celltype_counts.csv",
  "- results/task004d_cluster_annotation.csv",
  "- results/task004d_cluster_markers.csv.gz",
  "- results/task004d_xue_author_vs_v2.csv",
  "- results/task004d_task004_vs_v2.csv",
  "- results/task004d_uncertain_clusters.csv",
  "- results/task004d_neutrophil_marker_coherence.csv",
  "- figures/task004d_scdblfinder_qc.pdf",
  "- figures/task004e/task004d_broad_annotation_qc.pdf",
  "",
  "## Deviations / errors",
  if(nrow(sample_high)) paste0(
    "Doublet-rate review flags were retained for inspection and did not trigger threshold retuning: ",
    paste(sample_high$project_sample_id,collapse=", "),"."
  ) else "No sample-level doublet-rate review flag was triggered.",
  "No secondary UMI/nFeature doublet cutoff was used. No predicted doublet was rescued because of a preliminary cell-type label.",
  "",
  "## Hold point",
  "Task 005 remains paused. Review Task 004d QC and broad annotation before integration."
)
writeLines(lines,report)
cat("Final Task004d report written: ",report,"\n",sep="")
