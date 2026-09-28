#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(SingleCellExperiment); library(scDblFinder); library(BiocParallel)
})
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name, default=NULL) {
  k <- paste0("--",name); i <- which(args==k)
  if(!length(i) || i==length(args)) return(default)
  args[[i+1L]]
}
root <- normalizePath(arg_value("project-root","."), mustWork=TRUE)
source_qs <- arg_value("source-qs",file.path(root,"objects","merge","HCC_TA_8datasets_merged_review_v1.qs"))
out_qs <- arg_value("out-qs",file.path(root,"objects","merge","HCC_TA_8datasets_singlets_v1.qs"))
results_dir <- arg_value("results-dir",file.path(root,"results"))
report_path <- arg_value("report",file.path(root,"reports","task_004d_scdblfinder_report.md"))
seed_base <- as.integer(arg_value("seed","44000"))
dir.create(dirname(out_qs),recursive=TRUE,showWarnings=FALSE)
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(dirname(report_path),recursive=TRUE,showWarnings=FALSE)

dim.task004b_chunked_counts <- function(x) c(length(x$feature_names),length(x$cell_names))
dimnames.task004b_chunked_counts <- function(x) list(x$feature_names,x$cell_names)
assign("[.task004b_chunked_counts", function(x,i,j,drop=FALSE){
  d <- dim(x)
  if(missing(i)) i <- seq_len(d[[1L]])
  if(missing(j)) j <- seq_len(d[[2L]])
  if(is.character(i)) i <- match(i,x$feature_names)
  if(is.character(j)) j <- match(j,x$cell_names)
  if(anyNA(i)||anyNA(j)) stop("Unknown feature/cell in chunked subset")
  if(!length(i)||!length(j)) return(Matrix::Matrix(0,nrow=length(i),ncol=length(j),sparse=TRUE,
    dimnames=list(x$feature_names[i],x$cell_names[j])))
  pieces <- vector("list",length(x$chunks)); use <- logical(length(x$chunks))
  for(k in seq_along(x$chunks)){
    hit <- which(j>x$offsets[[k]] & j<=x$offsets[[k+1L]])
    if(!length(hit)) next
    jj <- j[hit]-x$offsets[[k]]
    pieces[[k]] <- x$chunks[[k]][i,jj,drop=FALSE]
    colnames(pieces[[k]]) <- x$cell_names[j[hit]]
    use[[k]] <- TRUE
  }
  z <- do.call(cbind,pieces[use]); rownames(z) <- x$feature_names[i]; z
}, envir=.GlobalEnv)

new_chunked <- function(chunks,features,cells) structure(
  list(chunks=chunks,feature_names=features,cell_names=cells,
       offsets=c(0L,cumsum(vapply(chunks,ncol,integer(1L))))),
  class="task004b_chunked_counts"
)

if(!file.exists(source_qs)) stop("Source QS not found: ",source_qs)
message("Reading ",source_qs)
obj <- qs::qread(source_qs,use_alt_rep=FALSE,nthreads=8L)
if(!inherits(obj,"Seurat")||!"RNA"%in%Assays(obj)) stop("Invalid Seurat object")
if(!identical(Layers(obj[["RNA"]]),"counts")) stop("Expected one RNA counts layer")
counts_layer <- obj[["RNA"]]@layers[["counts"]]
if(!inherits(counts_layer,"task004b_chunked_counts")) stop("Expected task004b_chunked_counts")
md <- obj[[]]
if(!"project_sample_id"%in%colnames(md)){
  if(!all(c("dataset","sample_id")%in%colnames(md))) stop("Missing sample identifiers")
  md$project_sample_id <- paste(md$dataset,md$sample_id,sep="::")
}
req <- c("dataset","sample_id","patient_id","tissue","project_sample_id")
if(length(setdiff(req,colnames(md)))) stop("Missing required metadata")
capture_ids <- as.character(md$project_sample_id)
capture_levels <- unique(capture_ids)

calls_list <- vector("list",length(capture_levels))
sample_list <- vector("list",length(capture_levels))
for(s in seq_along(capture_levels)){
  cap <- capture_levels[[s]]; idx <- which(capture_ids==cap); n <- length(idx)
  message(sprintf("[%d/%d] %s: %d cells",s,length(capture_levels),cap,n))
  base <- data.table(cell_id=colnames(obj)[idx],global_index=idx,
    dataset=as.character(md$dataset[idx]),sample_id=as.character(md$sample_id[idx]),
    patient_id=as.character(md$patient_id[idx]),tissue=as.character(md$tissue[idx]),
    project_sample_id=cap)
  x <- counts_layer[,idx,drop=FALSE]
  if(!inherits(x,"dgCMatrix")) x <- as(x,"dgCMatrix")
  if(any(Matrix::colSums(x)<=0)) stop("Zero-count cell in ",cap)
  x <- x[Matrix::rowSums(x>0)>0,,drop=FALSE]
  sce <- SingleCellExperiment(assays=list(counts=x))
  colnames(sce) <- base$cell_id
  set.seed(seed_base+s)
  cluster_mode <- if(n >= 500L) TRUE else NULL
  sce <- scDblFinder(sce,clusters=cluster_mode,dbr=NULL,dbr.per1k=0.008,
                     dbr.sd=NULL,iter=2,BPPARAM=SerialParam(),verbose=FALSE)
  cd <- as.data.frame(colData(sce))
  if(!all(c("scDblFinder.score","scDblFinder.class")%in%colnames(cd)))
    stop("Missing scDblFinder outputs for ",cap)
  base[,scDblFinder.score:=as.numeric(cd$scDblFinder.score)]
  base[,scDblFinder.class:=as.character(cd$scDblFinder.class)]
  base[,scDblFinder.status:="SCORED"]
  extras <- setdiff(grep("^scDblFinder\\.",colnames(cd),value=TRUE),
                    c("scDblFinder.score","scDblFinder.class"))
  for(nm in extras){
    v <- cd[[nm]]; if(is.factor(v)) v <- as.character(v); base[[nm]] <- v
  }
  nd <- sum(base$scDblFinder.class=="doublet",na.rm=TRUE)
  sample_list[[s]] <- data.table(project_sample_id=cap,dataset=base$dataset[[1]],
    sample_id=base$sample_id[[1]],tissue=base$tissue[[1]],n_cells=n,n_doublet=nd,
    n_singlet=n-nd,doublet_fraction=nd/n,expected_dbr=min(1,0.008*n/1000),
    cluster_mode=ifelse(n>=500L,"cluster_based","random"),status="SCORED")
  calls_list[[s]] <- base
  rm(x,sce,cd,base); gc(verbose=FALSE)
}
calls <- rbindlist(calls_list,fill=TRUE); samples <- rbindlist(sample_list,fill=TRUE)
setorder(calls,global_index)
if(nrow(calls)!=ncol(obj)||!identical(calls$cell_id,colnames(obj))) stop("Call table mismatch")
calls[,keep_after_scdblfinder:=scDblFinder.class!="doublet"|is.na(scDblFinder.class)]
keep <- calls$keep_after_scdblfinder
prelim <- if("project_broad_celltype"%in%colnames(md)) as.character(md$project_broad_celltype) else rep(NA_character_,nrow(md))
calls[,preliminary_task004_broad:=prelim]
broad_audit <- calls[,.(n_cells=.N,n_doublet=sum(scDblFinder.class=="doublet",na.rm=TRUE),
  doublet_fraction=mean(scDblFinder.class=="doublet",na.rm=TRUE)),
  by=.(dataset,tissue,preliminary_task004_broad)]
fwrite(calls,file.path(results_dir,"task004d_scdblfinder_cell_calls.csv.gz"),compress="gzip")
fwrite(samples,file.path(results_dir,"task004d_scdblfinder_by_sample.csv"))
fwrite(broad_audit,file.path(results_dir,"task004d_scdblfinder_by_preliminary_broad.csv"))
high <- samples[status=="SCORED"&doublet_fraction>0.30]
if(nrow(high)) stop("Filtering halted: >30% doublets in ",paste(high$project_sample_id,collapse=", "))

n_before <- length(keep); n_after <- sum(keep); n_removed <- n_before-n_after
neut_before <- sum(tolower(prelim)=="neutrophil",na.rm=TRUE)
neut_removed <- sum(tolower(prelim)=="neutrophil"&!keep,na.rm=TRUE)
neut_ret <- if(neut_before) (neut_before-neut_removed)/neut_before else NA_real_

new_chunks <- vector("list",length(counts_layer$chunks)); new_cells <- character()
for(k in seq_along(counts_layer$chunks)){
  a <- counts_layer$offsets[[k]]+1L; b <- counts_layer$offsets[[k+1L]]
  z <- counts_layer$chunks[[k]][,keep[a:b],drop=FALSE]
  new_chunks[[k]] <- z; new_cells <- c(new_cells,colnames(z))
}
if(length(new_cells)!=n_after) stop("Filtered chunk total mismatch")
new_counts <- new_chunked(new_chunks,counts_layer$feature_names,new_cells)
new_md <- md[keep,,drop=FALSE]
new_md$scDblFinder.score <- calls$scDblFinder.score[keep]
new_md$scDblFinder.class <- calls$scDblFinder.class[keep]
new_md$scDblFinder.status <- calls$scDblFinder.status[keep]
new_md$doublet_filter_status <- "retained_singlet"
rownames(new_md) <- new_cells
assay <- obj[["RNA"]]; assay@layers$counts <- new_counts
assay@cells <- SeuratObject:::LogMap(new_cells)
obj@assays$RNA <- assay; obj@meta.data <- new_md
obj@active.ident <- factor(rep("unassigned",n_after),levels="unassigned")
names(obj@active.ident) <- new_cells
obj@misc$task004d_scdblfinder <- list(
  source_qs=source_qs,version=as.character(packageVersion("scDblFinder")),
  detection_unit="project_sample_id",n_samples=length(capture_levels),
  min_cells_for_scoring=min_cells,dbr="automatic",dbr_per1k=0.008,
  nfeatures=1352,dims=20,seed_base=seed_base,n_cells_before=n_before,
  n_doublets_removed=n_removed,n_cells_after=n_after,
  preliminary_neutrophils_before=neut_before,
  preliminary_neutrophils_removed=neut_removed,
  preliminary_neutrophil_retention=neut_ret
)
ds <- calls[,.(n_cells_before=.N,n_doublet=sum(scDblFinder.class=="doublet",na.rm=TRUE),
  n_retained=sum(keep_after_scdblfinder),
  doublet_fraction=mean(scDblFinder.class=="doublet",na.rm=TRUE)),by=dataset]
fwrite(ds,file.path(results_dir,"task004d_scdblfinder_by_dataset.csv"))

message("Writing ",out_qs)
qs::qsave(obj,out_qs,preset="high",check_hash=TRUE,nthreads=8L)
expected_ids <- colnames(obj); expected_features <- nrow(obj)
rm(obj,assay,new_md,new_counts,new_chunks); gc(verbose=FALSE)
chk <- qs::qread(out_qs,use_alt_rep=FALSE,nthreads=8L)
if(ncol(chk)!=n_after||nrow(chk)!=expected_features||!identical(colnames(chk),expected_ids))
  stop("Filtered QS validation failed")
if(any(as.character(chk$scDblFinder.class)=="doublet")) stop("Doublets remain in filtered object")
val <- data.table(status="VALIDATED",source_qs=source_qs,filtered_qs=out_qs,
  n_cells_before=n_before,n_doublets_removed=n_removed,n_cells_after=ncol(chk),
  n_features=nrow(chk),n_samples=length(capture_levels),
  preliminary_neutrophils_before=neut_before,
  preliminary_neutrophils_removed=neut_removed,
  preliminary_neutrophil_retention=neut_ret,
  scDblFinder_version=as.character(packageVersion("scDblFinder")))
fwrite(val,file.path(results_dir,"task004d_scdblfinder_validation.csv"))
writeLines(c("# Task 004d report - scDblFinder",
  "","## Status","COMPLETED","",
  paste0("Cells before: ",format(n_before,big.mark=",")),
  paste0("Doublets removed: ",format(n_removed,big.mark=",")),
  paste0("Cells retained: ",format(n_after,big.mark=",")),
  paste0("Samples/captures: ",length(capture_levels)),
  paste0("Preliminary neutrophil retention: ",
         ifelse(is.na(neut_ret),"NA",sprintf("%.3f%%",100*neut_ret))),
  "","Detection was performed independently by project_sample_id on raw counts.",
  "Preliminary Task 004 labels were used only for retention auditing."),
  report_path)
message("Task 004d COMPLETED")
