#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(SingleCellExperiment); library(scDblFinder); library(BiocParallel); library(ggplot2)
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
sdf_version <- packageVersion("scDblFinder")
if(requireNamespace("xgboost",quietly=TRUE) &&
   packageVersion("xgboost") >= package_version("3.1.0") &&
   sdf_version < package_version("1.24.8")){
  stop("scDblFinder >=1.24.8 is required with xgboost >=3.1; found scDblFinder ",
       as.character(sdf_version)," and xgboost ",as.character(packageVersion("xgboost")))
}
sdf_has_dbr_per1k <- "dbr.per1k" %in% names(formals(scDblFinder::scDblFinder))
scdblfinder_rate_args <- function(n_cells) {
  if (sdf_has_dbr_per1k) {
    list(dbr=NULL,dbr.per1k=0.008,dbr.sd=NULL)
  } else {
    # Legacy scDblFinder defines dbr as a fraction for the current capture.
    # Scale the approved 0.8%-per-1,000-cells rate explicitly for this sample.
    list(dbr=min(1,0.008*n_cells/1000),dbr.sd=NULL)
  }
}
dbr_api_mode <- if(sdf_has_dbr_per1k) "native_dbr.per1k" else
  "legacy_dbr_explicit_per_1000_equivalent"
dir.create(dirname(out_qs),recursive=TRUE,showWarnings=FALSE)
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(dirname(report_path),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(root,"figures"),recursive=TRUE,showWarnings=FALSE)

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
if(!inherits(obj,"Seurat")) stop("Invalid Seurat object")
assay_names <- SeuratObject::Assays(obj)
if(!("RNA"%in%assay_names)) stop("Seurat object lacks the RNA assay")
if(!identical(SeuratObject::Layers(obj[["RNA"]]),"counts"))
  stop("Expected one RNA counts layer")
counts_layer <- obj[["RNA"]]@layers[["counts"]]
if(!inherits(counts_layer,"task004b_chunked_counts")) stop("Expected task004b_chunked_counts")
md <- as.data.frame(obj[[]])
if(ncol(obj)!=1490852L) stop("Preflight cell count mismatch: ",ncol(obj))
if(nrow(obj)!=38025L) stop("Preflight feature count mismatch: ",nrow(obj))
if(uniqueN(md$dataset)!=8L) stop("Preflight dataset count mismatch")
if("project_patient_id"%in%colnames(md) && uniqueN(md$project_patient_id)!=132L)
  stop("Preflight project_patient_id count mismatch")
if(!"project_sample_id"%in%colnames(md)){
  if(!all(c("dataset","sample_id")%in%colnames(md))) stop("Missing sample identifiers")
  md$project_sample_id <- paste(md$dataset,md$sample_id,sep="::")
}
req <- c("dataset","sample_id","patient_id","project_patient_id","tissue",
         "project_sample_id","percent.mt","project_broad_celltype")
if(length(setdiff(req,colnames(md)))) stop("Missing required metadata")
capture_ids <- as.character(md$project_sample_id)
capture_levels <- unique(capture_ids)
if(length(capture_levels)!=194L) stop("Preflight project_sample_id count mismatch: ",length(capture_levels))

calls_file <- file.path(results_dir,"task004d_scdblfinder_cell_calls.csv.gz")
sample_file <- file.path(results_dir,"task004d_scdblfinder_by_sample.csv")
reuse_call_table <- tolower(arg_value("reuse-call-table","false")) %in% c("true","1","yes")
if(reuse_call_table){
  if(!file.exists(calls_file)||!file.exists(sample_file))
    stop("Requested call-table reuse, but the all-cell calls or sample summary file is missing")
  calls <- fread(cmd=paste("gzip -dc --",shQuote(calls_file)))
  samples <- fread(sample_file)
  required_call_cols <- c("cell_id","global_index","dataset","sample_id","patient_id",
    "tissue","project_sample_id","nCount_RNA","nFeature_RNA","percent.mt",
    "preliminary_task004_broad","scDblFinder.score","scDblFinder.class","scDblFinder.status")
  if(length(setdiff(required_call_cols,names(calls))))
    stop("Saved all-cell call table is missing required fields")
  if(nrow(calls)!=ncol(obj)) stop("Saved call table cell count does not match the approved input")
  setorder(calls,global_index)
  if(!identical(calls$global_index,seq_len(ncol(obj)))||
     !identical(as.character(calls$cell_id),colnames(obj))||
     !identical(as.character(calls$dataset),as.character(md$dataset))||
     !identical(as.character(calls$sample_id),as.character(md$sample_id))||
     !identical(as.character(calls$project_sample_id),capture_ids))
    stop("Saved call table identifiers/order do not match the approved input metadata")
  if(anyNA(calls$scDblFinder.score)||anyNA(calls$scDblFinder.class)||
     anyNA(calls$scDblFinder.status)||any(calls$scDblFinder.status!="SCORED"))
    stop("Saved call table contains unscored or incomplete cells")
  required_sample_cols <- c("project_sample_id","n_cells","expected_dbr","dbr_api_mode","status")
  if(length(setdiff(required_sample_cols,names(samples)))||
     anyDuplicated(samples$project_sample_id)||
     !setequal(samples$project_sample_id,capture_levels)||
     anyNA(samples$status)||any(samples$status!="SCORED")||
     anyNA(samples$dbr_api_mode)||any(samples$dbr_api_mode!=dbr_api_mode))
    stop("Saved sample summary is incomplete or uses a different rate API")
  observed_sizes <- calls[,.(observed_n=.N),by=project_sample_id]
  obs_i <- match(samples$project_sample_id,observed_sizes$project_sample_id)
  if(anyNA(obs_i)||any(samples$n_cells!=observed_sizes$observed_n[obs_i])||
     any(abs(samples$expected_dbr-pmin(1,0.008*samples$n_cells/1000))>1e-12))
    stop("Saved sample summary does not agree with all-cell scores or approved rate")
  old_score_cols <- intersect(c("score_median","score_q25","score_q75",
    "singlet_nCount_median","doublet_nCount_median","singlet_nFeature_median",
    "doublet_nFeature_median","review_flag"),names(samples))
  if(length(old_score_cols)) samples[,(old_score_cols):=NULL]
  message("Reusing validated all-cell calls for ",nrow(calls)," cells across ",
          uniqueN(calls$project_sample_id)," captures; no doublet models are rerun")
} else {
calls_list <- vector("list",length(capture_levels))
sample_list <- vector("list",length(capture_levels))
for(s in seq_along(capture_levels)){
  cap <- capture_levels[[s]]; idx <- which(capture_ids==cap); n <- length(idx)
  message(sprintf("[%d/%d] %s: %d cells",s,length(capture_levels),cap,n))
  count_col <- if("hgnc_nCount_RNA"%in%colnames(md)) "hgnc_nCount_RNA" else "nCount_RNA"
  feature_col <- if("hgnc_nFeature_RNA"%in%colnames(md)) "hgnc_nFeature_RNA" else "nFeature_RNA"
  base <- data.table(cell_id=colnames(obj)[idx],global_index=idx,
    dataset=as.character(md$dataset[idx]),sample_id=as.character(md$sample_id[idx]),
    patient_id=as.character(md$patient_id[idx]),tissue=as.character(md$tissue[idx]),
    project_sample_id=cap,
    nCount_RNA=as.numeric(md[[count_col]][idx]),
    nFeature_RNA=as.numeric(md[[feature_col]][idx]))
  base[,project_patient_id:=as.character(md$project_patient_id[idx])]
  base[,percent.mt:=as.numeric(md$percent.mt[idx])]
  base[,preliminary_task004_broad:=as.character(md$project_broad_celltype[idx])]
  if("source_author_annotation"%in%colnames(md))
    base[,source_author_annotation:=as.character(md$source_author_annotation[idx])]
  x <- counts_layer[,idx,drop=FALSE]
  if(!inherits(x,"dgCMatrix")) x <- as(x,"dgCMatrix")
  if(any(Matrix::colSums(x)<=0)) stop("Zero-count cell in ",cap)
  x <- x[Matrix::rowSums(x>0)>0,,drop=FALSE]
  sce <- SingleCellExperiment(assays=list(counts=x))
  colnames(sce) <- base$cell_id
  set.seed(seed_base+s)
  cluster_mode <- if(n >= 500L) TRUE else NULL
  sce <- do.call(scDblFinder::scDblFinder,c(
    list(sce,clusters=cluster_mode),scdblfinder_rate_args(n),
    list(iter=2,BPPARAM=SerialParam(),verbose=FALSE)
  ))
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
    dbr_api_mode=dbr_api_mode,
    cluster_mode=ifelse(n>=500L,"cluster_based","random"),status="SCORED")
  calls_list[[s]] <- base
  rm(x,sce,cd,base); gc(verbose=FALSE)
}
calls <- rbindlist(calls_list,fill=TRUE); samples <- rbindlist(sample_list,fill=TRUE)
}
setorder(calls,global_index)
if(nrow(calls)!=ncol(obj)||!identical(calls$cell_id,colnames(obj))) stop("Call table mismatch")
if(anyNA(calls$scDblFinder.score)||anyNA(calls$scDblFinder.class))
  stop("At least one original cell lacks a scDblFinder score or class")
prelim <- as.character(md$project_broad_celltype)
calls[,candidate_neutrophil:=tolower(trimws(preliminary_task004_broad))=="neutrophil"]
calls[,keep_after_scdblfinder:=scDblFinder.class!="doublet"|is.na(scDblFinder.class)]
keep <- calls$keep_after_scdblfinder
calls[,removal_reason:=fifelse(scDblFinder.class=="doublet","scDblFinder_doublet",
                               fifelse(keep_after_scdblfinder,"retained_singlet","other_removed"))]
broad_audit <- calls[,.(n_cells=.N,n_doublet=sum(scDblFinder.class=="doublet",na.rm=TRUE),
  doublet_fraction=mean(scDblFinder.class=="doublet",na.rm=TRUE)),
  by=.(dataset,tissue,preliminary_task004_broad)]
fwrite(calls,file.path(results_dir,"task004d_scdblfinder_cell_calls.csv.gz"),compress="gzip")
fwrite(samples,file.path(results_dir,"task004d_scdblfinder_by_sample.csv"))
fwrite(broad_audit,file.path(results_dir,"task004d_scdblfinder_by_preliminary_broad.csv"))
score_stats <- calls[, .(
  score_median=median(scDblFinder.score,na.rm=TRUE),
  score_q25=quantile(scDblFinder.score,0.25,na.rm=TRUE,names=FALSE),
  score_q75=quantile(scDblFinder.score,0.75,na.rm=TRUE,names=FALSE),
  singlet_nCount_median=median(nCount_RNA[scDblFinder.class=="singlet"],na.rm=TRUE),
  doublet_nCount_median=median(nCount_RNA[scDblFinder.class=="doublet"],na.rm=TRUE),
  singlet_nFeature_median=median(nFeature_RNA[scDblFinder.class=="singlet"],na.rm=TRUE),
  doublet_nFeature_median=median(nFeature_RNA[scDblFinder.class=="doublet"],na.rm=TRUE)
),by=project_sample_id]
samples <- merge(samples,score_stats,by="project_sample_id",all.x=TRUE,sort=FALSE)
samples[,review_flag:=fifelse(doublet_fraction>0.30,"REVIEW_HIGH_GT30PCT",
                      fifelse(doublet_fraction<0.001 & n_cells>=1000,"REVIEW_LOW_LT0.1PCT","none"))]
fwrite(samples,file.path(results_dir,"task004d_scdblfinder_by_sample.csv"))

summarize_call_class <- function(d,group_cols,level_name){
  z <- d[,.(n_cells=.N,
    score_q25=quantile(scDblFinder.score,0.25,names=FALSE),
    score_median=median(scDblFinder.score),
    score_q75=quantile(scDblFinder.score,0.75,names=FALSE),
    nCount_q25=quantile(nCount_RNA,0.25,names=FALSE),
    nCount_median=median(nCount_RNA),
    nCount_q75=quantile(nCount_RNA,0.75,names=FALSE),
    nFeature_q25=quantile(nFeature_RNA,0.25,names=FALSE),
    nFeature_median=median(nFeature_RNA),
    nFeature_q75=quantile(nFeature_RNA,0.75,names=FALSE),
    percent_mt_q25=quantile(percent.mt,0.25,names=FALSE),
    percent_mt_median=median(percent.mt),
    percent_mt_q75=quantile(percent.mt,0.75,names=FALSE)),
    by=c(group_cols,"scDblFinder.class")]
  z[,analysis_level:=level_name]
  z
}
call_qc <- rbindlist(list(
  summarize_call_class(calls,character(),"overall"),
  summarize_call_class(calls,"dataset","dataset"),
  summarize_call_class(calls,c("dataset","tissue"),"dataset_tissue"),
  summarize_call_class(calls,"project_sample_id","sample")
),fill=TRUE)
fwrite(call_qc,file.path(results_dir,"task004d_scdblfinder_cell_qc_distributions.csv"))
overall <- calls[,.(n_cells_before=.N,
  n_singlet=sum(scDblFinder.class=="singlet"),
  n_doublet=sum(scDblFinder.class=="doublet"),
  doublet_fraction=mean(scDblFinder.class=="doublet"),
  score_q25=quantile(scDblFinder.score,0.25,names=FALSE),
  score_median=median(scDblFinder.score),
  score_q75=quantile(scDblFinder.score,0.75,names=FALSE))]
fwrite(overall,file.path(results_dir,"task004d_scdblfinder_overall.csv"))

n_before <- length(keep); n_after <- sum(keep); n_removed <- n_before-n_after
neut_before <- sum(tolower(prelim)=="neutrophil",na.rm=TRUE)
neut_removed <- sum(tolower(prelim)=="neutrophil"&!keep,na.rm=TRUE)
neut_ret <- if(neut_before) (neut_before-neut_removed)/neut_before else NA_real_

neut_audit <- calls[,.(n_cells_before=.N,
  n_candidate_neutrophils_before=sum(candidate_neutrophil),
  n_candidate_neutrophils_after=sum(candidate_neutrophil&keep_after_scdblfinder),
  n_candidate_neutrophils_removed=sum(candidate_neutrophil&!keep_after_scdblfinder),
  candidate_neutrophil_retention_fraction=ifelse(sum(candidate_neutrophil)>0,
    sum(candidate_neutrophil&keep_after_scdblfinder)/sum(candidate_neutrophil),NA_real_),
  n_candidate_neutrophils_removed_doublet=sum(candidate_neutrophil&scDblFinder.class=="doublet"),
  n_other_cells_before=sum(!candidate_neutrophil),
  n_other_cells_after=sum(!candidate_neutrophil&keep_after_scdblfinder),
  n_other_cells_removed=sum(!candidate_neutrophil&!keep_after_scdblfinder),
  n_removed_doublet=sum(scDblFinder.class=="doublet"),
  n_unscored=sum(is.na(scDblFinder.class))),by=.(dataset,tissue,project_sample_id)]
fwrite(neut_audit,file.path(results_dir,"task004d_neutrophil_doublet_retention.csv"))

distribution_phase <- function(d,phase){
  z <- copy(d); z[,audit_phase:=phase]
  z[,candidate_group:=fifelse(candidate_neutrophil,"candidate_neutrophil","other_cells")]
  z[,.(n_cells=.N,
    nCount_q25=quantile(nCount_RNA,0.25,names=FALSE),nCount_median=median(nCount_RNA),
    nCount_q75=quantile(nCount_RNA,0.75,names=FALSE),
    nFeature_q25=quantile(nFeature_RNA,0.25,names=FALSE),nFeature_median=median(nFeature_RNA),
    nFeature_q75=quantile(nFeature_RNA,0.75,names=FALSE),
    percent_mt_q25=quantile(percent.mt,0.25,names=FALSE),percent_mt_median=median(percent.mt),
    percent_mt_q75=quantile(percent.mt,0.75,names=FALSE)),
    by=.(dataset,tissue,project_sample_id,audit_phase,candidate_group)]
}
neut_dist <- rbindlist(list(distribution_phase(calls,"before_filter"),
                            distribution_phase(calls[which(calls$keep_after_scdblfinder)],"after_filter")),
                       use.names=TRUE)
fwrite(neut_dist,file.path(results_dir,"task004d_neutrophil_retention_distributions.csv"))

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
  dbr=if(sdf_has_dbr_per1k) "automatic" else "min(1, 0.008*n_cells/1000)",
  dbr_per1k=0.008,dbr_api_mode=dbr_api_mode,dbr_sd=NULL,iter=2,
  cluster_rule="clusters=TRUE for n>=500; clusters=NULL for n<500",
  seed_base=seed_base,n_cells_before=n_before,
  n_doublets_removed=n_removed,n_cells_after=n_after,
  preliminary_neutrophils_before=neut_before,
  preliminary_neutrophils_removed=neut_removed,
  preliminary_neutrophil_retention=neut_ret,
  scDblFinder_rate_api=dbr_api_mode,
  calls_reused_from_validated_table=reuse_call_table
)
ds <- calls[,.(n_cells_before=.N,n_doublet=sum(scDblFinder.class=="doublet",na.rm=TRUE),
  n_singlet=sum(scDblFinder.class=="singlet",na.rm=TRUE),
  n_retained=sum(keep_after_scdblfinder),
  doublet_fraction=mean(scDblFinder.class=="doublet",na.rm=TRUE),
  score_median=median(scDblFinder.score,na.rm=TRUE),
  score_q25=quantile(scDblFinder.score,0.25,na.rm=TRUE,names=FALSE),
  score_q75=quantile(scDblFinder.score,0.75,na.rm=TRUE,names=FALSE)),
  by=dataset]
fwrite(ds,file.path(results_dir,"task004d_scdblfinder_by_dataset.csv"))
tissue_ds <- calls[,.(n_cells_before=.N,
  n_doublet=sum(scDblFinder.class=="doublet",na.rm=TRUE),
  n_singlet=sum(scDblFinder.class=="singlet",na.rm=TRUE),
  doublet_fraction=mean(scDblFinder.class=="doublet",na.rm=TRUE)),
  by=.(dataset,tissue)]
fwrite(tissue_ds,file.path(results_dir,"task004d_scdblfinder_by_dataset_tissue.csv"))

pdf(file.path(root,"figures","task004d_scdblfinder_qc.pdf"),width=12,height=8)
print(ggplot(samples,aes(x=reorder(project_sample_id,doublet_fraction),y=doublet_fraction,fill=dataset))+
  geom_col()+coord_flip()+theme_classic()+labs(x="project_sample_id",y="Predicted doublet fraction"))
print(ggplot(calls,aes(x=scDblFinder.class,y=log10(nCount_RNA+1),fill=scDblFinder.class))+
  geom_boxplot(outlier.shape=NA)+facet_wrap(~dataset,scales="free_y")+theme_classic()+
  labs(x=NULL,y="log10(nCount_RNA+1)"))
print(ggplot(calls,aes(x=scDblFinder.class,y=log10(nFeature_RNA+1),fill=scDblFinder.class))+
  geom_boxplot(outlier.shape=NA)+facet_wrap(~dataset,scales="free_y")+theme_classic()+
  labs(x=NULL,y="log10(nFeature_RNA+1)"))
dev.off()

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
  scDblFinder_version=as.character(sdf_version),
  scDblFinder_rate_api=dbr_api_mode,dbr_per1k=0.008,
  scDblFinder_calls_reused=reuse_call_table)
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
  if(sdf_has_dbr_per1k) "Expected rate parameter: dbr.per1k=0.008." else
    "Expected rate parameter: legacy explicit dbr=min(1, 0.008*n_cells/1000), equivalent to 0.8% per 1,000 cells.",
  "Preliminary Task 004 labels were used only for retention auditing."),
  report_path)
message("Task 004d COMPLETED")
