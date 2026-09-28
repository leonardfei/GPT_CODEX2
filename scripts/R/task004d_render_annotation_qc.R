#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[[1]],mustWork=TRUE) else normalizePath(".",mustWork=TRUE)
cache <- file.path(root,"objects","task004e_sketch")
ct_file <- file.path(root,"results","task004d_cluster_annotation.csv")
fig_dir <- file.path(root,"figures","task004e")
fig_file <- file.path(fig_dir,"task004d_broad_annotation_qc.pdf")
tmp_fig <- paste0(fig_file,".rendering.pdf")
review_copy <- file.path(root,"figures","task004d_broad_annotation_qc.pdf")
if(!file.exists(ct_file)) stop("Missing corrected cluster annotation table: ",ct_file)
ct <- fread(ct_file)
if(!all(c("dataset","cluster","final_label")%in%names(ct))) stop("Unexpected cluster annotation schema")
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)

render_pdf <- function(){
  pdf(tmp_fig,width=12,height=9,onefile=TRUE)
  on.exit(dev.off(),add=TRUE)
  for(dataset_name in unique(ct$dataset)){
    path <- file.path(cache,paste0(dataset_name,"_broad_reference.rds"))
    if(!file.exists(path)) stop("Missing saved Seurat cache: ",path)
    sk <- readRDS(path)
    clusters <- as.character(Idents(sk))
    label_map <- setNames(ct[dataset==dataset_name,final_label],ct[dataset==dataset_name,cluster])
    labels <- unname(label_map[clusters])
    if(anyNA(labels)) stop("Some saved cache clusters lack final labels: ",dataset_name)
    sk$broad_celltype_cluster <- labels

    print(DimPlot(sk,reduction="umap",group.by="seurat_clusters",label=TRUE,repel=TRUE)+
          ggtitle(paste(dataset_name,"clusters")))
    print(DimPlot(sk,reduction="umap",group.by="broad_celltype_cluster",label=TRUE,repel=TRUE)+
          ggtitle(paste(dataset_name,"broad cell type v2")))
    print(DimPlot(sk,reduction="umap",group.by="tissue")+
          ggtitle(paste(dataset_name,"tissue")))
    print(DimPlot(sk,reduction="pca",group.by="broad_celltype_cluster")+
          ggtitle(paste(dataset_name,"PCA broad cell type v2")))
    print(DimPlot(sk,reduction="pca",group.by="tissue")+
          ggtitle(paste(dataset_name,"PCA tissue")))

    if(dataset_name=="nature_xue" && "source_author_annotation"%in%colnames(sk[[]])){
      print(DimPlot(sk,reduction="umap",group.by="source_author_annotation",label=FALSE)+
            ggtitle("nature_xue source author annotation")+
            guides(color=guide_legend(ncol=5,byrow=TRUE))+
            theme(legend.position="bottom",legend.text=element_text(size=5)))
      print(DimPlot(sk,reduction="pca",group.by="source_author_annotation",label=FALSE)+
            ggtitle("nature_xue PCA source author annotation")+
            guides(color=guide_legend(ncol=5,byrow=TRUE))+
            theme(legend.position="bottom",legend.text=element_text(size=5)))
    }
    rm(sk); gc(verbose=FALSE)
    message("Rendered QC panels for ",dataset_name)
  }
}
render_pdf()
if(!file.exists(tmp_fig) || file.info(tmp_fig)$size<1e6) stop("Temporary review PDF is missing or implausibly small")
if(!file.copy(tmp_fig,fig_file,overwrite=TRUE)) stop("Could not promote temporary review PDF")
if(!file.copy(fig_file,review_copy,overwrite=TRUE)) stop("Could not create required root review copy")
file.remove(tmp_fig)
message("Updated broad annotation QC PDF: ",review_copy)
