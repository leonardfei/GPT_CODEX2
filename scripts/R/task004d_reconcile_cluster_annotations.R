#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleR)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[[1]],mustWork=TRUE) else normalizePath(".",mustWork=TRUE)
res <- file.path(root,"results")
cache <- file.path(root,"objects","task004e_sketch")
final_qs <- file.path(root,"objects","merge","HCC_TA_8datasets_singlets_broad_v1.qs")
tmp_qs <- paste0(final_qs,".marker-fix.tmp.qs")
backup_qs <- paste0(final_qs,".pre_marker_hit_fix.qs")
ann_file <- file.path(res,"task004d_cluster_annotation.csv")
markers_file <- file.path(res,"task004d_cluster_markers.csv.gz")
map_file <- file.path(root,"config","task004d_xue_author_to_broad.tsv")
fig_dir <- file.path(root,"figures","task004e")
fig_file <- file.path(fig_dir,"task004d_broad_annotation_qc.pdf")
review_fig <- file.path(root,"figures","task004d_broad_annotation_qc.pdf")

needed <- c(final_qs,ann_file,markers_file,map_file,
            file.path(cache,"nature_xue_broad_reference.rds"))
if(any(!file.exists(needed))) stop("Missing required cached inputs: ",paste(needed[!file.exists(needed)],collapse=", "))
if(file.exists(tmp_qs)) stop("Temporary output already exists; inspect before continuing: ",tmp_qs)
if(file.exists(backup_qs)) stop("Recovery backup already exists; inspect before continuing: ",backup_qs)

neutrophil_core <- c("FCGR3B","CSF3R","FPR1","FCAR","ELANE","MPO")
cd3_core <- c("CD3D","CD3E","TRAC")
ct <- fread(ann_file)
de <- fread(cmd=paste("gzip -dc",shQuote(markers_file)))
if(!all(c("dataset","cluster","gene")%in%names(de))) stop("Marker file lacks dataset/cluster/gene")
if(!all(c("dataset","cluster","marker_label","marker_margin","marker_hits","anchor_label","anchor_purity")%in%names(ct)))
  stop("Cluster annotation table has an unexpected schema")

# Preserve the FindAllMarkers table's original row order: the original method
# defines support from the first 50 positive markers in each dataset/cluster.
de[,marker_key:=paste(dataset,as.character(cluster),sep="\r")]
top_by_cluster <- split(de$gene,de$marker_key)
ct[,marker_key:=paste(dataset,as.character(cluster),sep="\r")]
for(i in seq_len(nrow(ct))){
  genes <- unique(head(top_by_cluster[[ct$marker_key[[i]]]],50L))
  ct$neutrophil_core_hits[[i]] <- sum(neutrophil_core %in% genes)
  ct$cd3_core_hits[[i]] <- sum(cd3_core %in% genes)
}
ct[,marker_strong:=marker_hits>=2 & marker_margin>=0.15]
ct[marker_label=="Neutrophil",marker_strong:=marker_strong & neutrophil_core_hits>=1]
ct[marker_label=="NK_cell",marker_strong:=marker_strong & cd3_core_hits==0]

ct[,`:=`(broad_reference_label=NA_character_,broad_reference_score=NA_real_,
         final_label="Uncertain/Mixed",final_confidence="uncertain",
         final_basis="conflict_or_insufficient")]
ct[dataset=="nature_xue",`:=`(singler_label=NA_character_,singler_pruned=NA_character_,singler_delta=NA_real_)]

# Replay the original Xue anchor + marker rules without changing thresholds.
for(i in which(ct$dataset=="nature_xue")){
  anchor_label <- ct$anchor_label[[i]]
  anchor_purity <- ct$anchor_purity[[i]]
  anchor_ok <- !is.na(anchor_label) && is.finite(anchor_purity) && anchor_purity>=0.70 &&
    anchor_label==ct$marker_label[[i]] && ct$marker_hits[[i]]>=1
  if(ct$marker_label[[i]]=="Neutrophil")
    anchor_ok <- anchor_ok && ct$neutrophil_core_hits[[i]]>=1
  ct$broad_reference_label[[i]] <- anchor_label
  ct$broad_reference_score[[i]] <- anchor_purity
  if(anchor_ok){
    ct$final_label[[i]] <- anchor_label
    ct$final_confidence[[i]] <- if(anchor_purity>=0.85 && ct$marker_strong[[i]]) "high" else "medium"
    ct$final_basis[[i]] <- "xue_author_consensus_plus_markers"
  } else if(isTRUE(ct$marker_strong[[i]]) && (is.na(anchor_label) || anchor_purity<0.70)){
    ct$final_label[[i]] <- ct$marker_label[[i]]
    ct$final_confidence[[i]] <- if(ct$marker_hits[[i]]>=3 && ct$marker_margin[[i]]>=0.40) "high" else "medium"
    ct$final_basis[[i]] <- "canonical_cluster_markers"
  }
}

# Rebuild the Xue pseudobulk reference from its saved normalized sketch and
# the corrected cluster labels; then replay SingleR for each other dataset.
xue_path <- file.path(cache,"nature_xue_broad_reference.rds")
xue <- readRDS(xue_path)
xue_clusters <- as.character(Idents(xue))
xue_map <- setNames(ct[dataset=="nature_xue",final_label],ct[dataset=="nature_xue",cluster])
xue$broad_celltype_cluster <- unname(xue_map[xue_clusters])
xue_norm <- LayerData(xue[["RNA"]],layer="data")
xue_broad <- as.character(xue$broad_celltype_cluster)
xue_valid <- !is.na(xue_broad) & xue_broad!="Uncertain/Mixed"
xue_group <- paste(as.character(xue$project_sample_id),xue_broad,sep="||")
xue_groups <- unique(xue_group[xue_valid])
if(length(xue_groups)<2L) stop("Corrected Xue reference has fewer than two sample-by-label groups")
ref <- do.call(cbind,lapply(xue_groups,function(g)
  Matrix::rowMeans(xue_norm[,xue_valid & xue_group==g,drop=FALSE])))
rownames(ref) <- rownames(xue_norm)
colnames(ref) <- xue_groups
ref_labels <- sub("^.*\\|\\|","",xue_groups)
rm(xue,xue_norm,xue_broad,xue_valid,xue_group); gc(verbose=FALSE)

datasets <- setdiff(unique(ct$dataset),"nature_xue")
for(dataset_name in datasets){
  path <- file.path(cache,paste0(dataset_name,"_broad_reference.rds"))
  if(!file.exists(path)) stop("Missing saved Seurat cache: ",path)
  sk <- readRDS(path)
  cl <- as.character(Idents(sk))
  norm <- LayerData(sk[["RNA"]],layer="data")
  clusters <- sort(unique(cl))
  test_avg <- do.call(cbind,lapply(clusters,function(k)
    Matrix::rowMeans(norm[,cl==k,drop=FALSE])))
  if(is.null(dim(test_avg))) test_avg <- matrix(test_avg,ncol=1L)
  rownames(test_avg) <- rownames(norm)
  colnames(test_avg) <- clusters
  common <- intersect(rownames(test_avg),rownames(ref))
  pred <- SingleR(test=test_avg[common,,drop=FALSE],ref=ref[common,,drop=FALSE],
                  labels=ref_labels,prune=TRUE)
  ix <- which(ct$dataset==dataset_name)
  match_cluster <- match(as.character(ct$cluster[ix]),rownames(pred))
  if(anyNA(match_cluster)) stop("SingleR did not return all cluster predictions for ",dataset_name)
  ct$singler_label[ix] <- as.character(pred$labels[match_cluster])
  ct$singler_pruned[ix] <- as.character(pred$pruned.labels[match_cluster])
  ct$singler_delta[ix] <- as.numeric(pred$delta.next[match_cluster])
  for(i in ix){
    pruned <- ct$singler_pruned[[i]]
    ref_use <- if(!is.na(pruned) && nzchar(pruned)) pruned else ct$singler_label[[i]]
    ct$broad_reference_label[[i]] <- ref_use
    ct$broad_reference_score[[i]] <- ct$singler_delta[[i]]
    ref_agree <- !is.na(ref_use) && nzchar(ref_use) && ref_use==ct$marker_label[[i]]
    compatible_ref <- ref_agree && ct$marker_hits[[i]]>=1 && ct$marker_margin[[i]]>=0.05
    if(ct$marker_label[[i]]=="Neutrophil")
      compatible_ref <- compatible_ref && ct$neutrophil_core_hits[[i]]>=1
    if(isTRUE(ct$marker_strong[[i]]) && (is.na(ref_use) || !nzchar(ref_use) || ref_agree)){
      ct$final_label[[i]] <- ct$marker_label[[i]]
      ct$final_confidence[[i]] <- if(ref_agree) "high" else "medium"
      ct$final_basis[[i]] <- if(ref_agree) "xue_reference_plus_canonical_markers" else "canonical_cluster_markers"
    } else if(compatible_ref){
      ct$final_label[[i]] <- ref_use
      ct$final_confidence[[i]] <- "medium"
      ct$final_basis[[i]] <- "xue_reference_compatible_markers"
    }
  }
  rm(sk,norm,test_avg,pred); gc(verbose=FALSE)
  message("Reconciled SingleR reference predictions: ",dataset_name)
}

ct[,marker_key:=NULL]
ct <- ct[order(dataset,as.numeric(cluster))]
fwrite(ct,ann_file)
fwrite(ct[final_label=="Uncertain/Mixed"],file.path(res,"task004d_uncertain_clusters.csv"))
fwrite(ct[marker_label=="Neutrophil" | final_label=="Neutrophil" |
            (!is.na(anchor_label)&anchor_label=="Neutrophil") |
            (!is.na(broad_reference_label)&broad_reference_label=="Neutrophil")],
       file.path(res,"task004d_neutrophil_marker_coherence.csv"))

message("Loading final singlet object for label transfer")
obj <- qs::qread(final_qs,use_alt_rep=FALSE,nthreads=8L)
before_ids <- colnames(obj)
before_dims <- dim(obj)
md <- obj[[]]
cell_key <- paste(as.character(md$dataset),as.character(md$broad_cluster_id),sep="\r")
ct_key <- paste(ct$dataset,as.character(ct$cluster),sep="\r")
mt <- match(cell_key,ct_key)
if(anyNA(mt)) stop("Some final-object cells do not map to a corrected dataset/cluster row")
valid_labels <- c("Hepatocyte/Epithelial","T_cell","NK_cell","B_cell","Plasma_cell",
                  "Monocyte/Macrophage","Neutrophil","Dendritic_cell","Mast_cell",
                  "Endothelial","Fibroblast/Mesenchymal","Uncertain/Mixed")
labels <- ct$final_label[mt]
if(anyNA(labels) || any(!labels%in%valid_labels)) stop("Invalid labels after cluster-to-cell mapping")
obj$broad_celltype_v2 <- unname(labels)
obj$broad_annotation_confidence <- unname(ct$final_confidence[mt])
obj$broad_annotation_basis <- unname(ct$final_basis[mt])
obj$broad_reference_label <- unname(ct$broad_reference_label[mt])
obj$broad_reference_score <- as.numeric(ct$broad_reference_score[mt])
obj$annotation_status <- "task004d_broad_annotation_v2_marker_hit_reconciled"

# Recompute the review tables from the corrected object metadata.
md <- obj[[]]
counts <- data.table(dataset=as.character(md$dataset),tissue=as.character(md$tissue),
                     broad_celltype=as.character(md$broad_celltype_v2))[
  ,.(n_cells=.N),by=.(dataset,tissue,broad_celltype)]
counts[,dataset_tissue_total:=sum(n_cells),by=.(dataset,tissue)]
counts[,fraction_within_dataset_tissue:=n_cells/pmax(dataset_tissue_total,1)]
counts[,dataset_total:=sum(n_cells),by=dataset]
counts[,fraction_within_dataset:=n_cells/pmax(dataset_total,1)]
fwrite(counts,file.path(res,"task004d_broad_celltype_counts.csv"))

old_prelim <- if("task004_preliminary_label"%in%names(md)) as.character(md$task004_preliminary_label) else rep(NA_character_,nrow(md))
cross <- data.table(dataset=as.character(md$dataset),task004_preliminary_label=old_prelim,
                    broad_celltype_v2=as.character(md$broad_celltype_v2))[
  ,.(n_cells=.N),by=.(dataset,task004_preliminary_label,broad_celltype_v2)]
fwrite(cross,file.path(res,"task004d_task004_vs_v2.csv"))

xue_map_dt <- fread(map_file)
xue_lookup <- setNames(xue_map_dt$broad_celltype_v2,xue_map_dt$source_author_annotation)
xue_idx <- which(as.character(md$dataset)=="nature_xue")
xue_author <- as.character(md$source_author_annotation[xue_idx])
xue_author_broad <- unname(xue_lookup[xue_author])
xue_author_broad[is.na(xue_author_broad)] <- "Uncertain/Mixed"
xue_cross <- data.table(source_author_annotation=xue_author,source_author_broad=xue_author_broad,
                        broad_celltype_v2=as.character(md$broad_celltype_v2[xue_idx]))[
  ,.(n_cells=.N),by=.(source_author_annotation,source_author_broad,broad_celltype_v2)]
xue_cross[,author_label_total:=sum(n_cells),by=source_author_annotation]
xue_cross[,fraction_within_author_label:=n_cells/pmax(author_label_total,1)]
fwrite(xue_cross,file.path(res,"task004d_xue_author_vs_v2.csv"))

n_neut <- sum(as.character(md$broad_celltype_v2)=="Neutrophil")
n_uncertain <- sum(as.character(md$broad_celltype_v2)=="Uncertain/Mixed")
xue_match <- xue_cross[source_author_broad!="Uncertain/Mixed",
  sum(n_cells[source_author_broad==broad_celltype_v2])/sum(n_cells)]
if(!length(xue_match) || !is.finite(xue_match)) xue_match <- NA_real_
sketch <- fread(file.path(res,"task004d_sketch_usage.csv"))
val <- data.table(status="VALIDATED",biological_review_flag=if(n_neut==0) "no_final_neutrophil_labels" else "none",
                  n_cells=ncol(obj),n_features=nrow(obj),n_datasets=uniqueN(md$dataset),
                  n_broad_types=uniqueN(md$broad_celltype_v2),n_neutrophils=n_neut,
                  n_uncertain=n_uncertain,uncertain_fraction=n_uncertain/ncol(obj),
                  xue_author_broad_concordance=xue_match,
                  any_doublet_remaining=any(as.character(md$scDblFinder.class)=="doublet"),
                  sketch_datasets=paste(sketch[used_sketch==TRUE,dataset],collapse=";"))
if(val$n_cells!=1375127L || val$n_features!=38025L || val$n_datasets!=8L || val$any_doublet_remaining)
  stop("Corrected-object structural validation failed")
if(nrow(val)!=1L) stop("Validation summary malformed")
fwrite(val,file.path(res,"task004d_broad_annotation_validation.csv"))

# Rebuild the multipage broad-label review figure from persisted UMAP/PCA caches.
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
render_review_figure <- function(){
  pdf(fig_file,width=12,height=9,onefile=TRUE)
  on.exit(dev.off(),add=TRUE)
  for(dataset_name in unique(ct$dataset)){
    sk <- readRDS(file.path(cache,paste0(dataset_name,"_broad_reference.rds")))
    sk_cluster <- as.character(Idents(sk))
    sk_map <- setNames(ct[dataset==dataset_name,final_label],ct[dataset==dataset_name,cluster])
    sk$broad_celltype_cluster <- unname(sk_map[sk_cluster])
    print(DimPlot(sk,reduction="umap",group.by="seurat_clusters",label=TRUE,repel=TRUE)+ggtitle(paste(dataset_name,"clusters")))
    print(DimPlot(sk,reduction="umap",group.by="broad_celltype_cluster",label=TRUE,repel=TRUE)+ggtitle(paste(dataset_name,"broad cell type v2")))
    print(DimPlot(sk,reduction="umap",group.by="tissue")+ggtitle(paste(dataset_name,"tissue")))
    print(DimPlot(sk,reduction="pca",group.by="broad_celltype_cluster")+ggtitle(paste(dataset_name,"PCA broad cell type v2")))
    print(DimPlot(sk,reduction="pca",group.by="tissue")+ggtitle(paste(dataset_name,"PCA tissue")))
    if(dataset_name=="nature_xue" && "source_author_annotation"%in%colnames(sk[[]])){
      print(DimPlot(sk,reduction="umap",group.by="source_author_annotation",label=FALSE)+ggtitle("nature_xue source author annotation")+
            guides(color=guide_legend(ncol=5,byrow=TRUE))+
            theme(legend.position="bottom",legend.text=element_text(size=5)))
      print(DimPlot(sk,reduction="pca",group.by="source_author_annotation",label=FALSE)+ggtitle("nature_xue PCA source author annotation")+
            guides(color=guide_legend(ncol=5,byrow=TRUE))+
            theme(legend.position="bottom",legend.text=element_text(size=5)))
    }
    rm(sk); gc(verbose=FALSE)
  }
}
render_review_figure()
if(!file.copy(fig_file,review_fig,overwrite=TRUE)) stop("Could not create required root review figure")

# Write and validate the new QS without touching the old version until complete.
qs::qsave(obj,tmp_qs,preset="high",check_hash=TRUE,nthreads=8L)
rm(obj); gc(verbose=FALSE)
check <- qs::qread(tmp_qs,use_alt_rep=FALSE,nthreads=8L)
if(dim(check)[[1]]!=before_dims[[1]] || dim(check)[[2]]!=before_dims[[2]] ||
   !identical(colnames(check),before_ids)) stop("Corrected QS dimensions/cell order validation failed")
if(any(is.na(check$broad_celltype_v2)) || any(as.character(check$scDblFinder.class)=="doublet"))
  stop("Corrected QS metadata validation failed")
if(!file.rename(final_qs,backup_qs)) stop("Could not preserve previous final QS as recovery backup")
if(!file.rename(tmp_qs,final_qs)){
  file.rename(backup_qs,final_qs)
  stop("Could not promote corrected QS; restored prior final QS")
}
message("Corrected object promoted. Prior object preserved at: ",backup_qs)
message("Final broad labels: ",paste(names(sort(table(check$broad_celltype_v2),decreasing=TRUE)),collapse=", "))
message("Neutrophil cells: ",n_neut,"; uncertain cells: ",n_uncertain)
message("Corrected review figure: ",review_fig)
