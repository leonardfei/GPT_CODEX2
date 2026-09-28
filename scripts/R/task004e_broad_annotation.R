#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix); library(SingleR)
})
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name, default=NULL){
  k <- paste0("--",name); i <- which(args==k)
  if(!length(i)||i==length(args)) return(default)
  args[[i+1L]]
}
root <- normalizePath(arg_value("project-root","."),mustWork=TRUE)
source_qs <- arg_value("source-qs",file.path(root,"objects","merge","HCC_TA_8datasets_singlets_v1.qs"))
out_qs <- arg_value("out-qs",file.path(root,"objects","merge","HCC_TA_8datasets_singlets_broad_v1.qs"))
shared_file <- arg_value("shared-features",file.path(root,"results","task004c_shared_hgnc_features_8of8.txt"))
results_dir <- arg_value("results-dir",file.path(root,"results"))
sketch_dir <- arg_value("sketch-dir",file.path(root,"objects","task004e_sketch"))
fig_dir <- arg_value("figures-dir",file.path(root,"figures","task004e"))
report_path <- arg_value("report",file.path(root,"reports","task_004d_report.md"))
xue_map_file <- arg_value("xue-map",file.path(root,"config","task004d_xue_author_to_broad.tsv"))
sketch_n <- as.integer(arg_value("sketch-cells","50000"))
block_n <- as.integer(arg_value("projection-block","2000"))
resolution <- as.numeric(arg_value("resolution","0.6"))
seed <- as.integer(arg_value("seed","40500"))
dir.create(dirname(out_qs),recursive=TRUE,showWarnings=FALSE)
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(sketch_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(dirname(report_path),recursive=TRUE,showWarnings=FALSE)

dim.task004b_chunked_counts <- function(x) c(length(x$feature_names),length(x$cell_names))
dimnames.task004b_chunked_counts <- function(x) list(x$feature_names,x$cell_names)
assign("[.task004b_chunked_counts",function(x,i,j,drop=FALSE){
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
},envir=.GlobalEnv)

marker_sets <- list(
  "Hepatocyte/Epithelial"=c("ALB","APOA1","APOA2","TTR","ASGR1","KRT8","KRT18","EPCAM","KRT19","KRT7"),
  "T_cell"=c("CD3D","CD3E","TRAC","CD2","LTB"),
  "NK_cell"=c("NKG7","GNLY","KLRD1","PRF1","CTSW"),
  "B_cell"=c("CD79A","CD79B","MS4A1","CD37","CD74","CD22","CD19"),
  "Plasma_cell"=c("MZB1","JCHAIN","XBP1","DERL3","SDC1","IGKC"),
  "Monocyte/Macrophage"=c("LST1","TYROBP","FCER1G","CTSS","LILRB1","C1QC","APOC1","SPP1"),
  "Neutrophil"=c("FCGR3B","CSF3R","CXCR2","FPR1","FCAR","S100A8","S100A9","NAMPT","MCEMP1","ELANE","MPO"),
  "Dendritic_cell"=c("FCER1A","CD1C","CLEC10A","CLEC9A","XCR1","GZMB","TCF4"),
  "Mast_cell"=c("TPSAB1","TPSB2","KIT","MS4A2","HDC","CPA3"),
  "Endothelial"=c("PECAM1","VWF","EMCN","KDR","RAMP2","ENG","ESM1"),
  "Fibroblast/Mesenchymal"=c("COL1A1","COL1A2","COL3A1","DCN","LUM","COL6A1","COL6A2","PDGFRA","FAP","THY1")
)
neutrophil_core <- c("FCGR3B","CSF3R","FPR1","FCAR","ELANE","MPO")

if(!file.exists(xue_map_file)) stop("Xue author crosswalk missing: ",xue_map_file)
xue_map <- fread(xue_map_file)
if(!all(c("source_author_annotation","broad_celltype_v2")%in%colnames(xue_map)))
  stop("Invalid Xue author crosswalk")
xue_lookup <- setNames(xue_map$broad_celltype_v2,xue_map$source_author_annotation)
map_author_broad <- function(x){
  y <- unname(xue_lookup[as.character(x)])
  y[is.na(y)] <- "Uncertain/Mixed"
  y
}

stratified_sample <- function(indices,strata,target,seed){
  if(length(indices)<=target) return(indices)
  set.seed(seed); st <- as.character(strata); u <- unique(st)
  per <- max(1L,ceiling(target/length(u)))
  take <- unlist(lapply(u,function(g){
    z <- indices[st==g]
    if(length(z)<=per) z else sample(z,per)
  }),use.names=FALSE)
  take <- unique(take)
  if(length(take)>target) take <- sample(take,target)
  if(length(take)<target){
    rem <- setdiff(indices,take)
    take <- c(take,sample(rem,min(target-length(take),length(rem))))
  }
  sort(take)
}

cluster_marker_annotation <- function(sketch,markers_de,author_broad=NULL){
  cl <- as.character(Idents(sketch)); clusters <- sort(unique(cl))
  present <- intersect(unique(unlist(marker_sets)),rownames(sketch))
  dat <- LayerData(sketch[["RNA"]],layer="data")[present,,drop=FALSE]
  avg <- sapply(clusters,function(k) Matrix::rowMeans(dat[,cl==k,drop=FALSE]))
  if(is.null(dim(avg))) avg <- matrix(avg,ncol=1,dimnames=list(present,clusters))
  rownames(avg) <- present; colnames(avg) <- clusters
  z <- t(scale(t(avg))); z[!is.finite(z)] <- 0
  score <- matrix(NA_real_,nrow=length(clusters),ncol=length(marker_sets),
                  dimnames=list(clusters,names(marker_sets)))
  hits <- matrix(0L,nrow=length(clusters),ncol=length(marker_sets),
                 dimnames=list(clusters,names(marker_sets)))
  top_by_cluster <- split(markers_de$gene,as.character(markers_de$cluster))
  for(k in clusters){
    topgenes <- unique(head(top_by_cluster[[k]],50))
    for(lb in names(marker_sets)){
      g <- intersect(marker_sets[[lb]],rownames(z))
      score[k,lb] <- if(length(g)) mean(z[g,k]) else -Inf
      hits[k,lb] <- sum(marker_sets[[lb]]%in%topgenes)
    }
  }
  out <- rbindlist(lapply(clusters,function(k){
    v <- score[k,]; ord <- order(v,decreasing=TRUE)
    top <- names(v)[ord[[1]]]; second <- if(length(ord)>1) v[ord[[2]]] else -Inf
    data.table(cluster=k,marker_label=top,marker_score=v[ord[[1]]],
               marker_margin=v[ord[[1]]]-second,marker_hits=hits[k,top])
  }))
  if(!is.null(author_broad)){
    anchor <- data.table(cluster=cl,author_broad=author_broad)[!is.na(author_broad),
      .(anchor_label=names(sort(table(author_broad),decreasing=TRUE))[1],
        anchor_purity=max(table(author_broad))/.N,n_anchor=.N),by=cluster]
    out <- merge(out,anchor,by="cluster",all.x=TRUE)
  } else {
    out[,c("anchor_label","anchor_purity","n_anchor"):=list(NA_character_,NA_real_,0L)]
  }
  out
}

normalize_sparse <- function(x,lib){
  sf <- 1e4/pmax(as.numeric(lib),1)
  y <- x%*%Diagonal(x=sf); y <- as(y,"dgCMatrix"); y@x <- log1p(y@x); y
}

project_to_clusters <- function(counts_layer,global_idx,shared_genes,sketch,hvg,cluster_table,block_n){
  hdat <- LayerData(sketch[["RNA"]],layer="data")[hvg,,drop=FALSE]
  mu <- Matrix::rowMeans(hdat); n <- ncol(hdat)
  ss <- Matrix::rowSums(hdat*hdat)
  sdv <- sqrt(pmax((ss-n*mu^2)/pmax(n-1,1),1e-8))
  load <- Loadings(sketch[["pca"]])[hvg,,drop=FALSE]
  dims <- min(30,ncol(load)); load <- load[,seq_len(dims),drop=FALSE]
  emb <- Embeddings(sketch,"pca")[,seq_len(dims),drop=FALSE]
  scl <- as.character(Idents(sketch)); clusters <- sort(unique(scl))
  cent <- do.call(rbind,lapply(clusters,function(k) colMeans(emb[scl==k,,drop=FALSE])))
  rownames(cent) <- clusters
  final_map <- setNames(cluster_table$final_label,cluster_table$cluster)
  conf_map <- setNames(cluster_table$final_confidence,cluster_table$cluster)
  ans_cluster <- character(length(global_idx)); ans_label <- character(length(global_idx))
  ans_conf <- character(length(global_idx)); ans_margin <- numeric(length(global_idx))
  for(a in seq(1,length(global_idx),by=block_n)){
    b <- min(length(global_idx),a+block_n-1L); gi <- global_idx[a:b]
    xs <- counts_layer[shared_genes,gi,drop=FALSE]; lib <- Matrix::colSums(xs)
    xh <- xs[hvg,,drop=FALSE]; nh <- normalize_sparse(xh,lib)
    dense <- as.matrix(nh); dense <- sweep(dense,1,mu,"-"); dense <- sweep(dense,1,sdv,"/")
    dense[dense>10] <- 10; dense[dense< -10] <- -10
    pcs <- t(dense)%*%load
    q2 <- rowSums(pcs^2); c2 <- rowSums(cent^2)
    d2 <- outer(q2,c2,"+")-2*(pcs%*%t(cent))
    best <- max.col(-d2,ties.method="first")
    bestd <- d2[cbind(seq_len(nrow(d2)),best)]
    d2b <- d2; d2b[cbind(seq_len(nrow(d2b)),best)] <- Inf
    second <- apply(d2b,1,min)
    cc <- rownames(cent)[best]
    ans_cluster[a:b] <- cc; ans_label[a:b] <- final_map[cc]; ans_conf[a:b] <- conf_map[cc]
    ans_margin[a:b] <- (second-bestd)/pmax(second,1e-8)
    rm(xs,xh,nh,dense,pcs,d2,d2b); gc(verbose=FALSE)
  }
  data.table(global_index=global_idx,reference_cluster=ans_cluster,
             projected_label=ans_label,cluster_confidence=ans_conf,
             projection_margin=ans_margin)
}

if(!file.exists(source_qs)) stop("Filtered source QS not found: ",source_qs)
if(!file.exists(shared_file)) stop("Shared HGNC feature list not found: ",shared_file)
shared_genes <- scan(shared_file,what=character(),quiet=TRUE)
message("Reading filtered object: ",source_qs)
obj <- qs::qread(source_qs,use_alt_rep=FALSE,nthreads=8L)
counts_layer <- obj[["RNA"]]@layers[["counts"]]
if(!inherits(counts_layer,"task004b_chunked_counts")) stop("Expected chunked counts")
md <- obj[[]]
if(any(as.character(md$scDblFinder.class)=="doublet")) stop("Filtered object still contains doublets")
shared_genes <- intersect(shared_genes,counts_layer$feature_names)
if(length(shared_genes)<10000) stop("Too few shared genes: ",length(shared_genes))
datasets <- unique(as.character(md$dataset))
if(!"nature_xue"%in%datasets) stop("nature_xue reference dataset missing")

new_label <- rep(NA_character_,nrow(md)); new_conf <- rep(NA_character_,nrow(md))
new_cluster <- rep(NA_character_,nrow(md)); proj_margin <- rep(NA_real_,nrow(md))
cluster_rows <- list(); count_rows <- list(); de_rows <- list()

run_dataset <- function(dataset_name,is_reference=FALSE,reference_matrix=NULL,reference_labels=NULL){
  idx <- which(as.character(md$dataset)==dataset_name)
  strata <- as.character(md$project_sample_id[idx])
  author_all <- if("source_author_annotation"%in%colnames(md)) map_author_broad(md$source_author_annotation[idx]) else rep(NA_character_,length(idx))
  if(is_reference&&any(!is.na(author_all))) strata <- paste(strata,ifelse(is.na(author_all),"unmapped",author_all),sep="||")
  sk_idx <- stratified_sample(idx,strata,min(sketch_n,length(idx)),seed+match(dataset_name,datasets))
  x <- counts_layer[shared_genes,sk_idx,drop=FALSE]
  sk <- CreateSeuratObject(counts=x,assay="RNA",project=dataset_name,meta.data=md[sk_idx,,drop=FALSE])
  sk <- NormalizeData(sk,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  sk <- FindVariableFeatures(sk,selection.method="vst",nfeatures=min(3000,nrow(sk)),verbose=FALSE)
  hvg <- VariableFeatures(sk)
  sk <- ScaleData(sk,features=hvg,verbose=FALSE)
  npcs <- min(50,length(hvg)-1L)
  sk <- RunPCA(sk,features=hvg,npcs=npcs,verbose=FALSE)
  nd <- min(30,npcs)
  sk <- FindNeighbors(sk,reduction="pca",dims=seq_len(nd),verbose=FALSE)
  sk <- FindClusters(sk,resolution=resolution,random.seed=seed,verbose=FALSE)
  sk <- RunUMAP(sk,reduction="pca",dims=seq_len(nd),seed.use=seed,verbose=FALSE)
  de <- FindAllMarkers(sk,only.pos=TRUE,min.pct=0.15,logfc.threshold=0.25,
                       max.cells.per.ident=3000,random.seed=seed,verbose=FALSE)
  if(nrow(de)){de$dataset <- dataset_name; de_rows[[dataset_name]] <<- as.data.table(de)}
  else de <- data.frame(gene=character(),cluster=character())

  author_sk <- if(is_reference&&"source_author_annotation"%in%colnames(sk[[]])) map_author_broad(sk$source_author_annotation) else NULL
  ct <- cluster_marker_annotation(sk,de,author_sk)

  if(is_reference){
    ct[,c("singler_label","singler_pruned","singler_delta"):=list(NA_character_,NA_character_,NA_real_)]
    ct[,final_label:=fifelse(!is.na(anchor_label)&anchor_purity>=0.65,anchor_label,
                      fifelse(marker_hits>=2&marker_margin>=0.15,marker_label,"Other/uncertain"))]
    ct[,final_confidence:=fifelse(!is.na(anchor_label)&anchor_purity>=0.80,"high_author_anchor",
                           fifelse(!is.na(anchor_label)&anchor_purity>=0.65,"medium_author_anchor",
                           fifelse(marker_hits>=3&marker_margin>=0.40,"high_marker",
                           fifelse(marker_hits>=2&marker_margin>=0.15,"medium_marker","low_uncertain"))))]
  } else {
    cl <- as.character(Idents(sk)); norm <- LayerData(sk[["RNA"]],layer="data")
    cls <- sort(unique(cl))
    test_avg <- sapply(cls,function(k) Matrix::rowMeans(norm[,cl==k,drop=FALSE]))
    if(is.null(dim(test_avg))) test_avg <- matrix(test_avg,ncol=1,dimnames=list(rownames(norm),cls))
    common <- intersect(rownames(test_avg),rownames(reference_matrix))
    pred <- SingleR(test=test_avg[common,,drop=FALSE],ref=reference_matrix[common,,drop=FALSE],
                    labels=reference_labels,prune=TRUE)
    pd <- data.table(cluster=rownames(pred),singler_label=as.character(pred$labels),
                     singler_pruned=as.character(pred$pruned.labels),
                     singler_delta=as.numeric(pred$delta.next))
    ct <- merge(ct,pd,by="cluster",all.x=TRUE)
    ct[,ref_use:=fifelse(!is.na(singler_pruned)&nzchar(singler_pruned),singler_pruned,singler_label)]
    ct[,agree:=!is.na(ref_use)&ref_use==marker_label]
    ct[,final_label:=fifelse(agree,marker_label,
                      fifelse((is.na(ref_use)|!nzchar(ref_use))&marker_hits>=2&marker_margin>=0.15,marker_label,
                      fifelse(marker_hits<2|marker_margin<0.15,ref_use,
                      fifelse(marker_hits>=3&marker_margin>=0.60,marker_label,"Other/uncertain"))))]
    ct[is.na(final_label)|!nzchar(final_label),final_label:="Other/uncertain"]
    ct[,final_confidence:=fifelse(agree&marker_hits>=2,"high_reference_marker_agreement",
                           fifelse(final_label=="Other/uncertain","low_conflict",
                           fifelse(final_label==marker_label,"medium_marker","medium_reference")))]
    ct[,c("ref_use","agree"):=NULL]
  }

  ct[,dataset:=dataset_name]; cluster_rows[[dataset_name]] <<- ct
  map_final <- setNames(ct$final_label,ct$cluster)
  sk$broad_celltype_cluster <- map_final[as.character(Idents(sk))]
  saveRDS(sk,file.path(sketch_dir,paste0(dataset_name,"_broad_sketch.rds")),compress=FALSE)

  if(is_reference){
    norm <- LayerData(sk[["RNA"]],layer="data"); broad <- as.character(sk$broad_celltype_cluster)
    valid <- broad!="Other/uncertain"&!is.na(broad)
    grp <- paste(as.character(sk$project_sample_id),broad,sep="||")
    groups <- unique(grp[valid])
    ref <- sapply(groups,function(g) Matrix::rowMeans(norm[,valid&grp==g,drop=FALSE]))
    if(is.null(dim(ref))) ref <- matrix(ref,ncol=1,dimnames=list(rownames(norm),groups))
    reference_matrix <- ref; reference_labels <- sub("^.*\\|\\|","",groups)
  }

  proj <- project_to_clusters(counts_layer,idx,shared_genes,sk,hvg,ct,block_n)
  if(is_reference&&any(!is.na(author_all))){
    z <- !is.na(author_all); proj$projected_label[z] <- author_all[z]
    proj$cluster_confidence[z] <- "author_anchor_cell"
  }
  new_label[idx] <<- proj$projected_label; new_conf[idx] <<- proj$cluster_confidence
  new_cluster[idx] <<- proj$reference_cluster; proj_margin[idx] <<- proj$projection_margin
  count_rows[[dataset_name]] <<- data.table(dataset=dataset_name,broad_celltype=proj$projected_label)[,
    .(n_cells=.N),by=.(dataset,broad_celltype)]

  pdf(file.path(fig_dir,paste0(dataset_name,"_sketch_umap.pdf")),width=9,height=7)
  print(DimPlot(sk,reduction="umap",group.by="broad_celltype_cluster",label=TRUE,repel=TRUE)+
        ggtitle(paste(dataset_name,"broad-cell sketch")))
  dev.off()
  list(reference_matrix=reference_matrix,reference_labels=reference_labels)
}

message("Building nature_xue reference")
ref <- run_dataset("nature_xue",TRUE)
xue_ref <- ref$reference_matrix; xue_ref_labels <- ref$reference_labels
if(is.null(xue_ref)||ncol(xue_ref)<2) stop("Failed to build nature_xue reference")
for(d in setdiff(datasets,"nature_xue")){
  message("Annotating ",d); run_dataset(d,FALSE,xue_ref,xue_ref_labels); gc(verbose=FALSE)
}

if(any(is.na(new_label))) stop("Some cells lack broad annotation")
old_prelim <- if("project_broad_celltype"%in%colnames(md)) as.character(md$project_broad_celltype) else rep(NA_character_,nrow(md))
obj$task004_preliminary_broad_celltype <- old_prelim
obj$project_broad_celltype_v2 <- new_label
obj$broad_annotation_confidence_v2 <- new_conf
obj$broad_reference_cluster_v2 <- new_cluster
obj$broad_projection_margin_v2 <- proj_margin
obj$annotation_status_v2 <- "task004e_cluster_reference_annotation"
obj@misc$task004e_broad_annotation <- list(
  source_qs=source_qs,shared_feature_file=shared_file,n_shared_genes=length(shared_genes),
  sketch_cells_per_dataset=sketch_n,projection_block=block_n,resolution=resolution,
  reference_dataset="nature_xue",reference_source="source_author_annotation",
  method="dataset-stratified sketch clustering + cluster markers + nature_xue reference transfer + PCA-centroid projection",
  seed=seed)

all_clusters <- rbindlist(cluster_rows,fill=TRUE); all_counts <- rbindlist(count_rows,fill=TRUE)
fwrite(all_clusters,file.path(results_dir,"task004e_broad_cluster_annotations.csv"))
fwrite(all_counts,file.path(results_dir,"task004e_broad_cell_counts.csv"))
if(length(de_rows)) fwrite(rbindlist(de_rows,fill=TRUE),file.path(results_dir,"task004e_sketch_cluster_markers.csv.gz"),compress="gzip")
cross <- data.table(dataset=as.character(md$dataset),old=old_prelim,new=new_label)[,
  .(n_cells=.N),by=.(dataset,old,new)]
fwrite(cross,file.path(results_dir,"task004e_old_vs_new_broad_crosswalk.csv"))

message("Writing ",out_qs)
qs::qsave(obj,out_qs,preset="high",check_hash=TRUE,nthreads=8L)
expected_n <- ncol(obj); expected_ids <- colnames(obj)
rm(obj); gc(verbose=FALSE)
chk <- qs::qread(out_qs,use_alt_rep=FALSE,nthreads=8L)
if(ncol(chk)!=expected_n||!identical(colnames(chk),expected_ids)) stop("Annotated QS validation failed")
if(any(is.na(chk$project_broad_celltype_v2))) stop("Missing v2 broad labels")
val <- data.table(status="VALIDATED",source_qs=source_qs,annotated_qs=out_qs,
  n_cells=ncol(chk),n_features=nrow(chk),n_datasets=uniqueN(chk$dataset),
  n_broad_types=uniqueN(chk$project_broad_celltype_v2),
  n_uncertain=sum(chk$project_broad_celltype_v2=="Other/uncertain"),
  uncertain_fraction=mean(chk$project_broad_celltype_v2=="Other/uncertain"))
fwrite(val,file.path(results_dir,"task004e_broad_annotation_validation.csv"))
writeLines(c("# Task 004e report - broad cell annotation","",
  "## Status","COMPLETED","",
  paste0("Cells annotated: ",format(ncol(chk),big.mark=",")),
  paste0("Broad labels: ",uniqueN(chk$project_broad_celltype_v2)),
  paste0("Other/uncertain: ",format(sum(chk$project_broad_celltype_v2=="Other/uncertain"),big.mark=","),
         " (",sprintf("%.2f%%",100*mean(chk$project_broad_celltype_v2=="Other/uncertain")),")"),
  "","Annotation used dataset-specific sketch PCA/clustering, cluster markers, nature_xue author-label reference transfer, and PCA-centroid projection to all retained singlets.",
  "Previous Task 004 labels are preserved separately and were not used as the reference."),
  report_path)
message("Task 004e COMPLETED")
