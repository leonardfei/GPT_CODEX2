#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix); library(SingleR); library(ggplot2)
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

available_memory_bytes <- function(){
  if(!file.exists("/proc/meminfo")) return(NA_real_)
  line <- grep("^MemAvailable:",readLines("/proc/meminfo"),value=TRUE)
  if(!length(line)) return(NA_real_)
  as.numeric(sub("^MemAvailable:\\s*([0-9]+)\\s+kB.*$","\\1",line[[1]]))*1024
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
               marker_margin=v[ord[[1]]]-second,marker_hits=hits[k,top],
               neutrophil_core_hits=sum(neutrophil_core%in%topgenes),
               cd3_core_hits=sum(c("CD3D","CD3E","TRAC")%in%topgenes))
  }))
  if(!is.null(author_broad)){
    author_broad[author_broad=="Uncertain/Mixed"] <- NA_character_
    anchor <- data.table(cluster=cl,author_broad=author_broad)[!is.na(author_broad),
      .(anchor_label=names(sort(table(author_broad),decreasing=TRUE))[1],
        anchor_purity=max(table(author_broad))/.N,n_anchor=.N),by=cluster]
    out <- merge(out,anchor,by="cluster",all.x=TRUE)
  } else {
    out[,c("anchor_label","anchor_purity","n_anchor"):=list(NA_character_,NA_real_,0L)]
  }
  out
}

cluster_program_evidence <- function(sketch,dataset_name){
  cl <- as.character(Idents(sketch))
  clusters <- sort(unique(cl))
  norm <- LayerData(sketch[["RNA"]],layer="data")
  cnt <- LayerData(sketch[["RNA"]],layer="counts")
  out <- list(); z <- 1L
  for(k in clusters){
    cells <- which(cl==k)
    for(lb in names(marker_sets)){
      g <- intersect(marker_sets[[lb]],rownames(sketch))
      if(length(g)){
        mean_expr <- mean(Matrix::rowMeans(norm[g,cells,drop=FALSE]))
        det <- cnt[g,cells,drop=FALSE]>0
        mean_det <- mean(Matrix::rowMeans(det))
        ge2 <- mean(Matrix::colSums(det)>=2)
      } else {
        mean_expr <- NA_real_; mean_det <- NA_real_; ge2 <- NA_real_
      }
      out[[z]] <- data.table(dataset=dataset_name,cluster=k,program=lb,
                             n_cells=length(cells),n_program_genes=length(g),
                             mean_normalized_expression=mean_expr,
                             mean_marker_detection_fraction=mean_det,
                             fraction_cells_ge2_markers=ge2)
      z <- z+1L
    }
  }
  rbindlist(out,fill=TRUE)
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
  basis_map <- setNames(cluster_table$final_basis,cluster_table$cluster)
  ref_label_map <- setNames(cluster_table$broad_reference_label,cluster_table$cluster)
  ref_score_map <- setNames(cluster_table$broad_reference_score,cluster_table$cluster)
  ans_cluster <- character(length(global_idx)); ans_label <- character(length(global_idx))
  ans_conf <- character(length(global_idx)); ans_basis <- character(length(global_idx))
  ans_ref_label <- character(length(global_idx)); ans_ref_score <- numeric(length(global_idx))
  ans_margin <- numeric(length(global_idx))
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
    ans_basis[a:b] <- basis_map[cc]; ans_ref_label[a:b] <- ref_label_map[cc]
    ans_ref_score[a:b] <- as.numeric(ref_score_map[cc])
    ans_margin[a:b] <- (second-bestd)/pmax(second,1e-8)
    rm(xs,xh,nh,dense,pcs,d2,d2b); gc(verbose=FALSE)
  }
  data.table(global_index=global_idx,reference_cluster=ans_cluster,
             projected_label=ans_label,cluster_confidence=ans_conf,
             annotation_basis=ans_basis,reference_label=ans_ref_label,
             reference_score=ans_ref_score,projection_margin=ans_margin)
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
new_basis <- rep(NA_character_,nrow(md)); new_ref_label <- rep(NA_character_,nrow(md))
new_ref_score <- rep(NA_real_,nrow(md)); new_cluster <- rep(NA_character_,nrow(md))
cycling_state_global <- rep(NA_character_,nrow(md))
proj_margin <- rep(NA_real_,nrow(md))
cluster_rows <- list(); count_rows <- list(); de_rows <- list(); sketch_flags <- list(); program_rows <- list()

run_dataset <- function(dataset_name,is_reference=FALSE,reference_matrix=NULL,reference_labels=NULL){
  idx <- which(as.character(md$dataset)==dataset_name)
  author_all <- if("source_author_annotation"%in%colnames(md)) map_author_broad(md$source_author_annotation[idx]) else rep(NA_character_,length(idx))
  # The first nature_xue full-cell attempt exceeded 135 GB RSS before the
  # server OOM-killed it. Treat that observed failure as direct evidence that
  # its full-cell graph/marker workflow is infeasible at the available RAM.
  available_bytes <- available_memory_bytes()
  estimate_bytes <- length(idx)*min(3000L,length(shared_genes))*8*10 + length(idx)*32768
  use_sketch <- length(idx)>250000L &&
    (dataset_name=="nature_xue" || !is.finite(available_bytes) ||
       estimate_bytes>0.75*available_bytes)
  message(sprintf("%s: %s analysis (%d cells; conservative full working-set estimate %.1f GiB; MemAvailable %.1f GiB)",
    dataset_name,ifelse(use_sketch,"leverage-score sketch","full-cell"),
    length(idx),estimate_bytes/1024^3,available_bytes/1024^3))
  x <- counts_layer[shared_genes,idx,drop=FALSE]
  sk <- CreateSeuratObject(counts=x,assay="RNA",project=dataset_name,meta.data=md[idx,,drop=FALSE])
  sk <- NormalizeData(sk,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  sk <- FindVariableFeatures(sk,selection.method="vst",nfeatures=min(3000,nrow(sk)),verbose=FALSE)
  hvg <- VariableFeatures(sk)
  hvg <- hvg[!grepl("^MT-|^RPL|^RPS",hvg,ignore.case=TRUE)]
  if(length(hvg)<1000) stop("Too few non-mito/ribosomal HVGs for ",dataset_name)
  if(use_sketch){
    sk <- Seurat::SketchData(
      object=sk,assay="RNA",ncells=min(sketch_n,ncol(sk)),
      sketched.assay="task004d_sketch",method="LeverageScore",
      var.name="task004d_leverage_score",over.write=TRUE,
      seed=seed,verbose=TRUE,features=hvg
    )
    sketch_cells <- SeuratObject::Cells(sk[["task004d_sketch"]])
    if(length(sketch_cells)!=min(sketch_n,length(idx)))
      stop("Seurat leverage-score sketch returned an unexpected cell count for ",dataset_name)
    sketch_global_idx <- match(sketch_cells,colnames(obj))
    if(anyNA(sketch_global_idx)||anyDuplicated(sketch_global_idx))
      stop("Leverage-score sketch cell IDs did not map uniquely to the source object for ",dataset_name)
    rm(sk,x); gc(verbose=FALSE)
    x_sketch <- counts_layer[shared_genes,sketch_global_idx,drop=FALSE]
    sk <- CreateSeuratObject(counts=x_sketch,assay="RNA",project=dataset_name,
                             meta.data=md[sketch_global_idx,,drop=FALSE])
    rm(x_sketch); gc(verbose=FALSE)
    sk <- NormalizeData(sk,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
    VariableFeatures(sk) <- hvg
  } else {
    VariableFeatures(sk) <- hvg
  }
  sk <- ScaleData(sk,features=hvg,verbose=FALSE)
  npcs <- min(30,length(hvg)-1L)
  sk <- RunPCA(sk,features=hvg,npcs=npcs,verbose=FALSE)
  sk <- FindNeighbors(sk,reduction="pca",dims=seq_len(npcs),verbose=FALSE)
  sk <- FindClusters(sk,resolution=c(0.4,0.8),random.seed=seed,verbose=FALSE)
  res08 <- grep("res\\.0\\.8$|res\\.0.8$",colnames(sk[[]]),value=TRUE)
  if(!length(res08)) res08 <- grep("0.8$",colnames(sk[[]]),value=TRUE)
  if(!length(res08)) stop("Could not identify resolution 0.8 cluster column for ",dataset_name)
  Idents(sk) <- sk[[res08[[1]]]][,1]
  sk$broad_cluster_resolution <- "0.8"
  sk <- RunUMAP(sk,reduction="pca",dims=seq_len(npcs),seed.use=seed,verbose=FALSE)
  de <- FindAllMarkers(sk,only.pos=TRUE,min.pct=0.15,logfc.threshold=0.25,
                       max.cells.per.ident=3000,random.seed=seed,verbose=FALSE)
  if(nrow(de)){de$dataset <- dataset_name; de_rows[[dataset_name]] <<- as.data.table(de)}
  else de <- data.frame(gene=character(),cluster=character())

  author_sk <- if(is_reference&&"source_author_annotation"%in%colnames(sk[[]])) {
    map_author_broad(sk$source_author_annotation)
  } else NULL
  ct <- cluster_marker_annotation(sk,de,author_sk)
  pe <- cluster_program_evidence(sk,dataset_name)
  program_rows[[dataset_name]] <<- pe
  prop <- pe[order(cluster,-mean_normalized_expression)]
  top2 <- prop[, .(
    program_top=program[[1]],
    program_top_mean=mean_normalized_expression[[1]],
    program_second=if(.N>=2) program[[2]] else NA_character_,
    program_second_mean=if(.N>=2) mean_normalized_expression[[2]] else NA_real_
  ),by=cluster]
  ct <- merge(ct,top2,by="cluster",all.x=TRUE)
  cycling_genes <- c("MKI67","TOP2A","UBE2C","CENPF","TYMS")
  top_by_cluster <- split(de$gene,as.character(de$cluster))
  ct[,cycling_state:=vapply(cluster,function(k){
    tg <- unique(head(top_by_cluster[[k]],50))
    if(sum(cycling_genes%in%tg)>=2) "cycling" else "noncycling"
  },character(1))]
  ct[,marker_strong:=marker_hits>=2&marker_margin>=0.15]
  ct[marker_label=="Neutrophil",
     marker_strong:=marker_strong&neutrophil_core_hits>=1]
  ct[marker_label=="NK_cell",
     marker_strong:=marker_strong&cd3_core_hits==0]
  ct[,c("broad_reference_label","broad_reference_score","final_label",
        "final_confidence","final_basis") :=
       list(NA_character_,NA_real_,"Uncertain/Mixed","uncertain","conflict_or_insufficient")]

  if(is_reference){
    for(r in seq_len(nrow(ct))){
      anchor_ok <- !is.na(ct$anchor_label[r]) && ct$anchor_purity[r]>=0.70 &&
                   ct$anchor_label[r]==ct$marker_label[r] && ct$marker_hits[r]>=1
      if(ct$marker_label[r]=="Neutrophil")
        anchor_ok <- anchor_ok && ct$neutrophil_core_hits[r]>=1
      ct$broad_reference_label[r] <- ct$anchor_label[r]
      ct$broad_reference_score[r] <- ct$anchor_purity[r]
      if(anchor_ok){
        ct$final_label[r] <- ct$anchor_label[r]
        ct$final_confidence[r] <- ifelse(ct$anchor_purity[r]>=0.85 &&
                                        ct$marker_strong[r],"high","medium")
        ct$final_basis[r] <- "xue_author_consensus_plus_markers"
      } else if(ct$marker_strong[r] &&
                (is.na(ct$anchor_label[r]) || ct$anchor_purity[r]<0.70)){
        ct$final_label[r] <- ct$marker_label[r]
        ct$final_confidence[r] <- ifelse(ct$marker_hits[r]>=3 &&
                                        ct$marker_margin[r]>=0.40,"high","medium")
        ct$final_basis[r] <- "canonical_cluster_markers"
      }
    }
    ct[,c("singler_label","singler_pruned","singler_delta") :=
       list(NA_character_,NA_character_,NA_real_)]
  } else {
    cl <- as.character(Idents(sk)); norm <- LayerData(sk[["RNA"]],layer="data")
    cls <- sort(unique(cl))
    test_avg <- sapply(cls,function(k) Matrix::rowMeans(norm[,cl==k,drop=FALSE]))
    if(is.null(dim(test_avg)))
      test_avg <- matrix(test_avg,ncol=1,dimnames=list(rownames(norm),cls))
    common <- intersect(rownames(test_avg),rownames(reference_matrix))
    pred <- SingleR(test=test_avg[common,,drop=FALSE],
                    ref=reference_matrix[common,,drop=FALSE],
                    labels=reference_labels,prune=TRUE)
    pd <- data.table(cluster=rownames(pred),singler_label=as.character(pred$labels),
                     singler_pruned=as.character(pred$pruned.labels),
                     singler_delta=as.numeric(pred$delta.next))
    ct <- merge(ct,pd,by="cluster",all.x=TRUE)
    for(r in seq_len(nrow(ct))){
      ref_use <- if(!is.na(ct$singler_pruned[r])&&nzchar(ct$singler_pruned[r]))
        ct$singler_pruned[r] else ct$singler_label[r]
      ct$broad_reference_label[r] <- ref_use
      ct$broad_reference_score[r] <- ct$singler_delta[r]
      ref_agree <- !is.na(ref_use)&&nzchar(ref_use)&&ref_use==ct$marker_label[r]
      marker_ok <- isTRUE(ct$marker_strong[r])
      compatible_ref <- ref_agree && ct$marker_hits[r]>=1 && ct$marker_margin[r]>=0.05
      if(ct$marker_label[r]=="Neutrophil")
        compatible_ref <- compatible_ref && ct$neutrophil_core_hits[r]>=1
      if(marker_ok && (is.na(ref_use)||!nzchar(ref_use)||ref_agree)){
        ct$final_label[r] <- ct$marker_label[r]
        ct$final_confidence[r] <- ifelse(ref_agree,"high","medium")
        ct$final_basis[r] <- ifelse(ref_agree,
                                   "xue_reference_plus_canonical_markers",
                                   "canonical_cluster_markers")
      } else if(compatible_ref){
        ct$final_label[r] <- ref_use
        ct$final_confidence[r] <- "medium"
        ct$final_basis[r] <- "xue_reference_compatible_markers"
      }
    }
  }

  ct[,dataset:=dataset_name]
  ct[,used_sketch:=use_sketch]
  cluster_rows[[dataset_name]] <<- ct
  sketch_flags[[dataset_name]] <<- data.table(dataset=dataset_name,n_cells=length(idx),
    n_analysis_cells=ncol(sk),used_sketch=use_sketch,
    estimated_full_working_set_gib=estimate_bytes/1024^3,
    mem_available_before_gib=available_bytes/1024^3,
    sampling_method=if(use_sketch) "Seurat::SketchData method=LeverageScore" else "none_full_cell")
  map_final <- setNames(ct$final_label,ct$cluster)
  map_conf <- setNames(ct$final_confidence,ct$cluster)
  map_basis <- setNames(ct$final_basis,ct$cluster)
  map_ref <- setNames(ct$broad_reference_label,ct$cluster)
  map_ref_score <- setNames(ct$broad_reference_score,ct$cluster)
  map_cycle <- setNames(ct$cycling_state,ct$cluster)
  sk$broad_celltype_cluster <- map_final[as.character(Idents(sk))]
  sk$broad_cycling_state <- map_cycle[as.character(Idents(sk))]
  saveRDS(sk,file.path(sketch_dir,paste0(dataset_name,"_broad_reference.rds")),compress=FALSE)

  if(is_reference){
    norm <- LayerData(sk[["RNA"]],layer="data")
    broad <- as.character(sk$broad_celltype_cluster)
    valid <- broad!="Uncertain/Mixed"&!is.na(broad)
    grp <- paste(as.character(sk$project_sample_id),broad,sep="||")
    groups <- unique(grp[valid])
    ref <- sapply(groups,function(g) Matrix::rowMeans(norm[,valid&grp==g,drop=FALSE]))
    if(is.null(dim(ref))) ref <- matrix(ref,ncol=1,dimnames=list(rownames(norm),groups))
    reference_matrix <- ref
    reference_labels <- sub("^.*\\|\\|","",groups)
  }

  if(use_sketch){
    proj <- project_to_clusters(counts_layer,idx,shared_genes,sk,hvg,ct,block_n)
  } else {
    cc <- as.character(Idents(sk))
    proj <- data.table(global_index=idx,reference_cluster=cc,
                       projected_label=map_final[cc],
                       cluster_confidence=map_conf[cc],
                       annotation_basis=map_basis[cc],
                       reference_label=map_ref[cc],
                       reference_score=as.numeric(map_ref_score[cc]),
                       cycling_state=map_cycle[cc],
                       projection_margin=1)
  }

  new_label[idx] <<- proj$projected_label
  new_conf[idx] <<- proj$cluster_confidence
  new_basis[idx] <<- proj$annotation_basis
  new_ref_label[idx] <<- proj$reference_label
  new_ref_score[idx] <<- proj$reference_score
  new_cluster[idx] <<- proj$reference_cluster
  if(!"cycling_state"%in%colnames(proj)) proj$cycling_state <- map_cycle[proj$reference_cluster]
  cycling_state_global[idx] <<- proj$cycling_state
  proj_margin[idx] <<- proj$projection_margin
  count_rows[[dataset_name]] <<- data.table(
    dataset=dataset_name,tissue=as.character(md$tissue[idx]),
    broad_celltype=proj$projected_label
  )[,.(n_cells=.N),by=.(dataset,tissue,broad_celltype)]

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
  if(is_reference&&"source_author_annotation"%in%colnames(sk[[]])){
    print(DimPlot(sk,reduction="umap",group.by="source_author_annotation",label=FALSE)+
          ggtitle("nature_xue source author annotation"))
    print(DimPlot(sk,reduction="pca",group.by="source_author_annotation",label=FALSE)+
          ggtitle("nature_xue PCA source author annotation"))
  }
  list(reference_matrix=reference_matrix,reference_labels=reference_labels)
}

pdf(file.path(fig_dir,"task004d_broad_annotation_qc.pdf"),width=12,height=9,onefile=TRUE)
message("Building nature_xue reference")
ref <- run_dataset("nature_xue",TRUE)
xue_ref <- ref$reference_matrix; xue_ref_labels <- ref$reference_labels
if(is.null(xue_ref)||ncol(xue_ref)<2) {
  dev.off()
  stop("Failed to build nature_xue reference")
}
for(d in setdiff(datasets,"nature_xue")){
  message("Annotating ",d)
  run_dataset(d,FALSE,xue_ref,xue_ref_labels)
  gc(verbose=FALSE)
}
dev.off()

if(any(is.na(new_label))) stop("Some cells lack broad annotation")
valid_labels <- c("Hepatocyte/Epithelial","T_cell","NK_cell","B_cell","Plasma_cell",
                  "Monocyte/Macrophage","Neutrophil","Dendritic_cell","Mast_cell",
                  "Endothelial","Fibroblast/Mesenchymal","Uncertain/Mixed")
if(any(!new_label%in%valid_labels))
  stop("Unexpected broad label(s): ",paste(setdiff(unique(new_label),valid_labels),collapse=", "))

old_prelim <- if("project_broad_celltype"%in%colnames(md)) {
  as.character(md$project_broad_celltype)
} else rep(NA_character_,nrow(md))

obj$task004_preliminary_label <- old_prelim
obj$broad_celltype_v2 <- new_label
obj$broad_annotation_confidence <- new_conf
obj$broad_annotation_basis <- new_basis
obj$broad_cluster_id <- new_cluster
obj$broad_cluster_resolution <- "0.8"
obj$broad_reference_label <- new_ref_label
obj$broad_reference_score <- new_ref_score
obj$broad_cycling_state <- cycling_state_global
obj$broad_projection_margin <- proj_margin
obj$annotation_status <- "task004d_broad_annotation_v2"
obj@misc$task004d_broad_annotation <- list(
  source_qs=source_qs,
  shared_feature_file=shared_file,
  xue_author_crosswalk=xue_map_file,
  n_shared_genes=length(shared_genes),
  sketch_threshold_cells=250000L,
  sketch_target_cells=sketch_n,
  projection_block=block_n,
  cluster_resolutions=c(0.4,0.8),
  annotation_resolution=0.8,
  reference_dataset="nature_xue",
  reference_source="source_author_annotation",
  sketch_method="Seurat::SketchData method=LeverageScore; up to 50,000 cells for nature_xue after observed full-cell OOM, or when another >250,000-cell dataset's conservative estimate exceeds 75% of available memory",
  projection_method="PCA-centroid nearest-cluster projection in blocks for all cells when sketching is used",
  method=paste(
    "within-dataset LogNormalize/HVG/PCA/clustering;",
    "cluster markers + canonical lineage programs;",
    "Xue author-label pseudobulk SingleR reference;",
    "Seurat leverage-score sketch up to 50k cells for nature_xue after observed full-cell OOM, or when another >250k-cell dataset's conservative estimate exceeds 75% of available memory;",
    "PCA-centroid projection to all cells in blocks"
  ),
  seed=seed
)

all_clusters <- rbindlist(cluster_rows,fill=TRUE,use.names=TRUE)
all_counts <- rbindlist(count_rows,fill=TRUE,use.names=TRUE)
sketch_summary <- rbindlist(sketch_flags,fill=TRUE,use.names=TRUE)

all_counts[,dataset_tissue_total:=sum(n_cells),by=.(dataset,tissue)]
all_counts[,fraction_within_dataset_tissue:=n_cells/pmax(dataset_tissue_total,1)]
all_counts[,dataset_total:=sum(n_cells),by=dataset]
all_counts[,fraction_within_dataset:=n_cells/pmax(dataset_total,1)]
fwrite(all_counts,file.path(results_dir,"task004d_broad_celltype_counts.csv"))
fwrite(all_clusters,file.path(results_dir,"task004d_cluster_annotation.csv"))
program_all <- rbindlist(program_rows,fill=TRUE,use.names=TRUE)
fwrite(program_all,file.path(results_dir,"task004d_cluster_program_evidence.csv"))
fwrite(sketch_summary,file.path(results_dir,"task004d_sketch_usage.csv"))
if(length(de_rows)){
  markers_all <- rbindlist(de_rows,fill=TRUE,use.names=TRUE)
  fwrite(markers_all,file.path(results_dir,"task004d_cluster_markers.csv.gz"),compress="gzip")
}

uncertain <- all_clusters[final_label=="Uncertain/Mixed"]
fwrite(uncertain,file.path(results_dir,"task004d_uncertain_clusters.csv"))

old_cross <- data.table(
  dataset=as.character(md$dataset),
  task004_preliminary_label=old_prelim,
  broad_celltype_v2=new_label
)[,.(
  n_cells=.N
),by=.(dataset,task004_preliminary_label,broad_celltype_v2)]
fwrite(old_cross,file.path(results_dir,"task004d_task004_vs_v2.csv"))

xue_idx <- which(as.character(md$dataset)=="nature_xue")
xue_author <- if("source_author_annotation"%in%colnames(md)) {
  as.character(md$source_author_annotation[xue_idx])
} else rep(NA_character_,length(xue_idx))
xue_author_broad <- map_author_broad(xue_author)
xue_cross <- data.table(
  source_author_annotation=xue_author,
  source_author_broad=xue_author_broad,
  broad_celltype_v2=new_label[xue_idx]
)[,.(
  n_cells=.N
),by=.(source_author_annotation,source_author_broad,broad_celltype_v2)]
xue_cross[,author_label_total:=sum(n_cells),by=source_author_annotation]
xue_cross[,fraction_within_author_label:=n_cells/pmax(author_label_total,1)]
fwrite(xue_cross,file.path(results_dir,"task004d_xue_author_vs_v2.csv"))

neut_cluster_audit <- all_clusters[
  marker_label=="Neutrophil" | final_label=="Neutrophil" |
    (!is.na(anchor_label)&anchor_label=="Neutrophil") |
    (!is.na(broad_reference_label)&broad_reference_label=="Neutrophil")
]
fwrite(neut_cluster_audit,file.path(results_dir,"task004d_neutrophil_marker_coherence.csv"))

message("Writing ",out_qs)
qs::qsave(obj,out_qs,preset="high",check_hash=TRUE,nthreads=8L)
expected_n <- ncol(obj); expected_ids <- colnames(obj); expected_features <- nrow(obj)
rm(obj); gc(verbose=FALSE)

chk <- qs::qread(out_qs,use_alt_rep=FALSE,nthreads=8L)
if(ncol(chk)!=expected_n||nrow(chk)!=expected_features||!identical(colnames(chk),expected_ids))
  stop("Annotated QS validation failed")
if(any(is.na(chk$broad_celltype_v2))) stop("Missing v2 broad labels after reload")
if(any(as.character(chk$scDblFinder.class)=="doublet"))
  stop("Doublet detected in final singlet broad object")
if(any(!as.character(chk$broad_celltype_v2)%in%valid_labels))
  stop("Invalid final broad labels after reload")

n_uncertain <- sum(as.character(chk$broad_celltype_v2)=="Uncertain/Mixed")
uncertain_fraction <- n_uncertain/ncol(chk)
xue_match <- xue_cross[source_author_broad!="Uncertain/Mixed",
  sum(n_cells[source_author_broad==broad_celltype_v2])/sum(n_cells)]
if(!length(xue_match)||!is.finite(xue_match)) xue_match <- NA_real_

val <- data.table(
  status="VALIDATED",
  source_qs=source_qs,
  annotated_qs=out_qs,
  n_cells=ncol(chk),
  n_features=nrow(chk),
  n_datasets=uniqueN(chk$dataset),
  n_broad_types=uniqueN(chk$broad_celltype_v2),
  n_uncertain=n_uncertain,
  uncertain_fraction=uncertain_fraction,
  xue_author_broad_concordance=xue_match,
  any_doublet_remaining=any(as.character(chk$scDblFinder.class)=="doublet"),
  sketch_datasets=paste(sketch_summary[used_sketch==TRUE,dataset],collapse=";")
)
fwrite(val,file.path(results_dir,"task004d_broad_annotation_validation.csv"))

writeLines(c(
  "# Task 004d broad-annotation stage report",
  "",
  "## Status",
  "COMPLETED",
  "",
  paste0("- Cells annotated: ",format(ncol(chk),big.mark=",")),
  paste0("- Features retained: ",format(nrow(chk),big.mark=",")),
  paste0("- Broad classes represented: ",uniqueN(chk$broad_celltype_v2)),
  paste0("- Uncertain/Mixed: ",format(n_uncertain,big.mark=","),
         " (",sprintf("%.2f%%",100*uncertain_fraction),")"),
  paste0("- Xue author broad-label concordance (excluding unmapped author labels): ",
         ifelse(is.na(xue_match),"NA",sprintf("%.2f%%",100*xue_match))),
  paste0("- Sketch used for: ",
         ifelse(any(sketch_summary$used_sketch),
                paste(sketch_summary[used_sketch==TRUE,dataset],collapse=", "),
                "none")),
  "",
  "Within-dataset broad annotation used LogNormalize, shared-HGNC HVGs, 30-PC PCA,",
  "0.4/0.8 clustering sensitivity, canonical cluster-marker programs and Xue",
  "author-label reference transfer. Resolution 0.8 was used for the working partition.",
  "No malignant-cell call and no cross-dataset batch integration were performed."
),report_path)

message("Task 004d broad annotation stage COMPLETED")
