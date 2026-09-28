# Task 004d/e — scDblFinder filtering followed by broad cell annotation

## Status
READY TO EXECUTE

## User-authorized source
/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_merged_review_v1.qs

## Task 004d — doublet detection and filtering
Run scDblFinder independently for each project_sample_id using raw counts from the merged QS.

Parameters:
- automatic expected doublet rate;
- dbr.per1k = 0.008;
- nfeatures = 1352;
- dims = 20;
- random artificial doublets (clusters = NULL);
- deterministic per-sample seeds;
- samples with fewer than 100 cells are retained and flagged as unscored rather than force-called.

Safety:
- preliminary Task 004 labels are not used for detection;
- they are used only for post hoc neutrophil/cell-type retention auditing;
- if any scored sample has >30% predicted doublets, stop before filtering for review.

Output:
objects/merge/HCC_TA_8datasets_merged_scdblfinder_filtered_v1.qs

Required audits:
- per-cell score/class table;
- per-sample and per-dataset doublet rates;
- preliminary broad-cell-type doublet rate;
- neutrophil retention;
- reload validation.

## Task 004e — broad cell annotation
Run only after Task 004d validation.

Do not reuse the old Task 004 per-cell marker-score assignment as the reference annotation.

Broad categories:
- Tumor/epithelial
- T/NK
- B
- Plasma
- Monocyte/macrophage
- Neutrophil
- Dendritic
- Fibroblast/mesenchymal
- Endothelial
- Mast
- Erythroid
- Other/uncertain

Method:
1. use the strict 8/8 HGNC shared-gene set from Task 004c;
2. process each dataset separately;
3. draw a stratified sketch of up to 50,000 singlets per dataset;
4. LogNormalize, select 3,000 HVGs, PCA, kNN, clustering at resolution 0.6, UMAP;
5. calculate cluster markers;
6. use exact Xue author prefixes as the primary broad-cell reference anchors:
   Neu_*, Mph_*/Mo_*/Mono-like_*, DC_*/MonoDC, EC_*, Fb_*/Mu_*,
   B_*, CD4T_*/CD8T_*/NK_*/gdT_*, Mast, and Tumor;
7. construct Xue sample-by-broad-type reference pseudobulks;
8. transfer labels to other datasets at cluster level with SingleR;
9. require agreement with canonical marker programs when possible and retain conflicts as Other/uncertain;
10. project all retained singlets to the nearest dataset-specific sketch PCA cluster in blocks, preserving every retained cell.

The previous Task 004 broad label is copied to task004_preliminary_broad_celltype.
The new annotation is written to project_broad_celltype_v2.

Output:
objects/merge/HCC_TA_8datasets_merged_scdblfinder_filtered_broad_v1.qs

## Execution
cd /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
bash scripts/bash/task004de_doublet_and_broad_annotation.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Stop rule
Stop after the new broad-annotated QS validates. Do not start Task 005 integration automatically.
