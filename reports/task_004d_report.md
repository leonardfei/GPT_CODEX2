# Task 004d report — scDblFinder filtering and corrected broad annotation

## Status
COMPLETED — VALIDATED

## Reproducibility
- Runtime: R version 4.3.3 (2024-02-29); Seurat 5.3.0; scDblFinder 1.16.0; SingleCellExperiment 1.24.0; BiocParallel 1.36.0; SingleR 2.4.1; xgboost 1.7.11.1.
- Phase A/B: `TASK004D_REUSE_CALL_TABLE=TRUE bash scripts/bash/task004de_doublet_and_broad_annotation.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas` (validated all-cell calls reused; no models rerun).
- Phase C: `Rscript scripts/R/task004e_broad_annotation.R --project-root /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas --source-qs objects/merge/HCC_TA_8datasets_singlets_v1.qs --out-qs objects/merge/HCC_TA_8datasets_singlets_broad_v1.qs --shared-features results/task004c_shared_hgnc_features_8of8.txt --xue-map config/task004d_xue_author_to_broad.tsv --sketch-cells 50000 --projection-block 2000 --seed 40500`.
- Final QC reconciliation: `Rscript scripts/R/task004d_reconcile_cluster_annotations.R /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`; reused the saved marker table and per-dataset Seurat caches, without rerunning doublet calls, clustering, or marker tests.
- Broad QC figure-only re-render: `Rscript scripts/R/task004d_render_annotation_qc.R /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`; reused saved Seurat reductions and corrected cluster labels, with a compact bottom legend for the many nature_xue author labels.
- Annotation parameters: LogNormalize scale factor 10,000; up to 3,000 HVGs after mitochondrial/ribosomal exclusion; 30-PC PCA; resolutions 0.4/0.8 with 0.8 used; leverage sketch target 50,000; projection blocks 2,000; seed 40500.
- H5AD interoperability export was not generated: the existing Task004c exporter targets an in-memory standard RNA counts matrix, while this final QS preserves counts in the project's custom chunked representation; safe direct reuse was not verified. The validated deliverable is the QS object with raw counts and v2 metadata.

## Input and output objects
- Input HGNC merged QS: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_merged_review_v1.qs (4,562,511,399 bytes; SHA256 1ac8c7ba371879509fa44b1cb2626f0e85fae5da12a035dee27879a097eb6352).
- Singlet-only QS: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_singlets_v1.qs (4,038,274,070 bytes; SHA256 f2544eabb572c9427242df29f2dd3d901a78586ae47aedc74fd789782ccf2b11).
- Final singlet+broad QS: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_singlets_broad_v1.qs (4,053,789,818 bytes; SHA256 85abf8b3e5dd442dba5220fb770470d5694811f0d77cfa24124df431b6b81793).

## scDblFinder
- scDblFinder version: 1.16.0.
- Detection unit: each project_sample_id independently.
- Rate API compatibility: this scDblFinder release lacks dbr.per1k; per sample it received the explicit equivalent dbr=min(1, 0.008*n_cells/1000), preserving the approved 0.8%-per-1,000-cells expectation (not the package default). dbr.sd=NULL, iter=2, SerialParam().
- Resume record: all-cell calls from the preceding 194-capture scoring pass were reused only after validating cell order, capture/sample identifiers, score completeness, package-rate API, and per-sample expected rates; models were not recomputed.
- Samples with >=500 cells used cluster-based mode; smaller samples used random artificial-doublet mode.
- Seeds: 44000 + sample_index.
- Cells before: 1,490,852.
- Predicted doublets removed: 115,725 (7.76%).
- Singlets retained: 1,375,127.
- Preliminary Task004 neutrophil retention: 90.70% (59,452/65,546).
- Samples with doublet-rate review flags: 0.
- CRA002308 source preprocessing reported doublet exclusion; the same project scDblFinder audit was still run because no per-cell known-doublet truth labels are available.
- nature_xue is an author-processed object: all author labels were preserved, and scDblFinder was applied as an additional project-level computational filter.

### Per-dataset doublet summary
- GSE282701: 13,626/137,518 doublets (9.91%); 123,892 singlets retained.
- GSE242889: 3,744/52,954 doublets (7.07%); 49,210 singlets retained.
- GSE326201: 7,139/92,964 doublets (7.68%); 85,825 singlets retained.
- GSE149614: 3,544/63,101 doublets (5.62%); 59,557 singlets retained.
- GSE299340: 6,805/76,319 doublets (8.92%); 69,514 singlets retained.
- CRA002308: 17,776/158,741 doublets (11.20%); 140,965 singlets retained.
- nature_xue: 38,684/675,539 doublets (5.73%); 636,855 singlets retained.
- in_house: 24,407/233,716 doublets (10.44%); 209,309 singlets retained.

## Broad annotation v2
- Annotation was performed within each dataset without cross-dataset batch integration.
- Feature universe started from the Task004c 8/8 shared HGNC set; mitochondrial and ribosomal HVGs were excluded from PCA features.
- LogNormalize (scale factor 10,000), ~3,000 HVGs, ScaleData, 30-PC PCA, neighbors and clustering at resolutions 0.4 and 0.8 were used; resolution 0.8 is the working partition. After a documented full-cell OOM in nature_xue, Seurat leverage-score sketching (up to 50,000 cells) was used for that dataset; the same method is available to any >250,000-cell dataset exceeding the conservative memory guard. Sketch cluster labels were projected to all cells by block-wise PCA-centroid distance.
- Cluster markers and canonical lineage programs were combined with Xue author-label pseudobulk reference transfer.
- Neutrophil calls required coherent granulocyte evidence including at least one core marker; S100A8/S100A9 alone were insufficient.
- No malignant-cell call and no cross-dataset integration were performed.
- Uncertain/Mixed: 225,604 (16.41%).
- Xue broad author-label concordance (excluding unmapped author labels): 83.31%.

### Final broad cell counts
- T_cell: 459,301 cells (33.40%).
- Uncertain/Mixed: 225,604 cells (16.41%).
- NK_cell: 175,184 cells (12.74%).
- Monocyte/Macrophage: 131,316 cells (9.55%).
- Hepatocyte/Epithelial: 123,432 cells (8.98%).
- Endothelial: 115,262 cells (8.38%).
- Fibroblast/Mesenchymal: 46,886 cells (3.41%).
- Neutrophil: 37,948 cells (2.76%).
- B_cell: 34,158 cells (2.48%).
- Plasma_cell: 22,968 cells (1.67%).
- Mast_cell: 3,068 cells (0.22%).

### Scalability / sketch use
- nature_xue: 50,000 cells used for clustering of 636,855 retained cells; sketch=YES; method=Seurat::SketchData method=LeverageScore.
- GSE282701: 123,892 cells used for clustering of 123,892 retained cells; sketch=NO; method=none_full_cell.
- GSE242889: 49,210 cells used for clustering of 49,210 retained cells; sketch=NO; method=none_full_cell.
- GSE326201: 85,825 cells used for clustering of 85,825 retained cells; sketch=NO; method=none_full_cell.
- GSE149614: 59,557 cells used for clustering of 59,557 retained cells; sketch=NO; method=none_full_cell.
- GSE299340: 69,514 cells used for clustering of 69,514 retained cells; sketch=NO; method=none_full_cell.
- CRA002308: 140,965 cells used for clustering of 140,965 retained cells; sketch=NO; method=none_full_cell.
- in_house: 209,309 cells used for clustering of 209,309 retained cells; sketch=NO; method=none_full_cell.

## Key QC outputs
- results/task004d_scdblfinder_cell_calls.csv.gz
- results/task004d_scdblfinder_by_sample.csv
- results/task004d_scdblfinder_by_dataset.csv
- results/task004d_scdblfinder_by_dataset_tissue.csv
- results/task004d_scdblfinder_overall.csv
- results/task004d_scdblfinder_cell_qc_distributions.csv
- results/task004d_neutrophil_doublet_retention.csv
- results/task004d_neutrophil_retention_distributions.csv
- results/task004d_broad_celltype_counts.csv
- results/task004d_cluster_annotation.csv
- results/task004d_cluster_markers.csv.gz
- results/task004d_xue_author_vs_v2.csv
- results/task004d_task004_vs_v2.csv
- results/task004d_uncertain_clusters.csv
- results/task004d_neutrophil_marker_coherence.csv
- figures/task004d_scdblfinder_qc.pdf
- figures/task004d_broad_annotation_qc.pdf (required review copy; generated source: figures/task004e/task004d_broad_annotation_qc.pdf)

## Deviations / errors
No sample-level doublet-rate review flag was triggered.
Resolved runtime issues: the server's scDblFinder 1.16 API required the explicit legacy expected-rate equivalent; Seurat accessor validation was namespaced; saved gzip calls were streamed without optional R.utils; data.table singlet selection and grouped median output types were corrected. The validated first all-cell scoring pass was preserved and reused. A full-cell nature_xue attempt reached >135 GB RSS and was OOM-killed; the approved Seurat leverage-score sketch (50,000 cells) and block-wise PCA-centroid projection were then used. Sketch cells/counts/metadata were aligned explicitly by cell ID, and mapped cluster-label vectors were made unnamed before Seurat metadata assignment to prevent its names-as-cell-IDs overlap error. Seurat's optional `presto` acceleration was not installed; the standard Wilcoxon marker test completed.
During final QC, a scoping error was found in the original core-marker audit: neutrophil-core and CD3 hits had been copied from the last cluster into every cluster row. Per-cluster hit counts and all downstream reference decisions were recomputed from the preserved marker table and Seurat analysis caches using the original thresholds; doublet calls, cluster definitions, and raw counts were not changed.
No secondary UMI/nFeature doublet cutoff was used. No predicted doublet was rescued because of a preliminary cell-type label.

## Hold point
Task 005 remains paused. Review Task 004d QC and broad annotation before integration.
