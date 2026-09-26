# Task 004b report — eight-cohort unintegrated Seurat merge for annotation review

## Status

Task 004b status: COMPLETED
Final merge validation: VALIDATED
QS export: VALIDATED
H5AD export: VALIDATED

## Purpose

This object was created for manual inspection of the current Task 004 annotations. No RPCA/Harmony integration, batch correction, normalization, PCA, UMAP, or reclustering was performed.

## Input

- 103 Task 004 annotated Seurat objects
- 8 cohorts
- Total cells: 1,490,852
- Samples: 194 (validated using globally unique project_sample_id)
- Patients: 132 (validated using globally unique project_patient_id)
- Explicit paired Tumor-Adjacent patients: 59 (metadata design reference)
- Tumor cells: 1,039,293
- Adjacent cells: 451,559

## Merge content

- RNA counts and cell-level metadata were retained.
- The final object contains 68,394 features from the union of cohort feature sets; cohort-specific feature order was aligned and absent features were represented as sparse zeros.
- The single RNA `counts` layer is stored as chunked sparse blocks because one standard `dgCMatrix` exceeded Matrix's 32-bit cumulative non-zero index limit; the layer coordinates remain globally aligned to the 68,394 features and 1,490,852 cells.
- Pre-existing reductions/graphs were intentionally discarded because they are not directly comparable across independently processed source objects.
- Current project_broad_celltype and neutrophil_confidence are retained only for review and are marked annotation_status=preliminary_unvalidated_task004.
- Xue author labels remain available through source_author_annotation.
- The known Task 004 rescued-status bug was corrected for CRA002308 and in_house in this review object only.

## Output

Seurat QS target: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs (4,754,019,531 bytes; server)
AnnData H5AD target: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad (38,153,553,374 bytes; server)
QS status: VALIDATED
H5AD status: VALIDATED
QS validation record: [task004b_qs_validation.csv](../results/task004b_qs_validation.csv)
H5AD parts manifest: [task004b_h5ad_parts_manifest.csv](../results/task004b_h5ad_parts_manifest.csv)
H5AD validation record: [task004b_h5ad_validation.json](../results/task004b_h5ad_validation.json)
Recoverable merge checkpoint: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds

The final merge checkpoint is retained at approximately 36.6 GB for rollback. The eight temporary H5AD parts were removed only after final validation. Raw/source data were not modified.

## Export reproducibility and anomalies

- Server R runtime: R 4.3.3; Seurat 5.3.0.
- Export R packages: qs 0.27.3, stringfish 0.18.0, anndataR 0.2.0, rhdf5 2.46.1, rhdf5filters 1.14.1, Rhdf5lib 1.24.2.
- H5AD runtime: Python 3.14; anndata 0.13.4, h5py 3.16.0, scipy 1.18.1, pandas 3.0.6, numpy 2.5.3.
- H5AD validation: `X` is `csr_matrix`; `nnz=2,885,352,168`, matching the sum of the eight part manifests; `join=outer`, `fill_value=0`.
- Compatibility repairs were implementation-only: anndataR 0.2.0 does not accept the newer `chunk_size` argument; its HDF5AnnData dimension accessor is `shape()` rather than `dim()`; and the anndata temporary concat path must retain a `.h5ad` suffix to select the HDF5 backend. The final validator uses HDF5 axis encodings directly and does not require optional `xarray` lazy loading.

## Important limitation

This is a pure merge object, not an integrated atlas. Dataset-driven structure is expected if the merged counts are normalized/PCA/UMAPed without batch correction. Task 005 remains unexecuted.
