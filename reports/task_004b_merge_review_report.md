# Task 004b report — eight-cohort unintegrated Seurat merge for annotation review

## Status

Task 004b status: MERGE_COMPLETED_EXPORT_PARTIAL
Final merge validation: VALIDATED
QS export: BLOCKED_MISSING_qs
H5AD export: BLOCKED_MISSING_anndataR+rhdf5

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

Seurat QS target: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs (not created)
AnnData H5AD target: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad (not created)
QS status: BLOCKED_MISSING_qs
H5AD status: BLOCKED_MISSING_anndataR+rhdf5
QS SHA256: NOT_CREATED
H5AD SHA256: NOT_CREATED
Recoverable merge checkpoint: /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds

The final merge checkpoint is retained at approximately 36.6 GB for export-only resume. QS export requires `qs`; H5AD export requires `anndataR` and `rhdf5`. None of these packages is installed in the server R environment, so no format export was attempted.

## Important limitation

This is a pure merge object, not an integrated atlas. Dataset-driven structure is expected if the merged counts are normalized/PCA/UMAPed without batch correction. Task 005 remains unexecuted.
