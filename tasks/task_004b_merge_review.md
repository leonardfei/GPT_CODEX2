# Task 004b — Best-practice eight-cohort merge and QS/H5AD export

## Status
READY TO RUN / RESUMABLE

## Goal
Create two equivalent unintegrated review objects containing all eight cohorts:
1. Seurat `.qs`;
2. AnnData `.h5ad`.

## Engineering strategy

The expensive 1,490,852-cell merge is decoupled from format export.

### Phase A — merge, always possible with current core R environment
Requirements:
- Seurat
- data.table
- Matrix

Actions:
- read the 103 Task 004 annotated objects;
- rebuild each as a counts+metadata review object;
- merge to 8 cohorts and then one 1,490,852-cell Seurat object;
- preserve RNA raw counts and metadata;
- discard source-specific reductions/graphs;
- write a recoverable checkpoint:

`objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds`

If this checkpoint already exists and validates, subsequent runs skip the expensive 103-object merge and perform export only.

### Phase B — QS export
Requires only:
- `qs`

Output:
`objects/HCC_TA_8datasets_merged_review_v1.qs`

Reload with `qs::qread()` and verify class, dimensions and ordered cell IDs.

### Phase C — H5AD export
Use native R AnnData interoperability:
- `anndataR`
- `rhdf5`

Do **not** use SingleCellExperiment, zellkonverter, reticulate, or a Python export environment.

`anndataR::write_h5ad()` writes directly from the merged Seurat object with:
- assay_name = RNA
- x_mapping = counts
- layers_mapping = FALSE
- obs_mapping = TRUE
- reductions/graphs/misc excluded
- gzip compression

Output:
`objects/HCC_TA_8datasets_merged_review_v1.h5ad`

Validate by reopening as an `HDF5AnnData` object using `anndataR::read_h5ad(..., as="HDF5AnnData")` and checking observations/features.

After both QS and H5AD validate, delete the temporary checkpoint. If either export dependency is missing, retain the checkpoint and report a partial/blocker state; do not redo the merge after packages become available.

## Biological/data rules
No RPCA/Harmony/CCA, normalization, PCA, UMAP, reclustering, annotation change, cell removal or downsampling.

Keep:
- 1,490,852 cells;
- all 8 datasets;
- RNA raw counts;
- harmonized metadata;
- preliminary Task 004 labels for inspection;
- nature_xue source author annotation.

Known rescued-status metadata correction:
- CRA002308/in_house: rescued_by_corrected_qc=FALSE; qc_status=corrected_qc_pass
- nature_xue: rescued_by_corrected_qc=FALSE; qc_status=author_processed
- preserve premerge values for audit.

## Dependencies
Core merge:
- Seurat
- data.table
- Matrix

Final exports:
- QS: `qs`
- H5AD: `anndataR` + `rhdf5`

These are the only additional format dependencies.

If additional packages are unavailable, do not silently install them. The merge checkpoint allows export-only continuation later.

## Primary outputs
- `objects/HCC_TA_8datasets_merged_review_v1.qs`
- `objects/HCC_TA_8datasets_merged_review_v1.h5ad`

## Validation outputs
- `results/task004b_merge_review_cohort_summary.csv`
- `results/task004b_merge_review_dataset_summary.csv`
- `results/task004b_merge_review_metadata_fields.csv`
- `results/task004b_merge_review_validation.csv`
- `reports/task_004b_merge_review_report.md`

Python validation is no longer required.

## Execution
`Rscript scripts/R/task004b_merge_review.R --project-root /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`

Stop after both exports validate or after a partial export state with checkpoint retained. Do not execute Task 005.
