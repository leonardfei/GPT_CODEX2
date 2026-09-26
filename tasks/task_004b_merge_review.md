# Task 004b — Final eight-cohort merge and QS/H5AD export

## Status
READY TO EXECUTE — cohort-level merge already available; final merge/export pending.

## Current server state
Already completed:
- 8 cohort intermediate RDS objects exist under `objects/task004_merge_review/cohort_merged/`;
- the 8 cohorts collectively represent 103 Task 004 source objects;
- expected retained cell total is 1,490,852;
- source Task 004 objects must not be modified.

Do not rebuild the 103 source objects unless an existing cohort checkpoint fails structural validation.

## Goal
Produce one validated unintegrated review object in two formats:
1. Seurat QS:
   `objects/HCC_TA_8datasets_merged_review_v1.qs`
2. AnnData H5AD:
   `objects/HCC_TA_8datasets_merged_review_v1.h5ad`

This task is still pre-integration. Do not execute Task 005.

## Phase A — final merge (mandatory)
Use the 8 existing cohort intermediate RDS files first.

For every cohort checkpoint:
- verify it is a Seurat object;
- verify expected dataset identity and expected cell count;
- verify globally unique cell IDs within the cohort;
- retain RNA counts and metadata;
- if multiple `counts.*` layers exist, join only count layers and rebuild a clean counts-only Seurat object;
- do not retain PCA/UMAP/graphs.

Merge the 8 validated cohort objects.

The final object must satisfy all of the following before any format export is considered:
- exactly 1,490,852 cells;
- exactly 8 datasets;
- exactly 194 sample IDs;
- exactly 132 patient IDs;
- 1,039,293 Tumor cells;
- 451,559 Adjacent cells;
- zero duplicated cell IDs;
- RNA assay exists;
- final RNA assay contains exactly one layer named `counts`;
- `project_broad_celltype`, `neutrophil_confidence`, `source_author_annotation`, dataset/sample/patient/tissue and QC provenance metadata remain present;
- Task 004 labels remain marked `annotation_status=preliminary_unvalidated_task004`.

After validation, write the recoverable final merge checkpoint:

`objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds`

Use `compress=FALSE` for this temporary checkpoint to avoid avoidable serialization memory peaks.

If this checkpoint already exists on rerun and validates, skip cohort/final merge and go directly to export.

## Phase B — QS export
If R package `qs` is available:
- export the checkpoint/final merged Seurat object to:
  `objects/HCC_TA_8datasets_merged_review_v1.qs`
- reload using `qs::qread()`;
- require Seurat class;
- require exactly 1,490,852 cells;
- require ordered cell IDs identical to the merged object.

If `qs` is unavailable:
- do not fail Phase A;
- record `QS status = BLOCKED_MISSING_qs`;
- retain the final merge checkpoint.

## Phase C — H5AD export
Preferred native R route:
- `anndataR`
- `rhdf5`

Use `anndataR::write_h5ad()` directly from the Seurat object:
- `assay_name="RNA"`;
- `x_mapping="counts"`;
- `layers_mapping=FALSE`;
- `obs_mapping=TRUE`;
- do not copy reductions, graphs or misc;
- gzip compression.

Expected AnnData structure:
- `X` = RNA raw counts;
- `obs` = cell metadata;
- `obs_names` = globally unique cell IDs;
- `var_names` = gene/features.

Validate with:
`anndataR::read_h5ad(path, as="HDF5AnnData", mode="r")`

Require:
- 1,490,852 observations;
- feature count equal to the Seurat/QS object.

If `anndataR` or `rhdf5` is unavailable:
- do not fail Phase A;
- record the missing dependency;
- retain the final merge checkpoint.

## Important metadata correction
In the review object only:
- CRA002308 and in_house: `rescued_by_corrected_qc=FALSE`, `qc_status=corrected_qc_pass`;
- nature_xue: `rescued_by_corrected_qc=FALSE`, `qc_status=author_processed`;
- preserve pre-correction Task 004 values in `task004_qc_status_premerge` and `task004_rescued_premerge`.

Do not otherwise change cell annotations or cell inclusion.

## Status definitions
### COMPLETED
Only when:
- final merge validation passes;
- QS is written and reload-validates;
- H5AD is written and backed-validates.

Then the temporary final checkpoint may be removed.

### MERGE_COMPLETED_EXPORT_PARTIAL
When:
- final 1,490,852-cell merge validation passes;
- checkpoint is successfully written;
- one or both requested export formats are blocked by missing packages.

In this state, retain the checkpoint. A later rerun must skip the expensive merge and perform export only.

### FAILED
Any failure in final merge integrity, cell counts, dataset/sample/patient counts, cell ID uniqueness, metadata integrity, or counts-layer integrity.

Do not continue to Task 005.

## Required tracked outputs
- `results/task004b_merge_review_cohort_summary.csv`
- `results/task004b_merge_review_dataset_summary.csv`
- `results/task004b_merge_review_metadata_fields.csv`
- `results/task004b_merge_review_validation.csv`
- `reports/task_004b_merge_review_report.md`

The validation CSV must include:
- `task004b_status`;
- `merge_validation_status`;
- QS/H5AD statuses and file sizes;
- cells/features/datasets/samples/patients;
- duplicate cell count;
- RNA layers;
- annotation/provenance field checks.

## Execution
From the project root:

```bash
cd /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

Rscript scripts/R/task004b_merge_review.R \
  --project-root /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
```

## Stop rule
Stop after Task 004b reaches either:
- `COMPLETED`, or
- `MERGE_COMPLETED_EXPORT_PARTIAL`.

Do not execute Task 005.
