# Task 004b — Merge eight cohorts and export QS + H5AD review objects

## Status
READY TO RUN

## Goal
Create two equivalent unintegrated review objects containing all eight cohorts so the user can inspect the current Task 004 annotations:
1. a Seurat object in `.qs` format;
2. an AnnData object in `.h5ad` format.

## Scope
This task is merge/export only. Do not run RPCA, Harmony, CCA, normalization, PCA, UMAP, reclustering, annotation changes, cell removal, or downsampling.

Use the 103 Task 004 annotated Seurat objects as inputs. Retain all 1,490,852 cells, RNA counts, and metadata. Source-specific PCA/UMAP/graphs must be discarded because they are not directly comparable across independently processed objects.

The current `project_broad_celltype` and `neutrophil_confidence` fields must be retained strictly as preliminary review labels. Add:
`annotation_status = preliminary_unvalidated_task004`.

Preserve `source_author_annotation` from nature_xue.

Correct the known Task 004 rescued-status metadata bug in the exported review objects only:
- CRA002308 and in_house: `rescued_by_corrected_qc=FALSE`; `qc_status=corrected_qc_pass`.
- nature_xue: `rescued_by_corrected_qc=FALSE`; `qc_status=author_processed`.
- Preserve the original Task 004 values in `task004_qc_status_premerge` and `task004_rescued_premerge`.

## Required software
Before the expensive merge, verify these R packages are already available:
- Seurat
- data.table
- Matrix
- qs
- SingleCellExperiment
- zellkonverter

For H5AD validation, Python must have:
- anndata
- pandas

If export packages are missing, stop and report the missing packages. Do not silently install from the compute server.

## Primary outputs

Seurat:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs`

AnnData:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`

No final `.rds` output is required.

The H5AD must contain:
- RNA counts as `adata.X`;
- genes/features as `adata.var_names`;
- globally unique cell IDs as `adata.obs_names`;
- harmonized cell metadata as `adata.obs`.

## Required tracked outputs
- `results/task004b_merge_review_cohort_summary.csv`
- `results/task004b_merge_review_dataset_summary.csv`
- `results/task004b_merge_review_metadata_fields.csv`
- `results/task004b_merge_review_validation.csv`
- `results/task004b_h5ad_validation.json`
- `reports/task_004b_merge_review_report.md`

## Validation

### QS
Reload with `qs::qread()` and require:
- object inherits from Seurat;
- exactly 1,490,852 cells;
- ordered cell IDs unchanged;
- exactly 8 datasets;
- RNA counts retained.

### H5AD
Open using Python `anndata.read_h5ad(..., backed="r")` and require:
- exactly 1,490,852 observations;
- same feature count as the QS object;
- exactly 8 datasets;
- unique `obs_names`;
- required `obs` metadata fields present;
- RNA counts stored in `X`.

## Execution
Run:

`Rscript scripts/R/task004b_merge_review.R --project-root /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`

Then:

`python scripts/python/task004b_validate_h5ad.py --h5ad /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad --out /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/results/task004b_h5ad_validation.json`

Stop after validation. Do not execute Task 005.
