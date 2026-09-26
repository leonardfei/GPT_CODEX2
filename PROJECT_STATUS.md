# Project Status

## Current task
Task 004b — eight-cohort unintegrated merge/export to QS + H5AD — BLOCKED

## Purpose
At user request, create two equivalent review objects containing all eight cohorts before correcting Task 004 annotation or running Task 005 integration.

## Target outputs

Seurat QS:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs`

AnnData H5AD:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`

Expected:
- 8 cohorts
- 103 Task 004 annotated source objects
- 1,490,852 cells

No final RDS review object is required.

## Export semantics
- QS contains the full unintegrated Seurat counts+metadata review object.
- H5AD contains the same cells/features with RNA counts in `X` and metadata in `obs`.
- Neither output is an integrated atlas.
- No normalization, PCA, UMAP, reclustering, RPCA/Harmony/CCA, removal or downsampling is performed.

## Task 004 annotation status
Task 004 broad annotation remains preliminary/unvalidated and is retained only so the user can inspect it.

## Blocker
Task 004b cannot start because the server-side export/validation environment is missing `qs`, `SingleCellExperiment`, `zellkonverter`, `anndata` and `pandas`. The task requires stopping before the expensive merge when these packages are unavailable; no package installation was attempted.

## Pending
1. Provide or enable the required export/validation environment
2. Execute Task 004b
3. Validate QS by reload and H5AD with backed AnnData
4. User manually reviews annotation
5. Rebuild/correct annotation as needed
6. Run Task 005 integration only after annotation review

## Next execution command
After the environment blocker is resolved: `Execute task_004b.`

## Last update
2026-09-26 — Task 004b preflight blocked before merge/export; see `reports/task_004b_merge_review_report.md`.
