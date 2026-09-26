# Project Status

## Current task
Task 004b — final eight-cohort merge and QS/H5AD export — MERGE_COMPLETED_EXPORT_PARTIAL

## Current completion state
Phase A is complete and validated; requested QS/H5AD exports remain blocked by missing server dependencies.

Completed:
- 8 cohort intermediate RDS objects have already been generated;
- these represent 103 Task 004 source objects;
- expected retained total: 1,490,852 cells;
- source objects remain unchanged;
- final recoverable checkpoint was written on the server;
- merge validation passed for all mandatory cell, metadata, ID and layer invariants.

Pending:
1. install/provide `qs` on the server and rerun export-only;
2. install/provide `anndataR + rhdf5` on the server and rerun export-only.

## Mandatory final-merge invariants
- 1,490,852 cells
- 8 datasets
- 194 globally unique project samples
- 132 globally unique project patients
- 1,039,293 Tumor cells
- 451,559 Adjacent cells
- 0 duplicated cell IDs
- exactly one final RNA layer: `counts` (validated; stored as chunked sparse blocks because a single standard `dgCMatrix` exceeds Matrix's 32-bit non-zero index limit)
- preserve source-local IDs and add globally unique project sample/patient/paired IDs

## Final target outputs
Seurat:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs`

AnnData:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`

Recoverable checkpoint:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds`

## Export dependencies
- QS: `qs`
- H5AD: `anndataR` + `rhdf5`

Missing export dependencies must not block or invalidate the final merge. If they are absent, finish Phase A, retain the checkpoint and report `MERGE_COMPLETED_EXPORT_PARTIAL`.

## Completion rule
Task 004b is `COMPLETED` only if final merge + QS validation + H5AD validation all pass.

Task 005 remains paused.

## Next execution command
`Execute task_004b.`

## Last update
2026-09-26 — Task 004b Phase A completed; status is `MERGE_COMPLETED_EXPORT_PARTIAL` because `qs`, `anndataR` and `rhdf5` are unavailable on the server.
