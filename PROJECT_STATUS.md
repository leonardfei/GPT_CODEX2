# Project Status

## Current task
Task 004b — final eight-cohort merge and QS/H5AD export — READY TO EXECUTE

## Current completion state
Approximately 80% complete.

Completed:
- 8 cohort intermediate RDS objects have already been generated;
- these represent 103 Task 004 source objects;
- expected retained total: 1,490,852 cells;
- source objects remain unchanged;
- resumable merge logic is implemented.

Pending:
1. reuse and validate the 8 cohort intermediate objects;
2. perform final eight-cohort merge;
3. join any `counts.*` layers to one final RNA `counts` layer;
4. validate final merged object;
5. write recoverable final checkpoint;
6. export/validate QS if `qs` is available;
7. export/validate H5AD if `anndataR + rhdf5` are available.

## Mandatory final-merge invariants
- 1,490,852 cells
- 8 datasets
- 194 globally unique project samples
- 132 globally unique project patients
- 1,039,293 Tumor cells
- 451,559 Adjacent cells
- 0 duplicated cell IDs
- exactly one final RNA layer: `counts`
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
2026-09-26 — Task 004b rewritten for direct execution from the existing 8 cohort intermediates.
