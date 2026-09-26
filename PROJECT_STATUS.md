# Project Status

## Current task
Task 004b — resumable eight-cohort merge with native QS + H5AD export — READY TO RUN

## Revised strategy
The previous preflight blocker has been removed from the expensive merge stage.

Task 004b now separates:
1. core merge;
2. QS export;
3. H5AD export.

The 1,490,852-cell merge requires only the already available core R packages (Seurat, data.table, Matrix). After merge, a temporary checkpoint RDS is written. If an export package is missing, the checkpoint is retained so that a later rerun performs export only rather than repeating the 103-object merge.

## Final target outputs
Seurat:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs`

AnnData:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`

## H5AD strategy
Use `anndataR` + `rhdf5` directly from the merged Seurat object. No SingleCellExperiment, zellkonverter, reticulate, Python anndata or pandas environment is required.

H5AD mapping:
- AnnData X = RNA counts
- obs = Seurat cell metadata
- obs_names = globally unique cell IDs
- var_names = genes/features
- no reductions/graphs copied

## Additional package needs
- QS: qs
- H5AD: anndataR + rhdf5

If these are missing, merge still proceeds and the checkpoint is retained. A rerun automatically detects the checkpoint and skips the expensive 103-object merge. The checkpoint is removed only after both QS and H5AD validate successfully.

## Expected data
- 8 cohorts
- 103 annotated source objects
- 1,490,852 cells
- Task 004 labels retained as preliminary/unvalidated

## Next execution command
`Execute task_004b.`

## Last update
2026-09-26 — export architecture revised to resumable merge + native anndataR H5AD.
