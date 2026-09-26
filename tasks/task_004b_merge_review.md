# Task 004b — final QS + H5AD export

## Status
MERGE VALIDATED; PACKAGE DOWNLOADS AUTHORIZED; EXPORT READY

The user explicitly authorizes downloading and installing the packages needed for export.

Do not redo the 103 source objects or the validated 1,490,852-cell final merge.

### Frozen merge invariants
- 1,490,852 cells
- 68,394 union features
- 8 datasets
- 194 globally unique project samples
- 132 globally unique project patients
- 1,039,293 Tumor cells
- 451,559 Adjacent cells
- 0 duplicated cell IDs
- checkpoint: `objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds`

## QS export
The user requested legacy `.qs`.
Install archived `qs 0.27.3` into `.task004b_Rlib`; do not substitute `qs2`, because its `.qs2` format is not compatible with legacy `.qs`.

Write and reload-validate:
`objects/HCC_TA_8datasets_merged_review_v1.qs`

## H5AD export
Do not convert the custom final chunked R counts layer directly.

Use:
1. `anndataR + rhdf5` in an isolated project R library to write the 8 standard cohort Seurat checkpoints as temporary H5AD files;
2. global cell IDs `dataset::cell_id`;
3. RNA counts -> AnnData X;
4. all harmonized cell metadata -> obs;
5. no reductions/graphs;
6. Python AnnData `experimental.concat_on_disk`, `axis="obs"`, `join="outer"`, `fill_value=0`, `index_unique=None`.

This avoids creating a single R `dgCMatrix` and permits the final sparse H5AD to use 64-bit on-disk sparse pointers when required.

Final:
`objects/HCC_TA_8datasets_merged_review_v1.h5ad`

Validate dimensions, global IDs, dataset/sample/patient counts, tissue cell counts, annotation status, unique obs/var names, sparse X, and X nnz equality to the eight input matrices.

## Isolated environments
- R: `.task004b_Rlib`
- Python: `.task004b_pyenv`

Do not update the main analysis environment. If system Python is <3.11, automatically create an isolated Python 3.12 environment with an available conda installation (including likely project/local Miniconda paths).

## Disk
Require >=60 GB free before export. Remove temporary cohort H5AD parts only after final H5AD validates. Retain the 36.6-GB merge checkpoint for rollback until the user explicitly approves deletion.

## Execute
```bash
cd /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
bash scripts/bash/task004b_install_and_export.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
```

## Required final outputs
- `objects/HCC_TA_8datasets_merged_review_v1.qs`
- `objects/HCC_TA_8datasets_merged_review_v1.h5ad`
- `results/task004b_qs_validation.csv`
- `results/task004b_h5ad_parts_manifest.csv`
- `results/task004b_h5ad_validation.json`

Task 004b is COMPLETED only when both final files validate. Task 005 remains paused.
