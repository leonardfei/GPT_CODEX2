# Project Status

## Current task
Task 004b — export validated eight-cohort merge to QS + H5AD — COMPLETED

The final merge is already VALIDATED:
- 1,490,852 cells
- 68,394 features
- 8 datasets
- 194 samples
- 132 patients
- 1,039,293 Tumor cells
- 451,559 Adjacent cells
- zero duplicate cell IDs

The user has explicitly authorized downloading export dependencies.

Completed export validation:
- QS: `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs`, 4,754,019,531 bytes, qs 0.27.3.
- H5AD: `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`, 38,153,553,374 bytes, anndata 0.13.4.
- H5AD shape: 1,490,852 observations × 68,394 variables; `X` is CSR with 2,885,352,168 non-zero entries.
- Cohort/sample/patient counts: 8/194/132; Tumor/Adjacent counts: 1,039,293/451,559.
- Validation records committed under `results/task004b_*.csv` and `results/task004b_*.json`.
- Eight temporary H5AD parts were removed after final validation; the merge checkpoint remains retained.

## Export architecture
QS:
- install archived qs 0.27.3 in isolated `.task004b_Rlib`
- write/reload validate `objects/HCC_TA_8datasets_merged_review_v1.qs`

H5AD:
- install anndataR + rhdf5 in isolated R library
- write eight standard cohort H5AD parts
- install AnnData in isolated `.task004b_pyenv`
- concatenate on disk with outer gene union and sparse zero fill
- validate final `objects/HCC_TA_8datasets_merged_review_v1.h5ad`

This avoids converting the custom 1.49M-cell chunked counts layer into one R dgCMatrix.

## Next execution
No further Task 004b execution is required. Task 005 remains paused.

Equivalent:
`bash scripts/bash/task004b_install_and_export.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`

Task 005 remains paused.
