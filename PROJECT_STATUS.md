# Project Status

## Current task
Task 004b — export validated eight-cohort merge to QS + H5AD — READY TO EXECUTE

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
`Execute task_004b export.`

Equivalent:
`bash scripts/bash/task004b_install_and_export.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`

Task 005 remains paused.
