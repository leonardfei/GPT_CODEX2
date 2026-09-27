# Project Status

## Current task
Task 004c — in-place HGNC harmonisation of merged QS/H5AD — READY TO EXECUTE

Task 004b is COMPLETED:
- 1,490,852 cells
- 68,394 raw union features before harmonisation
- 8 datasets
- 194 samples
- 132 patients

## Task 004c strategy
The user explicitly requested direct replacement of the existing merged files.

Target paths remain:
- objects/HCC_TA_8datasets_merged_review_v1.qs
- objects/HCC_TA_8datasets_merged_review_v1.h5ad

Task 004c first creates temporary HGNC-standardised versions, validates both, then atomically overwrites these two original paths. If either validation fails, the originals remain untouched.

After successful replacement, the merged v1 paths themselves will contain HGNC-approved gene symbols with duplicate mappings collapsed by exact raw-count summation.

Unmapped/custom/non-human features will no longer be present in the overwritten merged objects, but remain recoverable from upstream cohort/source objects and the retained Task 004b merge checkpoint.

## Next command
Execute task_004c.

Equivalent:
bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

Task 005 remains paused pending Task 004c review.