# Project Status

## Current task
Task 004c — harmonise gene identifiers directly in merged QS/H5AD — READY TO EXECUTE

Task 004b remains COMPLETED and its validated raw-union outputs are frozen:
- 1,490,852 cells
- 68,394 union features
- 8 datasets
- 194 samples
- 132 patients

## Task 004c strategy
Do not create eight separate harmonised cohort objects.

Instead, use the validated merged files directly:
- source QS: objects/HCC_TA_8datasets_merged_review_v1.qs
- source H5AD: objects/HCC_TA_8datasets_merged_review_v1.h5ad

Create new HGNC-standardised analysis objects without overwriting the sources:
- objects/HCC_TA_8datasets_merged_review_HGNC_v1.qs
- objects/HCC_TA_8datasets_merged_review_HGNC_v1.h5ad

Gene mapping uses current HGNC-approved symbols, with unique Ensembl/Entrez/previous-symbol/alias resolution and exact count summation for multiple source features mapping to the same HGNC symbol.

The QS is transformed chunk-by-chunk. The H5AD is transformed directly from its CSR matrix in 10,000-cell blocks and written on disk, so the full 38 GB matrix is never loaded into memory.

The original 68,394-feature merged v1 files remain preserved as the raw audit layer. Unmapped/custom/non-human features remain there but are excluded from the HGNC analysis objects.

## Next command
Execute task_004c.

Equivalent:
bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

Task 005 remains paused pending review of Task 004c.