# Project Status

## Current task
Task 004c — in-place HGNC harmonisation of merged QS/H5AD — COMPLETED

Task 004b is COMPLETED:
- 1,490,852 cells
- 68,394 source-union features before harmonisation
- 8 datasets, 194 samples, 132 patients

Task 004c is COMPLETED:
- Both existing merged v1 paths were replaced only after temporary QS and H5AD outputs passed validation.
- Final objects contain 1,490,852 cells and 38,025 unique HGNC-approved symbols; cell order and metadata were preserved.
- Raw-count retention after removal of unmapped/ambiguous features: 99.3047% (9,890,904,647 / 9,960,157,768).
- Shared HGNC features: 17,680 in all 8 cohorts; 21,169 in at least 7 cohorts.
- Candidate neutrophils (preliminary Task 004 label): 65,546 before/after, 100% retention; no cells were filtered.
- Review details: `reports/task_004c_report.md` and `results/task004c_*`.

Final server paths:
- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs`
- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`

Task 005 remains PAUSED pending review of Task 004c; it has not been executed.
