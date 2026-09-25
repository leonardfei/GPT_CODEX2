# Project Status

## Current task
Task 004 broad annotation and neutrophil confirmation — COMPLETED

## Last completed task
Task 004 harmonized broad annotation and neutrophil confirmation — COMPLETED; corrected cell set and annotation layer are frozen for review before Task 005.

Task 002 corrected neutrophil-preserving QC for the original five cohorts remains accepted for downstream use.

## Correction approved by Web GPT/user
Two metadata issues identified during Task 003 review must be corrected without changing cell sets or rerunning QC:

1. **nature_xue patient/pairing**
   - 92 HCC Tumor/AL sample records are currently stored.
   - Sample suffixes `_HCC`, `_HCC_N`, and `_HCC_IM<number>` must map to the same base A-number patient.
   - Expected: 79 unique patients, 10 Tumor–Adjacent paired patients, 69 Tumor-only patients.
   - A074 and A119 contain multiple Tumor/IM sample records; these remain separate sample IDs but share patient IDs.
   - A119_HCC, A119_HCC_IM1, A119_HCC_IM2 and A119_HCC_N must all map to `nature_xue_A119`.

2. **CRA002308 abundance eligibility**
   - Change from `NO` to `CONDITIONAL`.
   - Live nucleated-cell flow sorting after doublet exclusion can alter composition, but no lineage-specific immune enrichment is documented.
   - CRA002308 may be used for paired within-cohort Tumor–Adjacent abundance sensitivity analyses.
   - It must not be pooled as an unbiased absolute whole-tissue fraction dataset across cohorts.

## Cell-set preservation
The correction is metadata-only:
- nature_xue must remain exactly 675,539 cells.
- CRA002308 must remain exactly 158,741 corrected-QC cells.
- No Task 002/003 QC threshold or cell pass/fail decision may change.
- Original raw/source objects remain unchanged.

## Original five cohorts
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Authoritative inputs:
`objects/task002_corrected_seurat/`

## Extension cohorts
- CRA002308 — 7 paired patients; 158,741 cells; atlas YES; abundance CONDITIONAL after correction.
- nature_xue — 675,539 cells; expected 79 patients with 10 paired after correction; atlas YES; abundance CONDITIONAL.
- in_house — 10 paired patients; 233,716 cells; atlas YES; abundance CONDITIONAL pending sampling provenance.

## On hold
- GSE202642
- GSE290298

## Planned atlas
Eight cohorts:
GSE282701, GSE242889, GSE326201, GSE149614, GSE299340, CRA002308, nature_xue, in_house.

Final target:
`objects/HCC_TA_8datasets_integrated_v1.rds`

## Pending tasks
1. Web GPT/user review Task 004 annotation, neutrophil confirmation, rescued-cell audit, and abundance-role summaries
2. Task 005 — final eight-cohort Seurat v5 integration, only after review approval

## Next execution command
`Execute task_005.` after review and approval of the Task 004 hold point.

## Last update
2026-09-25
