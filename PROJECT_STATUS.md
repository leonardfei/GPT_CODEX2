# Project Status

## Current task
Task 003 — prepare three newly uploaded cohorts — READY TO RUN

## Last completed task
Task 002 corrected neutrophil-preserving QC for the original five cohorts — COMPLETED and accepted for downstream use.

## Repository status
GPT_CODEX2 control files are synchronized to `origin/main`; raw matrices and large RDS objects remain server-side.

## Original five cohorts
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Authoritative inputs:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task002_corrected_seurat/`

Corrected Task 002 retained 422,856 of 447,030 source cells. The corrected objects supersede the initial Task 002 objects for downstream work. Keep `qc_status`/rescued-cell provenance available for Task 004 annotation QC.

## Newly uploaded cohorts
The user added the following under:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/raw_data/`

1. `CRA002308`
   - expected 7 HCC Tumor and 7 matched adjacent/normal liver samples;
   - file format and cell-selection/composition bias must be audited before abundance eligibility is assigned.

2. `nature_xue`
   - author-processed Seurat object from Xue et al., Nature 2022;
   - retain human HCC only;
   - retain Tumor and adjacent liver (AL) only;
   - preserve author annotations and original assay/layer provenance;
   - do not reapply Task 002 QC blindly to this author-processed object.

3. `in_house`
   - in-house HCC cohort;
   - YJCA = Tumor;
   - YJP = Adjacent;
   - sample pairing and matrix structure must be verified from the uploaded files.

## On hold
- GSE202642
- GSE290298

## Planned atlas
The target integration expands from five to eight cohorts:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340
- CRA002308
- nature_xue
- in_house

Final target:
`objects/HCC_TA_8datasets_integrated_v1.rds`

## Pending tasks
1. Task 003 — audit/subset/QC-prepare CRA002308, nature_xue and in_house
2. Web GPT/user review extension-cohort preparation
3. Task 004 — harmonized broad annotation and neutrophil confirmation across the eight-cohort atlas
4. Task 005 — final eight-cohort Seurat v5 integration

## Next execution command
`Execute task_003.`

## Last update
2026-09-24
