# Project Status

## Current task
Task 002 — corrected neutrophil-preserving QC — READY TO RUN

## Last completed task
Task 002 initial QC — COMPLETED but superseded for downstream use pending correction.

## Repository status
GPT_CODEX2 control files are synchronized to `origin/main`; no raw matrices or large data are committed.

## Phase-1 analysis cohort
Five datasets remain approved:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Excluded for now:
- GSE202642
- GSE290298

## Task 002 review decision
The initial QC run successfully processed 68 samples and produced valid Seurat objects, but its QC threshold selection is considered circular for neutrophil-retention validation because candidate-neutrophil distributions directly determined the lower feature/count thresholds and mitochondrial ceilings.

Observed initial thresholds included:
- min_nFeature_RNA up to 2,064;
- min_nCount_RNA up to 7,825;
- max_percent_mt up to 60%.

Therefore the initial Task 002 objects must not be used for Task 003/004 until the corrected QC comparison is reviewed.

## Corrected Task 002 rule
QC thresholds must be selected independently from neutrophil identity:
- min_nFeature_RNA: global-sample 1st percentile, bounded to 100–300;
- min_nCount_RNA: global-sample 1st percentile, bounded to 200–500;
- max_percent_mt: global-sample 98th percentile, bounded to 20–30%.

High-confidence and broad granulocyte candidate definitions are used only for audit and must not determine thresholds.

## Preservation
Initial Task 002 report is archived as:
`reports/task_002_initial_report.md`

Initial server RDS objects under `objects/task002_seurat/` must remain unchanged.

Corrected objects must be written under:
`objects/task002_corrected_seurat/`

## Pending tasks
1. Execute corrected Task 002 and compare with initial QC
2. Web GPT/user review corrected QC
3. Task 003 — broad annotation and neutrophil confirmation
4. Task 004 — five-dataset Seurat integration

## Next execution command
`Execute task_002.`

## Last update
2026-09-24
