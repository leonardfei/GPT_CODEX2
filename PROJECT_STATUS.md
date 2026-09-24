# Project Status

## Current task
Task 002 — corrected neutrophil-preserving QC — COMPLETED; HOLD FOR REVIEW

## Last completed task
Task 002 corrected QC — COMPLETED; corrected-vs-initial comparison is ready for review.

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

## Corrected Task 002 execution record
- Corrected QC ran after the latest GitHub Task 002 decision was synchronized and merged locally.
- Processed 68 samples from GSE282701, GSE242889, GSE326201, GSE149614, and GSE299340; GSE202642 and GSE290298 remained excluded.
- Thresholds were independent of neutrophil identity: all-source-cell 1st percentiles bounded to 100–300 features and 200–500 counts, and all-source-cell 98th percentile mitochondrial fraction bounded to 20–30%.
- Source/initial/corrected cells: 447,030 / 372,882 / 422,856.
- Corrected rescued cells: 53,023; initial-only cells: 3,049.
- Corrected high-confidence candidate retention: 92.6% minimum, 99.5% median; corrected broad-granulocyte retention: 90.7% minimum, 99.5% median.
- 16/68 samples carry mandatory review flags: 15 for >20 percentage-point total-cell retention gain and 1 for more than half of rescued broad-granulocyte-like cells lacking a core marker.
- 68 corrected Seurat RDS objects were created under `objects/task002_corrected_seurat/`; raw counts and sample metadata were validated for 68/68 objects. Initial objects remain under `objects/task002_seurat/` for comparison.
- Corrected artifacts: `results/task002_corrected_qc_thresholds_by_sample.csv`, `results/task002_corrected_neutrophil_retention_audit.csv`, `results/task002_initial_vs_corrected_qc.csv`, `results/task002_rescued_cells_summary.csv`, `results/task002_corrected_rds_validation.csv`, and `figures/task002_corrected_qc_review.pdf`.

## Pending tasks
1. Web GPT/user review corrected QC and the 16 flagged samples
2. Task 003 — broad annotation and neutrophil confirmation
3. Task 004 — five-dataset Seurat integration

## Next execution command
Await Web GPT/user review of corrected Task 002 before executing Task 003.

## Last update
2026-09-24
