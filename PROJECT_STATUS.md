# Project Status

## Current task
Task 002 — sample-level neutrophil-friendly QC — COMPLETED; HOLD FOR REVIEW

## Last completed task
Task 002 — sample-level neutrophil-friendly QC on the five approved datasets — COMPLETED; review hold active.

## Repository status
GPT_CODEX2 control files are synchronized to `origin/main`; no raw matrices or large data are committed.

## Phase-1 analysis cohort
The user has approved a five-dataset phase-1 cohort:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

The following two datasets are explicitly excluded from Tasks 002–004 for now:
- GSE202642 — FACS/composition bias and unresolved patient pairing
- GSE290298 — normalized-only matrix in GEO

Their existing files/metadata must be retained unchanged but not processed further unless explicitly reactivated.

## Important constraints
- Server does not rely on GitHub/general internet connectivity.
- Large data remain on the server and are not committed.
- Mature neutrophil retention is a mandatory preprocessing QC endpoint.
- Patient/sample, not cell, is the biological replicate for abundance analyses.

## Task 001 acquisition record
- Nine required files for six count-based datasets are present under `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/raw_data/`.
- Server-side `sha256sum -c` passed for all 9/9 required files; TAR/GZIP validation passed.
- GSE290298 was intentionally not downloaded because GEO exposes only normalized expression.

## Task 002 execution record
- Local repository was updated from `origin/main` before execution; the run started from commit `dade84f`.
- Processed 68 samples from GSE282701, GSE242889, GSE326201, GSE149614, and GSE299340.
- GSE202642 and GSE290298 remained explicitly excluded.
- All 68 input structures passed validation, including recovery of feature IDs, feature names, and barcodes for all GSE242889 samples.
- Cells before/after QC: 447,030 / 372,882.
- Candidate neutrophils before/after QC: 71,150 / 66,323.
- Candidate neutrophil retention: 89.6% minimum, 93.6% median; 0/68 suspiciously low-retention flags.
- 68 per-sample Seurat RDS objects were created on the server with raw counts preserved for QC-passing cells. RDS validation passed for 68/68 objects.
- Review artifacts: `results/task002_neutrophil_retention_audit.csv`, `results/task002_qc_thresholds_by_sample.csv`, `results/task002_input_structure_audit.csv`, `results/task002_rds_validation.csv`, and `figures/task002_qc_review.pdf`.
- Annotation and integration are paused pending review of the retention audit, sample-specific thresholds, and QC plots.

## Pending tasks
1. Review Task 002 retention audit, adaptive thresholds, and QC plots; hold point is active
2. Task 003 — broad annotation and neutrophil confirmation on five approved datasets
3. Task 004 — integrate the five approved datasets into one Seurat object
4. Task 005 — optional raw-data reprocessing if neutrophil recovery is implausibly poor

## Next execution command
Await Web GPT/user review of Task 002 before executing Task 003.

## Last update
2026-09-24
