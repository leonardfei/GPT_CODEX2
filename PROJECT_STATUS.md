# Project Status

## Current task
Task 002 — sample-level neutrophil-friendly QC — READY TO RUN

## Last completed task
Task 001 — public dataset audit, server-side acquisition, and integrity verification — COMPLETED.

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

## Pending tasks
1. Task 002 — neutrophil-friendly QC on five approved datasets
2. Task 003 — broad annotation and neutrophil confirmation on five approved datasets
3. Task 004 — integrate the five approved datasets into one Seurat object
4. Task 005 — optional raw-data reprocessing if neutrophil recovery is implausibly poor

## Next execution command
`Execute task_002.`

## Last update
2026-09-24
