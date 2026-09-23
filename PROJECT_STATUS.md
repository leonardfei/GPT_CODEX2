# Project Status

## Current task
Task 001 — Public dataset audit and acquisition — COMPLETED

## Last completed task
Task 001 — public dataset audit, server-side acquisition, and integrity verification.

## Repository status
GPT_CODEX2 control files are present locally and synchronized to `origin/main`; no raw matrices or large data were committed.

## Scientific scope
Seven public human HCC scRNA-seq accessions with tumour/adjacent tissue are scheduled for validation and download.

## Important constraints
- Server does not rely on GitHub/general internet connectivity.
- Large data remain on the server and are not committed.
- Mature neutrophil retention is a mandatory preprocessing QC endpoint.
- Cohort counts/pairing remain provisional until Task 001 validation.

## Pending tasks
1. Task 001 — dataset audit and acquisition
2. Task 002 — sample-level neutrophil-friendly QC
3. Task 003 — broad annotation and neutrophil confirmation
4. Task 004 — Seurat v5 integration
5. Task 005 — optional raw-data reprocessing for neutrophil recovery

## Task 001 outputs
- `config/sample_manifest.tsv` — 92 GEO sample/library rows.
- `config/download_manifest.tsv` — 10 public matrix/metadata files with resumable server targets.
- `results/task001_dataset_inventory.csv` — seven dataset-level inventory/QC rows.
- `scripts/python/task001_build_manifests.py` — deterministic SOFT parser and manifest builder.
- `scripts/shell/task001_download_public_matrices.sh` — resumable/retryable server downloader with logs/checksums.
- `reports/task_001_report.md` — audit report.

## Task 001 acquisition record
- Nine required phase-1 files for six count-based datasets are present under `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/raw_data/`.
- Server-side `sha256sum -c` passed for all 9/9 required files; TAR/GZIP validation also passed.
- GSE290298 remains explicitly `DOWNLOAD_BLOCKED` by design because GEO exposes only an optional normalized matrix, not a preferred raw/filtered count input.
- Raw matrices, server logs, and the checksum manifest remain server-side and are not committed.

## Next execution command
Review the completed Task 001 inventory and scientific exclusions before executing Task 002.

## Last update
2026-09-24
