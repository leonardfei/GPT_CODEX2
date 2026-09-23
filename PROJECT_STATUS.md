# Project Status

## Current task
Task 001 — Public dataset audit and acquisition — PARTIAL

## Last completed task
None. Task 001 metadata audit is complete, but server-side acquisition is blocked.

## Repository status
GPT_CODEX2 control files are present locally. A local Git repository was initialized and committed at `407280a`; push to GitHub remains pending explicit authorization.

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

## Task 001 blocker
No configured compute server was reachable from this run, so all seven dataset downloads are explicitly `DOWNLOAD_BLOCKED`. Do not begin Task 002 until the matrices are downloaded and server-side integrity checks pass.

## Next execution command
After server access is restored, run `bash scripts/shell/task001_download_public_matrices.sh` on the compute server, then review Task 001 before executing any later task.

## Last update
2026-09-23
