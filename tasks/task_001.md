# Task 001 — Public dataset audit and acquisition

## Status
PENDING

## Scientific question
Which samples from the seven candidate human HCC scRNA-seq datasets are genuinely usable for tumour–adjacent analyses, and which downloadable matrix representation should be used as the phase-1 input?

## Datasets
GSE282701, GSE242889, GSE326201, GSE149614, GSE299340, GSE290298, GSE202642.

## Required work
1. Validate accession-level and sample-level metadata from authoritative public records.
2. Build one row per downloadable scRNA library/sample in `config/sample_manifest.tsv`.
3. Record patient ID, tissue, pairing, etiology, MVI when available, platform, cell-selection/enrichment strategy, and data representation.
4. Determine availability of raw-count matrix, filtered-count matrix and FASTQ/SRA data.
5. Prefer public count matrices for phase 1; do not bulk-download all FASTQs.
6. Create resumable server download scripts under `scripts/shell/`.
7. Download eligible phase-1 matrices to `raw_data/<dataset>/` on the server.
8. Generate SHA256 checksums and server-side download logs.
9. Mark inaccessible sources `DOWNLOAD_BLOCKED` and give an exact local-download/upload plan rather than silently substituting data.
10. Inspect downloaded files sufficiently to verify dimensions/file integrity without performing Task 002 preprocessing.

## Required outputs
- config/sample_manifest.tsv
- config/download_manifest.tsv
- results/task001_dataset_inventory.csv
- scripts/shell/task001_download_public_matrices.sh
- reports/task_001_report.md
- updated PROJECT_STATUS.md

Server-only outputs:
- raw_data/<dataset>/
- logs/task001_*
- checksums/task001_download_checksums.sha256 or equivalent documented location

## Required QC
Report per dataset:
- number of patients;
- tumour samples;
- adjacent/non-tumour samples;
- confirmed paired patients;
- matrix format and approximate size;
- raw/filtered/normalized status;
- FASTQ availability;
- enrichment/FACS/remixing bias;
- suitability for whole-tissue neutrophil abundance analysis;
- metadata ambiguities.

## Do not
- Do not begin QC filtering.
- Do not merge datasets.
- Do not annotate cells.
- Do not download all FASTQs unless explicitly required for a blocked matrix.
- Do not treat provisional sample counts in prior discussion as verified.

## Completion criterion
All seven datasets have a traceable sample inventory and each has either a verified server-side phase-1 matrix download or an explicit DOWNLOAD_BLOCKED plan. Stop after reporting.
