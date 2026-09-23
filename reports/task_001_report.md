# Task 001 — Public dataset audit and acquisition

## Status

COMPLETED WITH EXPLICIT EXCLUSION — accession/sample audit, server-side acquisition, and integrity verification completed. Six count-based datasets were acquired; the optional normalized-only GSE290298 file was intentionally not downloaded.

## Inputs actually used

- Project instructions: `AGENTS.md`, `PROJECT_CONTEXT.md`, `SCIENTIFIC_DECISIONS.md`, `PROJECT_STATUS.md`, `tasks/task_001.md` and `config/metadata_schema.tsv`.
- NCBI Entrez GEO `esearch`/`esummary` records for GSE282701, GSE242889, GSE326201, GSE149614, GSE299340, GSE290298 and GSE202642.
- GEO family SOFT metadata and GEO supplementary directory/filelist records retrieved on 2026-09-23 UTC.
- Representative matrix/file-header checks were performed without Task 002 preprocessing. No raw data were written to the repository.
- Required public matrices were staged outside the repository, uploaded to the configured server data disk, and verified in place. The temporary staging area was removed after transfer.
- NCBI SRA accession-term search and SRA summaries. The 18 returned runs all belong to GSE326201; no SRA hit was returned for the other six accession terms. GSE149614's GEO design notes EGA raw-data access under EGAS00001004468.

Authoritative accession records: [GSE282701](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282701), [GSE242889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242889), [GSE326201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326201), [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614), [GSE299340](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299340), [GSE290298](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290298), [GSE202642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202642).

## Data inspection

The audit produced 92 GEO sample rows. GSE290298 contains 10 GEO records, of which two are explicitly scTCR-seq-only; they are retained in the sample manifest but excluded from the phase-1 expression matrix. GSE149614 contains 21 records, including 18 tumour/adjacent target libraries and three non-target libraries (two PVTT and one metastatic lymph-node library).

| Dataset | Patients / tissue records | Confirmed paired patients | Phase-1 representation | Approx. compressed size | FASTQ/SRA | Abundance suitability |
|---|---:|---:|---|---:|---|---|
| GSE282701 | 6 / 6 T + 6 A | 6 | 12 per-sample Matrix Market matrices in RAW.tar; filtered counts | 856,524,800 B | No SRA hit | Yes, unselected paired design; MVI field absent |
| GSE242889 | 5 / 5 T + 5 A | 5 | 10 per-sample tar.gz archives; author-described raw count matrices | 1,376,245,760 B | No SRA hit | Yes, patient-level paired design; retain MVI strata |
| GSE326201 | 10 / 10 T + 8 A | 8 | 18 Cell Ranger filtered HDF5 matrices in RAW.tar | 244,039,680 B | Yes, 18 public SRA runs | Yes for the eight paired patients; two tumour-only patients |
| GSE149614 | 10 / 10 T + 8 A + 3 other | 8 | Combined raw count text matrix plus updated cell metadata | 165,349,783 B | EGA cited; no SRA hit | Conditional; exclude PVTT/lymph-node libraries |
| GSE299340 | 5 / 5 T + 5 A | 5 | 10 integer Matrix Market matrices with features/barcodes in RAW.tar | 843,857,920 B | No SRA hit | Yes with source-QC caveat |
| GSE290298 | 4 / 4 T + 4 A + 2 TCR-only | 4 | Normalized floating-point processed CSV only | 923,544,964 B | No SRA hit | Conditional/not preferred for raw-count phase 1 |
| GSE202642 | Patient pairing not verified / 7 T + 4 A | 0 confirmed | Combined integer Matrix Market, 36,601 × 115,732; features/barcodes supplied | 698,853,827 B matrix | No SRA hit | No for unbiased whole-tissue abundance: FACS-sorted |

### Sample-level metadata findings

- GSE282701: six paired patients, HBV-positive and HBV-negative groups; the accession titles encode tumour/adjacent and viral status. The Series design states no vascular invasion/metastatic features, but this was not converted into a sample-level MVI-negative value.
- GSE242889: five paired cases; tumour records identify three MVI-present and two MVI-absent cases. Adjacent records inherit the patient-level MVI assignment for pairing, with the source wording retained as an ambiguity note.
- GSE326201: ten patients, eight confirmed tumour/adjacent pairs, two tumour-only libraries; etiology is explicit as HBV, HCV or NBNC. The Series provides filtered `.h5` matrices.
- GSE149614: eight confirmed primary tumour/adjacent pairs; HCC07/HCC08 PVTT and HCC10 lymph-node libraries are marked `Other` and excluded from the tumour–adjacent phase-1 role. The updated metadata contains 71,915 public cell rows and provides patient, site, stage and HBV/HCV fields.
- GSE299340: five paired tumour/noncancerous patient pairs. GEO metadata explicitly records microvascular invasion as negative for the patients and records HBV/HCV status as a two-part field; the manifest normalizes this to HBV, HCV or NBNC while preserving the source accession.
- GSE290298: four paired gene-expression samples plus one N4 and one T4 scTCR-only record. The only Series-level matrix observed is normalized/processed, so it is not the default phase-1 raw-count input.
- GSE202642: seven tumour tissues and four adjacent liver tissues, all HBV-related. The protocol explicitly states FACS-sorted cells; patient-level pairing could not be verified from the GEO sample records and is marked unknown rather than inferred.

### Matrix integrity checks

The following checks were performed on public headers/file listings only:

- GSE282701 representative matrix: Matrix Market integer header, 36,601 genes × 18,329 barcodes.
- GSE299340 representative matrix: Cell Ranger 7.0.0 Matrix Market integer header, 36,601 genes × 7,388 barcodes.
- GSE202642 matrix: Matrix Market header reports 36,601 × 115,732 with 180,199,647 nonzeros; features and barcodes are separate gzip files.
- GSE290298 CSV header begins with sample-prefixed cell barcodes and floating-point expression values, supporting the `normalized` classification.
- GSE149614 count matrix is tab-delimited with cell-prefixed columns; GEO describes it as raw counts, with a separate normalized file.
- GSE326201 GEO filelist identifies 18 HDF5 filtered feature-barcode matrices.
- GSE242889 representative archive exposes a `matrix.mtx` member and a Matrix Market header; full server-side extraction and companion-file validation remain pending.

No QC filtering, merging, annotation, normalization, or cell-level transformation was performed.

## Implementation

1. Parsed the seven GEO family SOFT files with `scripts/python/task001_build_manifests.py`.
2. Normalized accession/sample titles into patient, tissue, pairing, etiology and MVI fields without inventing missing relationships.
3. Added the required sample-level manifest, file-level download manifest and dataset inventory.
4. Added a resumable/retryable parallel-range download script with per-file logging, SHA256 checksums and basic gzip/tar verification.
5. Uploaded the nine required phase-1 files to `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/` and verified them in place. The server checksum manifest passed `sha256sum -c` for all 9/9 files; temporary upload parts and the transfer test file were removed.

## Parameters and software

- Python 3.14.7; Bash 3.2.57; curl 7.88.1; bsdtar 3.5.3; gzip 403.100.6.
- NCBI/GEO records checked at 2026-09-23 UTC.
- Server download parameters: resumable HTTP Range requests with retry/backoff, 32 MiB chunks, and bounded parallelism; SHA256 via `sha256sum`; gzip validation via `gzip -t`; tar validation via `tar -tf`.
- No stochastic operation; no seed required.

## QC

- Accession inventory: 7/7 requested GEO Series resolved in NCBI GEO.
- Sample inventory: 92 rows generated; all requested accession/sample IDs are represented.
- Pairing: 6, 5, 8, 8, 5 and 4 confirmed paired patients for GSE282701, GSE242889, GSE326201, GSE149614, GSE299340 and GSE290298 respectively; GSE202642 remains unverified.
- Matrix availability: six datasets have a public raw/filtered count representation identified; GSE290298 is normalized-only in the observed Series supplementary record.
- Server acquisition: 6/7 dataset groups have their required phase-1 files downloaded and verified on the configured server; GSE290298 is explicitly `DOWNLOAD_BLOCKED` because its only observed matrix is optional normalized-only material.
- Checksums: `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/checksums/task001_download_checksums.sha256`; `sha256sum -c` passed for 9/9 required files.

## Neutrophil-specific QC

Not applicable to Task 001. No cell filtering or neutrophil candidate audit was performed. GSE202642's FACS selection bias is recorded because it affects downstream abundance interpretation.

## Outputs

- [config/sample_manifest.tsv](../config/sample_manifest.tsv)
- [config/download_manifest.tsv](../config/download_manifest.tsv)
- [results/task001_dataset_inventory.csv](../results/task001_dataset_inventory.csv)
- [scripts/python/task001_build_manifests.py](../scripts/python/task001_build_manifests.py)
- [scripts/shell/task001_download_public_matrices.sh](../scripts/shell/task001_download_public_matrices.sh)
- [PROJECT_STATUS.md](../PROJECT_STATUS.md)
- Server-only outputs produced: `raw_data/<dataset>/` for six count-based datasets, `logs/task001_*`, and `checksums/task001_download_checksums.sha256`.

## Unexpected findings / deviations

- The control workspace initially contained project files but no `.git` directory or usable remote refs. A local Git repository was initialized and the Task 001 control-layer files were committed; external publication remains pending explicit authorization.
- Direct server-side NCBI transfers were heavily throttled, so the verified public files were staged outside the repository and uploaded over the authenticated SSH connection; no raw data were committed.
- GSE242889's representative archive exposed `matrix.mtx` but did not expose companion feature/barcode files in the observed archive listing; this requires server-side confirmation before phase-1 ingestion.
- GSE290298 is not a raw/filtered count source in the observed GEO Series file; it is retained as optional normalized-only material.

## Scientific interpretation candidates

- The six unselected paired datasets are candidates for patient-level tumour–adjacent analyses, subject to source-QC and matrix validation.
- GSE202642 should be treated as composition-biased and restricted to cell-state analyses unless a separate scientific decision justifies FACS-derived abundance use.
- GSE149614's non-target PVTT/lymph-node libraries should not be pooled into tumour–adjacent abundance denominators.

These are implementation handoffs, not new scientific decisions; no existing decision file was changed.

## Questions for Web GPT

1. Confirm whether GSE290298 normalized-only material should be downloaded for exploratory cell-state work after the six count-based datasets.
2. Confirm whether GSE242889 should remain phase-1 eligible if its archive lacks feature/barcode companions for some per-sample members.

## Git

- Branch/commit: `main`; final Task 001 control-layer update is ready to commit and push.
- Commit: no raw matrices or large data included.
- Push: will be completed to `origin/main` after final local validation, under the user's explicit authorization.
