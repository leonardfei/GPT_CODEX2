# Task 002 report — sample-level neutrophil-friendly QC

**Status:** COMPLETED — HOLD FOR REVIEW

**Run date:** 2026-09-24

This report was regenerated after syncing local `origin/main` at commit `244a322`. The server-side Python QC, review PDF, and all 68 Seurat RDS objects were rerun; the rerun reproduced the committed QC table values exactly.

## Scope and synchronization

Before execution, the local repository was fast-forwarded from `origin/main` (`30db0f7` to `dade84f`). Task 002 was then run on the configured server under:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`

The approved phase-1 cohort was processed exactly as specified:

- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

GSE202642 and GSE290298 were not processed. No annotation or cross-dataset integration was performed.

## Input structure validation

All 68 planned samples passed the input structure audit.

- **GSE282701:** 12 per-sample gzipped 10x-style matrices.
- **GSE242889:** 10 nested POSIX tar archives containing matrix, feature/gene, and barcode files. Feature IDs, feature names, and cell barcodes were recovered for all 10 samples.
- **GSE326201:** 18 Cell Ranger filtered feature-barcode H5 files.
- **GSE149614:** one combined gene-by-cell count table plus metadata. Only primary Tumor and Adjacent/NTL samples were retained; PVTT and LN samples were excluded. Metadata labels recorded as `Normal` were mapped to the approved Adjacent tissue role.
- **GSE299340:** 10 gzipped 10x-style matrices.

The detailed structure results are in `results/task002_input_structure_audit.csv`.

## QC implementation

Each sample was processed independently. The workflow:

1. Loaded raw count matrices while preserving gene identifiers, gene names, and cell barcodes.
2. Calculated `nCount_RNA`, `nFeature_RNA`, and mitochondrial percentage per cell.
3. Applied sample-aware thresholds derived from the sample distributions and neutrophil-candidate distributions; no universal high minimum was imposed.
4. Used a conservative candidate neutrophil screen for retention auditing, based on core markers `FCGR3B`, `CSF3R`, `CXCR2`, and `FPR1`, with supporting markers including `S100A8`, `S100A9`, `CTSG`, `ELANE`, `MPO`, `LYZ`, `MCEMP1`, `FCAR`, `FFAR2`, `NAMPT`, and `CXCR4`.
5. Evaluated potential doublet risk using a high `nCount_RNA` plus high `nFeature_RNA` proxy. This was retained as a review flag only; it was not used for hard removal because mature neutrophils can be unusually low-complexity and aggressive doublet filtering could remove true biology.
6. Wrote a per-sample retention audit and a review PDF with one QC page per sample.

The Python implementation is `scripts/python/task002_preprocess_qc.py`. Seurat object creation is in `scripts/R/task002_create_seurat_objects.R`.

## Retention summary

| Metric | Result |
|---|---:|
| Samples processed | 68 |
| Cells before QC | 447,030 |
| Cells after QC | 372,882 |
| Candidate neutrophils before QC | 71,150 |
| Candidate neutrophils after QC | 66,323 |
| Candidate neutrophil retention, minimum | 89.6% |
| Candidate neutrophil retention, median | 93.6% |
| Samples flagged for suspiciously low retention | 0/68 |
| Doublet-proxy review flags | 1,082 |

The lowest candidate retention was GSE326201 sample `GSM9625180`, with 60 of 67 candidates retained (89.6%). This did not meet the pre-specified suspicious-retention flag threshold of 80%.

By dataset:

| Dataset | Samples | Cells before | Cells after | Candidates before | Candidates after |
|---|---:|---:|---:|---:|---:|
| GSE149614 | 18 | 63,101 | 48,919 | 9,246 | 8,647 |
| GSE242889 | 10 | 60,340 | 48,239 | 17,901 | 16,418 |
| GSE282701 | 12 | 141,842 | 125,833 | 14,585 | 13,669 |
| GSE299340 | 10 | 81,801 | 72,756 | 17,145 | 16,083 |
| GSE326201 | 18 | 99,946 | 77,135 | 12,273 | 11,506 |

## Threshold review items

Observed sample-specific thresholds were:

- `min_nFeature_RNA`: 290–2,064; median 689.
- `min_nCount_RNA`: 557–7,825; median 1,576.5.
- `max_percent_mt`: 20–60%; median 20%.

Twenty-nine samples received an adaptive mitochondrial ceiling above 20%, and one reached the 60% ceiling. These settings were intentionally conservative for mature neutrophil retention and should be reviewed before annotation/integration. The complete values and threshold basis are in `results/task002_qc_thresholds_by_sample.csv`.

## Seurat objects and validation

The server contains 68 per-sample Seurat RDS objects under:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task002_seurat/<dataset>/<sample_id>_qc.rds`

Validation passed for all 68 objects:

- 68/68 have a non-negative raw `counts` layer.
- 68/68 have exact cell counts matching the Python retention audit.
- Spot checks across GSE149614, GSE242889, and GSE326201 loaded successfully and exposed the expected `counts` layer.
- Raw unnormalized counts for QC-passing cells are preserved in the RDS objects. Full pre-QC cell-level metrics and staged matrices remain on the server.

The validation table is `results/task002_rds_validation.csv`. Large matrices and RDS objects are intentionally not committed to GitHub.

## Review artifacts

- `results/task002_neutrophil_retention_audit.csv`
- `results/task002_qc_thresholds_by_sample.csv`
- `results/task002_input_structure_audit.csv`
- `results/task002_rds_validation.csv`
- `results/task002_run_metadata.json`
- `figures/task002_qc_review.pdf`

## Hold point

Task 002 is complete and ready for Web GPT/user review. The next authorized step is to review the per-sample retention audit, adaptive thresholds, and QC PDF. Task 003 broad annotation and Task 004 integration must remain paused until that review explicitly approves proceeding.
