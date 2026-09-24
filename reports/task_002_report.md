# Task 002 report — corrected neutrophil-preserving per-sample QC

**Status:** COMPLETED — HOLD FOR REVIEW

**Run date:** 2026-09-24

## Scope and provenance

Before the corrected run, the local repository was synchronized with GitHub and the latest Task 002 correction decision was merged. The corrected run processed the approved phase-1 cohort only:

- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

GSE202642 and GSE290298 remained excluded. No annotation or cross-dataset integration was performed.

The initial QC report is preserved at `reports/task_002_initial_report.md`. Initial objects remain under:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task002_seurat/`

Corrected objects were written separately under:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task002_corrected_seurat/`

## Corrected QC rule

The initial run used neutrophil-candidate distributions to select thresholds, which made candidate-retention validation circular. The corrected run selected thresholds from all source cells before applying any marker-based identity rule:

```text
min_nFeature_RNA = max(100, min(300, floor(global nFeature_RNA 1st percentile)))
min_nCount_RNA   = max(200, min(500, floor(global nCount_RNA 1st percentile)))
max_percent_mt   = min(30, max(20, ceil(global percent.mt 98th percentile)))
```

Neutrophil and broad granulocyte definitions were used only for retention auditing. High-count/high-feature doublet proxies were retained as review flags and were not hard-filtered.

## Input and sample structure

All 68 planned samples passed the structure audit. GSE242889 feature IDs, feature names, and barcodes were recovered for all samples. GSE149614 retained only primary Tumor and Adjacent/NTL samples; PVTT and LN samples were excluded. The full corrected structure table is `results/task002_corrected_input_structure_audit.csv`.

## Initial versus corrected totals

| Metric | Initial QC | Corrected QC |
|---|---:|---:|
| Source cells | 447,030 | 447,030 |
| Cells passing QC | 372,882 | 422,856 |
| High-confidence granulocyte/neutrophil candidates retained | 61,337 | 64,838 |
| Broad granulocyte-like candidates retained | 65,335 | 69,128 |
| Cells rescued by corrected QC | — | 53,023 |
| Cells retained only by initial QC | — | 3,049 |

The corrected run retained more cells overall by design, while remaining subject to bounded mitochondrial thresholds and explicit review flags.

## Corrected thresholds and retention

Observed corrected threshold ranges across 68 samples:

- `min_nFeature_RNA`: 100–300; median 256.5.
- `min_nCount_RNA`: 500 for all samples because the global 1st-percentile values reached the upper bound.
- `max_percent_mt`: 20–30%; median 30.

Corrected high-confidence candidate retention was 92.6% at minimum and 99.5% at the median. Corrected broad-granulocyte-like retention was 90.7% at minimum and 99.5% at the median.

The corrected audit generated 16 mandatory sample review flags:

- 15 samples retained more than 20 percentage points more total cells than the initial QC.
- 1 sample had more than half of its rescued broad-granulocyte-like cells without a core granulocyte marker.
- No sample had corrected high-confidence retention below 90%.
- No sample had rescued-cell median mitochondrial percentage above 25%.
- No sample exceeded its source-cell count or had a cell-identity mismatch.

These are review flags, not automatic rejection decisions. The sample-level details are in `results/task002_corrected_neutrophil_retention_audit.csv` and `results/task002_initial_vs_corrected_qc.csv`.

By dataset:

| Dataset | Samples | Source cells | Initial pass | Corrected pass | Rescued | Corrected high-confidence retained | Corrected broad retained |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE149614 | 18 | 63,101 | 48,919 | 63,101 | 14,182 | 7,282 | 8,895 |
| GSE242889 | 10 | 60,340 | 48,239 | 52,954 | 7,691 | 15,141 | 16,530 |
| GSE282701 | 12 | 141,842 | 125,833 | 137,518 | 11,685 | 13,870 | 14,395 |
| GSE299340 | 10 | 81,801 | 72,756 | 76,319 | 3,628 | 16,703 | 17,054 |
| GSE326201 | 18 | 99,946 | 77,135 | 92,964 | 15,837 | 11,842 | 12,254 |

The rescued-cell distributions and marker-context fractions are in `results/task002_rescued_cells_summary.csv`.

## Corrected Seurat objects and validation

The server contains 68 corrected per-sample Seurat RDS objects under:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task002_corrected_seurat/<dataset>/<sample_id>_qc.rds`

Validation passed for all corrected objects:

- 68/68 corrected RDS files exist.
- 68/68 have a non-negative raw `counts` layer.
- 68/68 corrected RDS cell counts exactly match the corrected QC audit.
- Spot checks across GSE149614, GSE242889, and GSE326201 loaded successfully with the expected `counts` layer.
- Raw unnormalized counts for corrected-QC cells and sample metadata are preserved.

The validation table is `results/task002_corrected_rds_validation.csv`. Large matrices and RDS objects remain server-side and are not committed to GitHub.

## Review artifacts

- `results/task002_corrected_qc_thresholds_by_sample.csv`
- `results/task002_corrected_neutrophil_retention_audit.csv`
- `results/task002_initial_vs_corrected_qc.csv`
- `results/task002_rescued_cells_summary.csv`
- `results/task002_corrected_rds_validation.csv`
- `results/task002_corrected_input_structure_audit.csv`
- `results/task002_corrected_run_metadata.json`
- `figures/task002_corrected_qc_review.pdf`

## Hold point

Corrected Task 002 is complete, but the 16 flagged samples and the initial-versus-corrected comparison require Web GPT/user review. Do not execute Task 003 annotation or Task 004 integration until the corrected QC is explicitly approved.
