# Task 004c Report — In-place HGNC Harmonisation

## Status
COMPLETED

## Inputs actually used
- Server-side merged Task 004b QS: `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs` (1,490,852 cells; 68,394 source features).
- Server-side merged Task 004b H5AD: `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad` (1,490,852 observations; 68,394 variables; CSR `X`).
- Eight cohort source objects referenced by `results/task004b_merge_review_cohort_summary.csv`: GSE282701, GSE242889, GSE326201, GSE149614, GSE299340, CRA002308, nature_xue, and in_house. These were read only to calculate feature presence; no cohort object was rewritten.
- Official HGNC complete set, retrieved 2026-09-27 from [HGNC](https://www.genenames.org/download/). Download URL: `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`; SHA256 `3e5da5b757afce6cae333d3e6c97a59ff3aadfb13da1c70f89d9a8f252b9e7ec`.

## Data inspection
- Confirmed the expected 8 datasets, 194 samples, 132 patients and 1,490,852 cells before processing.
- Tissue counts: Tumor 1,039,293; Adjacent 451,559.
- Required dataset, sample, patient, tissue, project sample/patient IDs and source `nCount_RNA`/`nFeature_RNA` were complete. H5AD core observation metadata had zero missing values. Optional/source-specific fields retain their original, expected missingness; details are in `results/task004c_input_metadata_missingness.csv`.
- Cell IDs were unique and in the same order across the merged objects. Both replacements retain all 1,490,852 cells.

## Implementation
- Applied exact-case HGNC mapping in the specified order: current approved symbol, version-stripped Ensembl gene ID, Entrez ID, unique previous symbol, then unique alias. Ambiguous mappings were not forced. Multiple source rows assigned to one symbol were summed in raw counts.
- For QS, transformed the custom sparse counts layer chunk-by-chunk, preserved cell order and metadata, and added source/HGNC count and feature-detection metrics plus the harmonisation label.
- For H5AD, transformed the existing CSR matrix in 10,000-cell blocks, preserved observation-side structures/order, wrote int64 `indptr`, and cross-validated the retained count total against QS.
- Wrote temporary QS/H5AD outputs, validated both, then atomically replaced the existing merged v1 paths under the D014-approved in-place replacement. No `*_HGNC_v1` duplicate was retained. The validated Task 004b source checkpoint remains on the server for recovery.
- No cell filtering, sample filtering, or biological reannotation was performed.

## Parameters and software
- R 4.3.3; Seurat 5.3.0; SeuratObject 5.2.0; qs 0.27.3; data.table 1.17.8; Matrix 1.6.5.
- Python 3.14.7; anndata 0.13.4; SciPy 1.18.1; NumPy 2.5.3; pandas 3.0.6; h5py 3.16.0.
- Scripts: `scripts/R/task004c_harmonize_features.R`, `scripts/python/task004c_harmonize_merged_h5ad.py`, `scripts/bash/task004c_harmonize_features.sh`, and `scripts/R/task004c_neutrophil_retention_audit.R`.
- Main command: `bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`; final promotion used `--promote-only` after both temporary outputs passed validation. H5AD block size was 10,000 cells; QS sparse processing used 8 threads. No random operation or seed was used.

## QC
- Final QS: VALIDATED; 1,490,852 cells × 38,025 unique approved symbols; 0 duplicate symbols; 4,562,511,399 bytes.
- Final H5AD: VALIDATED; 1,490,852 × 38,025; unique observation and variable names; CSR `X`; int64 `indptr`; 2,846,908,974 nonzero entries; 7,879,382,527 bytes. Cross-validation against QS passed.
- Source raw-count sum: 9,960,157,768. Retained HGNC count sum: 9,890,904,647 (99.3047%); 69,253,121 counts (0.6953%) belonged to unmapped/ambiguous features removed from the merged analysis objects. This is feature-level count removal, not cell loss. Counts mapping to retained symbols were summed exactly.
- Mapping status across 68,394 source feature rows: 36,620 exact approved symbols; 1,910 unique previous symbols; 427 unique aliases; 5 ambiguous previous symbols; 5 ambiguous aliases; 29,427 unmapped. No successful Ensembl-ID or Entrez-ID mapping rows occurred in this source feature union.
- Eight-cohort mapped feature-row / unique-symbol counts: GSE282701 24,208 / 24,098; GSE242889 37,538 / 37,436; GSE326201 24,561 / 24,268; GSE149614 19,795 / 19,743; GSE299340 24,208 / 24,098; CRA002308 24,208 / 24,098; nature_xue 19,516 / 19,513; in_house 24,208 / 24,098.
- Shared sets: 17,680 HGNC symbols present in all 8 cohorts; 21,169 present in at least 7 of 8.
- The final paths were reopened after promotion and their dimensions, identifiers, sparse structure and counts were checked. Temporary outputs and replacement backups were absent after successful promotion.

## Main observed results
The merged objects now use 38,025 unique HGNC-approved symbols at the original v1 paths, with all cells and their ordering preserved. Feature rows that could not be mapped unambiguously were excluded as authorized; the retained Task 004b checkpoint and upstream cohort objects remain the recovery sources.

## Neutrophil-specific QC
- The sample-aware audit uses the preliminary Task 004 broad label `project_broad_celltype == "neutrophil"`; this is inherited annotation, not a new identity call. Task 004b annotations remain preliminary and require review.
- Candidate neutrophils before/after: 65,546 / 65,546; retention 100%. Other cells: 1,425,306. Task 004c removed zero cells, so there were no cell-removal reasons.
- `results/task004c_neutrophil_retention_audit.csv` reports, by dataset/sample/patient/tissue and candidate-vs-other group, cell counts and q25/median/q75 for source and harmonised `nCount_RNA`, source and harmonised `nFeature_RNA`, and mitochondrial percentage. Since no cells were filtered, the candidate retention is 1 in every sample with candidates.
- This task does not estimate tissue abundance and does not support pooled-cell abundance comparisons.

## Outputs
- Required review tables: `results/task004c_merged_feature_mapping.csv`, `results/task004c_hgnc_var.csv`, `results/task004c_mapping_status_counts.csv`, `results/task004c_feature_audit_by_cohort.csv`, `results/task004c_hgnc_feature_presence.csv`, `results/task004c_shared_hgnc_features_8of8.txt`, `results/task004c_shared_hgnc_features_7plus.txt`, `results/task004c_qs_validation.csv`, and `results/task004c_h5ad_validation.json`.
- Additional audit/provenance: `results/task004c_input_metadata_missingness.csv`, `results/task004c_neutrophil_retention_audit.csv`, and `results/task004c_hgnc_reference_provenance.txt`.
- Final server objects: `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.qs` and `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.h5ad`.

## Unexpected findings / deviations
- Ten source feature rows had ambiguous previous/alias mappings and were left unmapped. The 29,427 unmapped features include custom/non-human/otherwise unmatched identifiers; no case-insensitive rescue was attempted.
- Several implementation/validation issues were caught and corrected before promotion (sparse chunk-slot preservation, mapping-table field handling, and H5AD categorical-index validation). Failed attempts did not modify the original merged paths; promotion occurred only after both outputs validated. No scientific mapping rule was changed.

## Scientific interpretation candidates
None. This is identifier harmonisation only; no biological cell-type, differential-expression or abundance conclusion is made.

## Questions for Web GPT
Review the mapping and QC outputs before authorising Task 005. Task 005 remains paused and was not run.

## Git
Pending this report/status update; commit and push to `origin/main` after final local QC.
