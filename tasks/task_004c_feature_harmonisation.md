# Task 004c — HGNC harmonisation directly in merged QS and H5AD

## Status
READY TO EXECUTE

## Goal
Standardise gene identifiers directly in the validated 1,490,852-cell merged review objects, while preserving the original merged v1 QS/H5AD unchanged for audit and rollback.

## Source objects
- objects/HCC_TA_8datasets_merged_review_v1.qs
- objects/HCC_TA_8datasets_merged_review_v1.h5ad

## New outputs
- objects/HCC_TA_8datasets_merged_review_HGNC_v1.qs
- objects/HCC_TA_8datasets_merged_review_HGNC_v1.h5ad

Do not overwrite the source v1 objects.

## HGNC reference and mapping
Use the current official HGNC complete set and record retrieval date and SHA256.

Mapping priority:
1. exact current HGNC-approved symbol;
2. version-stripped Ensembl gene ID;
3. Entrez/NCBI Gene ID;
4. unique previous HGNC symbol;
5. unique HGNC alias.

Ambiguous previous symbols/aliases/IDs are not force-mapped. Case-insensitive matching is not used.

If multiple original features map to the same approved HGNC symbol, their raw counts are summed exactly.

Unmapped/custom/non-human features are excluded from the HGNC analysis object but remain preserved in the original merged_review_v1 QS/H5AD.

## QS implementation
The merged QS contains a task004b_chunked_counts layer. Apply the feature mapping once to the 68,394-feature union, then transform each existing sparse chunk independently using a sparse aggregation matrix. Preserve all 1,490,852 cells and metadata.

Add:
- source_nCount_RNA_before_feature_harmonisation
- source_nFeature_RNA_before_feature_harmonisation
- hgnc_nCount_RNA
- hgnc_nFeature_RNA
- feature_harmonisation = HGNC_approved_symbol

Validate by qread after writing the new QS.

## H5AD implementation
Transform the validated merged H5AD directly, not via eight cohort H5AD files.

Read the source CSR matrix in cell blocks (default 10,000 cells), remap/collapse feature columns with a sparse old-feature-to-HGNC matrix, and append the harmonised CSR arrays directly to a new H5AD on disk.

Requirements:
- preserve obs and cell order exactly;
- replace var with unique HGNC-approved symbols plus HGNC metadata;
- use CSR X with int64 indptr;
- avoid loading the full 38 GB matrix into memory;
- cross-validate total source and retained HGNC counts against the QS transformation.

## Shared-gene audit for Task 005
Although the harmonisation is performed directly in the merged files, inspect the eight cohort source feature lists only to calculate:
- per-cohort mapped HGNC gene counts;
- strict 8/8 shared HGNC set;
- >=7/8 shared HGNC set.

No separate harmonised cohort objects are written.

## Required outputs
- results/task004c_merged_feature_mapping.csv
- results/task004c_hgnc_var.csv
- results/task004c_mapping_status_counts.csv
- results/task004c_feature_audit_by_cohort.csv
- results/task004c_hgnc_feature_presence.csv
- results/task004c_shared_hgnc_features_8of8.txt
- results/task004c_shared_hgnc_features_7plus.txt
- results/task004c_qs_validation.csv
- results/task004c_h5ad_validation.json

## Disk safety
Require at least 50 GB free before execution. Keep both original merged v1 files and the new HGNC v1 files.

## Execute
cd /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Stop rule
Stop after both harmonised merged files validate and the shared-gene audit is complete. Do not execute Task 005 until Task 004c results are reviewed.