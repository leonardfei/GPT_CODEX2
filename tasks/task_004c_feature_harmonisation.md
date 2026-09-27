# Task 004c — HGNC harmonisation in-place on merged QS and H5AD

## Status
READY TO EXECUTE

## Goal
Standardise gene identifiers directly in the validated 1,490,852-cell merged review objects and overwrite the existing merged v1 QS/H5AD paths after successful validation.

## Files to overwrite
- objects/HCC_TA_8datasets_merged_review_v1.qs
- objects/HCC_TA_8datasets_merged_review_v1.h5ad

No separate *_HGNC_v1 duplicate files should be retained.

## Safety rule
Do not write directly into the live source files while transformation is running.

Instead:
1. create temporary HGNC-standardised QS/H5AD files beside the originals;
2. fully validate both temporary files;
3. only if both validate, atomically replace the original merged v1 files with mv;
4. if either transformation or validation fails, leave the original merged v1 files untouched.

## HGNC mapping
Use the current official HGNC complete set and record retrieval date plus SHA256.

Mapping priority:
1. exact current HGNC-approved symbol;
2. version-stripped Ensembl gene ID;
3. Entrez/NCBI Gene ID;
4. unique previous HGNC symbol;
5. unique HGNC alias.

Ambiguous mappings are not force-resolved. Case-insensitive matching is not used.

If multiple source rows map to the same approved HGNC symbol, sum raw counts exactly.

Unmapped/custom/non-human features are removed from the overwritten merged analysis objects. They remain recoverable from upstream cohort/source objects and the retained Task 004b merge checkpoint.

## QS implementation
Transform the existing task004b_chunked_counts layer chunk-by-chunk with a sparse aggregation matrix.
Preserve all 1,490,852 cells, cell order, and metadata.

Add:
- source_nCount_RNA_before_feature_harmonisation
- source_nFeature_RNA_before_feature_harmonisation
- hgnc_nCount_RNA
- hgnc_nFeature_RNA
- feature_harmonisation = HGNC_approved_symbol

## H5AD implementation
Read the existing 38-GB CSR X matrix in 10,000-cell blocks and write a temporary HGNC-standardised CSR H5AD on disk.

Requirements:
- preserve obs and cell order exactly;
- replace var with unique HGNC-approved symbols and HGNC metadata;
- use int64 CSR indptr;
- cross-validate source and retained count totals against the QS transformation.

## Shared-gene audit
Inspect the eight cohort source feature lists to calculate:
- per-cohort mapped HGNC feature counts;
- strict 8/8 shared HGNC set;
- >=7/8 shared HGNC set.

No separate harmonised cohort objects are created.

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
Require >=50 GB free because temporary validated replacements must coexist with the current files until final atomic replacement.

## Execution
cd /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Stop rule
Stop after the original merged v1 paths have been replaced by validated HGNC-standardised files. Do not execute Task 005 until Task 004c results are reviewed.