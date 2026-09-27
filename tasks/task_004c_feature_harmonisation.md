# Task 004c — HGNC feature harmonisation before integration

## Status
READY TO EXECUTE

## Goal
Standardise human gene identifiers across the eight HCC scRNA-seq cohorts before Task 005 integration, while preserving all source/raw objects unchanged.

## Reference
Use the current HGNC complete dataset from the official HGNC public download bucket:
https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

Cache it under references/task004c/ and record retrieval date plus SHA256.

## Mapping rules
Apply in this order:
1. exact current HGNC-approved symbol -> keep;
2. Ensembl gene ID -> remove version suffix and map uniquely to HGNC-approved symbol;
3. numeric Entrez/NCBI Gene ID -> map uniquely;
4. previous HGNC symbol -> map only when unique;
5. HGNC alias symbol -> map only when unique;
6. ambiguous previous/alias/ID mapping -> do not force-map;
7. unmapped/custom/non-human features -> preserve in source objects only and exclude from the cross-cohort human integration feature universe.

Do not use case-insensitive matching.

## Duplicate mapping
If multiple source rows map to the same approved HGNC symbol, sum raw counts exactly and retain one row under the approved symbol. Never silently keep only the first duplicate.

## Harmonised cohort objects
Create independent counts+metadata Seurat objects containing HGNC-mapped human genes only:
objects/task004c_feature_harmonized/<dataset>_HGNC_harmonized.qs

Requirements:
- preserve all cells and cell metadata;
- preserve original sample/patient/tissue/QC/annotation provenance;
- preserve pre-harmonisation nCount/nFeature in dedicated metadata fields;
- exactly one RNA counts layer;
- no duplicated gene symbols;
- every row name must be an HGNC-approved symbol;
- no normalization, HVG selection, PCA, UMAP or integration.

Source objects must not be modified.

## Shared feature sets
Calculate the number of HGNC genes per cohort, the strict 8/8 intersection, the >=7/8 set, and pairwise overlaps.
The strict 8/8 HGNC intersection is the candidate universe for Task 005 HVG/integration-feature selection, not the final PCA feature set.
Stop if the strict 8/8 set is <10,000 genes.

## Required outputs
- results/task004c_feature_mapping.csv
- results/task004c_mapping_status_counts.csv
- results/task004c_feature_audit_by_cohort.csv
- results/task004c_hgnc_feature_presence.csv
- results/task004c_pairwise_feature_overlap.csv
- results/task004c_shared_hgnc_features_8of8.txt
- results/task004c_shared_hgnc_features_7plus.txt
- results/task004c_harmonized_objects.csv
- reports/task_004c_report.md

## Execution
cd /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas
bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Stop rule
Stop after harmonised objects and audit outputs validate. Do not execute Task 005 until Task 004c has been reviewed.