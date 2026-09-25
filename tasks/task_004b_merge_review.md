# Task 004b — Merge eight cohorts into one Seurat review object

## Status
READY TO RUN

## Goal
Create one Seurat object containing all eight currently prepared cohorts so the user can manually inspect the existing Task 004 annotations before final re-annotation or integration.

## Scope
This task is merge only. Do not run RPCA, Harmony, CCA, normalization, PCA, UMAP, reclustering, annotation changes, cell removal, or downsampling.

Use the 103 Task 004 annotated Seurat objects as inputs. Retain all 1,490,852 cells, RNA counts, and metadata. Source-specific PCA/UMAP/graphs must be discarded because they are not directly comparable across independently processed objects.

The current project_broad_celltype and neutrophil_confidence fields must be retained strictly as preliminary review labels. Add annotation_status=preliminary_unvalidated_task004.

Preserve source_author_annotation from nature_xue.

Correct the known Task 004 rescued-status metadata bug in the merged review object only:
- CRA002308 and in_house: rescued_by_corrected_qc=FALSE; qc_status=corrected_qc_pass.
- nature_xue: rescued_by_corrected_qc=FALSE; qc_status=author_processed.
- Preserve the original Task 004 values in task004_qc_status_premerge and task004_rescued_premerge.

## Primary output
/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.rds

This is not the final integrated atlas.

## Required outputs
- results/task004b_merge_review_cohort_summary.csv
- results/task004b_merge_review_dataset_summary.csv
- results/task004b_merge_review_metadata_fields.csv
- results/task004b_merge_review_validation.csv
- reports/task_004b_merge_review_report.md

## Validation
- exactly 1,490,852 cells
- exactly 8 datasets
- globally unique cell IDs
- RNA counts retained
- current preliminary annotations retained
- nature_xue author annotations retained
- no source object overwritten
- no integration or dimensionality reduction performed

## Execution
Rscript scripts/R/task004b_merge_review.R --project-root /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

Stop after validation. Do not execute Task 005.