# Task 004 — Five-dataset Seurat integration

## Status
PENDING

## Approved datasets
Integrate only:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Explicitly exclude:
- GSE202642
- GSE290298

## Goal
Integrate approved cells from the five datasets into one traceable Seurat v5 object without erasing raw counts or biologically relevant Tumour–Adjacent variation.

## Default integration
Seurat v5 RPCA, with parameters determined and documented during implementation/QC.

## Required checks
- retain original RNA counts and source sample layers/provenance;
- verify metadata completeness and uniqueness;
- inspect integration by dataset, patient, tissue and broad cell type;
- assess overcorrection and undercorrection;
- do not regress out tissue, etiology or neutrophil state merely to improve visual mixing;
- preserve an unintegrated representation alongside the integrated reduction when practical.

## Final object
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_5datasets_integrated_v1.rds`

## Required metadata
dataset, sample_id, patient_id, tissue, paired_status/paired_id, etiology, MVI where available, platform, selection_strategy, QC metrics, broad cell type and neutrophil-confidence fields.

## Outputs
- final Seurat object on server;
- integration QC figures/tables;
- `reports/task_004_report.md`;
- object checksum, size, cell count, feature count and Seurat/package versions recorded in the Git-tracked report.

Stop after final object validation.
