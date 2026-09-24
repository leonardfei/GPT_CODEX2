# Task 002 — Neutrophil-friendly per-sample preprocessing

## Status
PENDING

## Approved datasets
Process only:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Explicitly exclude from this task:
- GSE202642
- GSE290298

Do not preprocess excluded datasets unless Web GPT/user explicitly reactivates them.

## Goal
Create per-sample Seurat objects using defensible QC that maximizes retention of true mature neutrophils without retaining obvious empty droplets/low-quality cells.

## Preconditions
Task 001 is complete and the five-dataset phase-1 cohort has been approved.

## Required principles
- First validate the actual on-server file structure for all five datasets.
- For GSE242889, explicitly verify whether gene/feature and barcode identifiers are recoverable for every sample before ingestion.
- For GSE149614, retain only primary Tumor and Adjacent/NTL cells for the phase-1 cohort; do not include PVTT or lymph-node cells.
- Process each dataset/sample separately before integration.
- Preserve raw counts.
- Calculate nCount_RNA, nFeature_RNA, percent.mt and dataset-appropriate QC.
- Use sample-aware thresholds; do not impose a universal high minimum nFeature threshold.
- Define a conservative pre-QC candidate-neutrophil screen using coherent marker evidence (e.g. FCGR3B, CSF3R, CXCR2, S100A8, S100A9, FPR1 and related context; avoid one-marker rescue).
- Evaluate doublets where technically valid, but document risk of preferentially removing granulocytes.
- Do not finalize QC thresholds until pre-QC candidate-neutrophil distributions have been compared with other cells.

## Mandatory neutrophil-retention audit
For every sample report:
- total cells before QC;
- total cells after QC;
- candidate neutrophils before QC;
- candidate neutrophils after QC;
- neutrophil retention fraction;
- major removal reason(s);
- median/IQR nFeature_RNA for candidate neutrophils and non-neutrophils;
- median/IQR nCount_RNA for candidate neutrophils and non-neutrophils;
- median/IQR percent.mt for candidate neutrophils and non-neutrophils.

Flag any sample with suspiciously low neutrophil retention for review before proceeding.

## Required outputs
- per-sample QC tables and reviewable plots;
- per-sample Seurat RDS objects on server;
- `results/task002_neutrophil_retention_audit.csv`;
- `results/task002_qc_thresholds_by_sample.csv`;
- `reports/task_002_report.md`.

## Hold point
Do not proceed to annotation/integration until Web GPT reviews retention and dataset-specific QC thresholds.
