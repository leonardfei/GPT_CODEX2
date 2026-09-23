# Task 002 — Neutrophil-friendly per-sample preprocessing

## Status
PENDING

## Goal
Create per-sample Seurat objects using defensible QC that maximizes retention of true mature neutrophils without retaining obvious empty droplets/low-quality cells.

## Preconditions
Task 001 reviewed and approved by Web GPT.

## Required principles
- Process each dataset/sample separately before integration.
- Preserve raw counts.
- Calculate nCount_RNA, nFeature_RNA, percent.mt and dataset-appropriate QC.
- Use sample-aware thresholds; do not impose a universal high minimum nFeature threshold.
- Define a conservative pre-QC candidate-neutrophil screen using coherent marker evidence (e.g. FCGR3B, CSF3R, CXCR2, S100A8, S100A9, FPR1, MNDA as context; avoid one-marker rescue).
- Evaluate doublets where technically valid, but document risk of removing granulocytes.
- Generate a per-sample neutrophil-retention audit before accepting filters.

## Required outputs
- per-sample QC tables and plots;
- per-sample Seurat RDS objects on server;
- results/task002_neutrophil_retention_audit.csv;
- reports/task_002_report.md.

## Hold point
Do not proceed to annotation/integration until Web GPT reviews retention and dataset-specific QC thresholds.
