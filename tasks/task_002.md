# Task 002 — Corrected neutrophil-preserving per-sample QC

## Status
PENDING

## Purpose
Repeat Task 002 with a corrected QC strategy that is independent of neutrophil-candidate status.

The initial Task 002 run is preserved for comparison in:
`reports/task_002_initial_report.md`

Do not delete or overwrite the initial server-side objects. Write corrected objects to a separate directory.

## Approved datasets
Process only:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Explicitly exclude:
- GSE202642
- GSE290298

## Scientific rationale
The initial Task 002 selected sample-specific lower QC thresholds from the 5th percentile of the pre-QC neutrophil-candidate population and selected mitochondrial ceilings from its 99th percentile. This creates circular validation because candidate cells partly determine the thresholds used to measure their own retention.

The correction must:
1. decouple QC thresholds from neutrophil identity;
2. use permissive secondary QC because most source matrices are already cell-called/filtered;
3. preserve low-RNA granulocytes;
4. prevent extremely permissive mitochondrial thresholds;
5. compare the corrected result directly with the initial Task 002 result.

## Input
Reuse the already staged matrices and metadata from Task 002. Do not redownload data.

## Corrected QC thresholds
For each sample calculate QC metrics from all source cells before applying any neutrophil marker rule.

Use the following independent secondary-QC thresholds:

```text
min_nFeature_RNA =
    max(100, min(300, floor(global sample nFeature_RNA 1st percentile)))

min_nCount_RNA =
    max(200, min(500, floor(global sample nCount_RNA 1st percentile)))

max_percent_mt =
    min(30, max(20, ceil(global sample percent.mt 98th percentile)))
```

These thresholds must be derived from all cells, not from any candidate-neutrophil subset.

Do not raise the lower nFeature threshold above 300 or the lower nCount threshold above 500 without explicit Web GPT approval.

Do not allow mitochondrial thresholds above 30% in this correction run.

If a sample has an unusual distribution that makes these rules clearly invalid, do not silently substitute a new rule. Flag the sample and report it.

## Neutrophil audit definitions

### High-confidence neutrophil candidate
Use only for audit/retention evaluation, never for threshold selection.

Core markers:
- FCGR3B
- CSF3R
- CXCR2
- FPR1

Support markers may include:
- S100A8
- S100A9
- FCAR
- FFAR2
- MCEMP1
- NAMPT
- CTSG
- ELANE
- MPO
- CXCR4

High-confidence rule:
```text
(core_hits >= 2)
OR
(core_hits >= 1 AND support_hits >= 2)
```

### Broad granulocyte-like candidate
Use only as a sensitivity-analysis category.

```text
high_confidence
OR
(support_hits >= 4 AND S100A8 > 0 AND S100A9 > 0)
```

Do not call the broad category "neutrophil" in final biological conclusions before Task 003 annotation.

## Doublet handling
Do not hard-filter cells solely using a high-count/high-feature proxy in this correction run. Keep such cells as review flags. Formal doublet decisions can be revisited after broad annotation if needed.

## Required comparison with initial Task 002
For every sample report:
- source cells;
- cells passing initial Task 002 QC;
- cells passing corrected QC;
- cells rescued by corrected QC;
- cells retained only by initial QC;
- high-confidence neutrophils before QC;
- high-confidence neutrophils after initial QC;
- high-confidence neutrophils after corrected QC;
- broad granulocyte-like cells before and after both QC schemes;
- corrected high-confidence neutrophil retention;
- corrected broad-granulocyte retention.

For rescued cells report:
- median/IQR nFeature_RNA;
- median/IQR nCount_RNA;
- median/IQR percent.mt;
- fraction with >=1 core granulocyte marker;
- fraction satisfying high-confidence neutrophil rule;
- fraction satisfying broad granulocyte-like rule;
- fraction with only S100A8/S100A9 inflammatory evidence and no core marker.

## Mandatory sample flags
Flag a sample for Web GPT review if any of the following occur:
- corrected high-confidence neutrophil retention < 90%;
- corrected QC retains >20 percentage points more total cells than initial QC;
- rescued-cell median percent.mt >25%;
- >50% of rescued granulocyte-like cells have no core granulocyte marker;
- corrected total retained cells exceed source count or any identity mismatch occurs.

## Required outputs
Git-tracked:
- `results/task002_corrected_qc_thresholds_by_sample.csv`
- `results/task002_corrected_neutrophil_retention_audit.csv`
- `results/task002_initial_vs_corrected_qc.csv`
- `results/task002_rescued_cells_summary.csv`
- `results/task002_corrected_rds_validation.csv`
- `figures/task002_corrected_qc_review.pdf`
- `reports/task_002_report.md`

Server-side:
- corrected objects under:
  `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task002_corrected_seurat/<dataset>/<sample_id>_qc.rds`

Do not overwrite:
`objects/task002_seurat/`

## Validation
All corrected RDS objects must:
- preserve raw non-negative RNA counts;
- match the corrected cell list exactly;
- preserve sample/patient/tissue/dataset metadata;
- load successfully in Seurat.

## Hold point
After producing the corrected comparison, stop.

Do not execute Task 003 until Web GPT/user explicitly approves the corrected QC.
