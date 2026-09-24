# Task 003 — Broad annotation and neutrophil confirmation

## Status
PENDING

## Approved datasets
Only:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Do not include GSE202642 or GSE290298.

## Goal
Establish robust broad cell classes within each approved dataset and confirm neutrophils before cross-dataset integration.

## Broad classes
At minimum: hepatocyte/tumour epithelial, T/NK, B, plasma, monocyte/macrophage, neutrophil, endothelial, fibroblast/mesenchymal, and other/uncertain where needed.

## Requirements
- Annotate using coherent marker programs, not one-gene rules.
- Distinguish neutrophils from inflammatory monocytes and low-quality myeloid cells.
- Record annotation confidence and uncertainty.
- Report neutrophil counts and fractions by patient, tissue and dataset.
- Preserve the patient-paired structure for downstream abundance testing.
- Freeze the approved cell set and broad annotations entering Task 004.

## Outputs
- `results/task003_celltype_counts.csv`
- `results/task003_neutrophil_by_sample.csv`
- server-side annotated per-dataset objects
- `reports/task_003_report.md`

Stop for review before Task 004.
