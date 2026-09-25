# Task 004 — Harmonized broad annotation and neutrophil confirmation across the expanded atlas

## Status
COMPLETED

## Preconditions
Task 003 extension preparation **and Task 003 metadata correction** reviewed and approved.

Required corrected metadata before this task:
- nature_xue: 92 sample records -> 79 unique patients, including 10 Tumor–Adjacent paired patients;
- nature_xue multiregion samples share base patient IDs (including A074 and A119);
- CRA002308: abundance_eligible = CONDITIONAL, with paired within-cohort abundance sensitivity use only.

## Candidate cohorts
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340
- CRA002308
- nature_xue
- in_house

Only cohorts successfully prepared/approved in Task 003 may enter this task.

## Inputs
For the original five cohorts use corrected Task 002 objects only:
`objects/task002_corrected_seurat/`

For new cohorts use:
`objects/task003_extension/`

## Goal
Create a harmonized broad cell-type layer and confirm neutrophil identities while preserving source annotations and QC provenance.

## Broad classes
At minimum:
- hepatocyte/tumor epithelial
- T/NK
- B
- plasma
- monocyte/macrophage
- neutrophil
- dendritic
- mast
- endothelial
- fibroblast/mesenchymal
- other/uncertain

## Requirements
1. Normalize/analyze each cohort in a way appropriate to its source representation before cross-cohort integration.
2. Use coherent marker programs, not one-gene rules.
3. Distinguish neutrophils from inflammatory monocytes and low-quality myeloid cells.
4. Preserve `source_author_annotation`; never overwrite nature_xue author labels.
5. Create `project_broad_celltype` and `neutrophil_confidence`.
6. Preserve `qc_status`:
   - original five/in_house/CRA count cohorts: distinguish standard pass vs rescued-by-corrected-QC where available;
   - nature_xue: `author_processed`.
7. For every cluster, report the fraction of rescued cells. Flag clusters dominated by rescued cells for marker/QC review.
8. Do not delete a rescued-cell-dominated cluster solely because it was rescued; determine whether it has a coherent biological identity.
9. Report neutrophil counts/fractions by dataset, patient and tissue.
10. Separate abundance roles in all abundance summaries:
   - `YES`: primary abundance analysis if otherwise methodologically eligible;
   - `CONDITIONAL`: report only under its permitted cohort-specific sensitivity rule;
   - `NO`: state/atlas only.
11. Never count multiple nature_xue regions from the same patient as independent patient-level replicates. For A119/A074 and any other multi-region patient, preserve sample-level data but use patient-aware statistics.

## Required outputs
- `results/task004_celltype_counts.csv`
- `results/task004_neutrophil_by_sample.csv`
- `results/task004_cluster_rescued_cell_audit.csv`
- `results/task004_author_vs_project_annotation_crosswalk.csv`
- reviewable UMAP/dotplot/QC figures
- server-side annotated cohort objects
- `reports/task_004_report.md`

## Hold point
Freeze the approved cell set and harmonized broad annotations, then stop. Do not execute Task 005 integration until Web GPT/user review.
