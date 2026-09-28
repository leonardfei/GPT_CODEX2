# Task 005 — Eight-cohort Seurat v5 integration

## Status
PENDING — BLOCKED UNTIL TASK 004d DOUBLET FILTERING AND BROAD ANNOTATION ARE REVIEWED

## Goal
Integrate all cohorts approved after Tasks 003–004 into one traceable HCC Tumor–Adjacent atlas while preserving raw/source expression, biological provenance and abundance-eligibility metadata.

## Intended cohorts
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340
- CRA002308
- nature_xue
- in_house

If Task 003 or Task 004 marks a cohort technically incompatible or blocked, do not silently omit it. Record the blocker and stop for review.

## Task 004d prerequisite
Do not start integration from the unfiltered Task 004c merged object.

Use the validated Task 004d singlet-only object with `broad_celltype_v2` annotation as the biological starting point. Task 004 preliminary labels are audit-only.

Expected Task 004d input for this task:
- `objects/merge/HCC_TA_8datasets_singlets_broad_v1.qs`
- and, when validated, `objects/merge/HCC_TA_8datasets_singlets_broad_v1.h5ad`.

Task 005 must preserve scDblFinder provenance and the v2 broad annotation.

## Feature harmonisation prerequisite
Task 004d inherits the validated HGNC-standardised feature space from Task 004c. Do not restart Task 005 from the older unfiltered Task 004c merged object.

Candidate integration feature universe:
`results/task004c_shared_hgnc_features_8of8.txt`

Select HVGs/integration features from this shared HGNC universe while using the Task 004d singlet-only object as the cell-level input. For raw-union/custom-feature provenance checks, use upstream cohort/source objects or the retained Task 004b checkpoint.

## Default approach
Use Seurat v5 RPCA integration as the default, with parameters chosen after inspecting object sizes, shared features and normalization state.

Do not downsample cells silently.

If the full nature_xue HCC Tumor/AL subset makes full-cell RPCA infeasible under available RAM, stop and report resource requirements plus a scientifically defensible alternative (e.g. reference mapping/sketch strategy). Do not substitute an approximation without approval.

## Required checks
- globally unique cell IDs;
- common metadata schema;
- source counts preserved where available;
- author-processed expression preserved for nature_xue;
- counts availability explicitly recorded;
- dataset/sample/patient/tissue provenance retained;
- `source_author_annotation`, `broad_celltype_v2`, `broad_annotation_confidence`, scDblFinder provenance, `qc_provenance`, `qc_status`, and `abundance_eligible` retained;
- unintegrated representation preserved alongside integrated reduction where practical;
- integration inspected by dataset, patient, tissue, broad cell type and QC provenance;
- overcorrection/undercorrection assessed;
- tissue, etiology and neutrophil state must not be regressed out merely to improve mixing.

## Abundance rule
The integrated object may contain atlas-only cohorts. Whole-tissue Tumor–Adjacent abundance analyses must subset to `abundance_eligible == YES` and use patient/sample as the biological replicate.

## Final object
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_integrated_v1.rds`

## Required outputs
- final integrated Seurat object;
- object checksum;
- object size, total cells/features and counts by dataset/tissue;
- integration QC figures/tables;
- `results/task005_integrated_cell_counts.csv`;
- `results/task005_integration_qc_summary.csv`;
- `reports/task_005_report.md`.

Stop after final object validation.
