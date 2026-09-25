# Task 004 report - harmonized broad annotation and neutrophil confirmation

## Execution status

Task 004 completed on the corrected Task 002 objects and corrected Task 003 extension objects. The source objects were not overwritten. Annotated copies were written server-side under `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task004_annotated/`. Task 005 integration was not executed.

Processed 103 source objects and 1,490,852 annotated cells. The project broad layer classified 65,546 cells as neutrophil; 445,480 cells were marked `rescued_by_corrected_qc`.

## Method and provenance

Each object was read and annotated independently using its RNA counts layer when available, otherwise its RNA data layer. Coherent multi-marker programs were scored for the ten requested lineage classes. Neutrophil confidence required coherent neutrophil marker hits and a score margin against the monocyte/macrophage program; ambiguous myeloid cells were not silently promoted to neutrophils.

`source_author_annotation` was preserved. For `nature_xue`, `qc_status` is `author_processed`; source author annotations remain in the annotated objects and are summarized in the crosswalk output. For corrected count cohorts, cells rescued by corrected QC are explicitly labeled and retained.

Objects without a pre-existing `seurat_clusters` field were audited at the sample review-unit level; `nature_xue` objects with source author labels were audited at the author-annotation-unit level. This is explicitly recorded in `task004_cluster_type` and does not imply a new clustering analysis.

## Dataset summary

- `CRA002308`: 158,741 cells;  9,804 neutrophil-classified; abundance role `CONDITIONAL`.
- `GSE149614`:  63,101 cells;  1,496 neutrophil-classified; abundance role `YES`.
- `GSE242889`:  52,954 cells;  2,811 neutrophil-classified; abundance role `YES`.
- `GSE282701`: 137,518 cells;  5,556 neutrophil-classified; abundance role `YES`.
- `GSE299340`:  76,319 cells;  4,303 neutrophil-classified; abundance role `YES`.
- `GSE326201`:  92,964 cells;  1,573 neutrophil-classified; abundance role `YES`.
- `in_house`: 233,716 cells; 23,768 neutrophil-classified; abundance role `CONDITIONAL`.
- `nature_xue`: 675,539 cells; 16,235 neutrophil-classified; abundance role `CONDITIONAL`.

## Rescued-cell review

There are 34 rescued-dominated review units (rescued fraction > 0.5). None were deleted solely because of rescue status. The cluster audit CSV records the review flag, rescued fraction, neutrophil count, and top broad type.

## Abundance rules

The five original cohorts retain abundance role `YES` for primary analyses subject to their study design. `CRA002308`, `nature_xue`, and `in_house` retain `CONDITIONAL` roles and must not be pooled as unbiased whole-tissue absolute fractions. Patient-level summaries use `patient_id`; multiple `nature_xue` regions remain sample-level records but are not treated as independent patient replicates.

## Required outputs

- `results/task004_celltype_counts.csv`
- `results/task004_neutrophil_by_sample.csv` (sample and patient-tissue summaries)
- `results/task004_cluster_rescued_cell_audit.csv`
- `results/task004_author_vs_project_annotation_crosswalk.csv`
- `figures/task004_marker_program_umap.pdf`, `figures/task004_marker_program_dotplot.pdf`, and `figures/task004_qc_annotation_review.pdf`
- server-side annotated objects under `objects/task004_annotated/`

## Hold point

The corrected cell set and harmonized broad annotation layer are frozen at this Task 004 hold point. Do not execute Task 005 integration until Web GPT/user review.
