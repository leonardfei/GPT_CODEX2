# Project Status

## Current task
Task 004b — eight-cohort unintegrated Seurat merge for manual annotation review — READY TO RUN

## Purpose
At user request, create one unintegrated Seurat object containing all eight cohorts before correcting the Task 004 annotation or running Task 005 integration.

Target object:
/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_merged_review_v1.rds

Expected:
- 8 cohorts
- 103 Task 004 annotated source objects
- 1,490,852 cells

## Task 004 annotation status
Task 004 computational execution completed, but its broad annotation is preliminary/unvalidated. The merge object must preserve these labels for user inspection without treating them as final.

Known issues retained for review:
1. Task 004 rescued status was wrong for CRA002308/in_house; Task 004b corrects only this metadata in the merged review object.
2. Task 004 marker-program annotation disagrees systematically with the Xue author reference for several lineages, especially neutrophil, dendritic and monocyte/macrophage.

## Pending
1. Execute Task 004b
2. User manually reviews annotation
3. Rebuild/correct annotation as needed
4. Run Task 005 integration only after annotation review

## Next execution command
Execute task_004b.

## Last update
2026-09-25