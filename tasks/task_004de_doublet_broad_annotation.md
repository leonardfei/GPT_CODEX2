# Task 004d/e helper note

## Status
SUPERSEDED AS A SCIENTIFIC SPECIFICATION

The authoritative specification is:

`tasks/task_004d.md`

This file is retained only because the execution wrapper and implementation are split internally into a scDblFinder stage (`task004d_scdblfinder_filter.R`) and a broad-annotation stage (`task004e_broad_annotation.R`).

All scientific parameters, output names, annotation classes, QC requirements and stop rules must follow `tasks/task_004d.md` and Scientific Decision D015.

Execution wrapper:

`bash scripts/bash/task004de_doublet_and_broad_annotation.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas`

Task 005 must remain paused until the canonical Task 004d outputs have been reviewed.
