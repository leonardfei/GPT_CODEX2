# Project Status

## Current task
Task 004d — per-sample scDblFinder filtering and corrected broad cell-type annotation — BLOCKED / APPROVED

The canonical implementation is ready, and the read-only input preflight passed on 2026-09-28. Execution is currently blocked before Phase A: the compute-server dependency installation had not completed at the last verified check, and server access must be re-established through a safe authentication method. No doublet calls, singlet object, or corrected broad annotation have been generated. See `reports/task_004d_report.md`.

Task 004c is COMPLETED:
- 1,490,852 cells;
- 38,025 unique HGNC-approved symbols;
- 8 datasets, 194 samples, 132 patients;
- 17,680 HGNC genes shared across all 8 cohorts;
- Task 005 has not been executed.

## Task 004d approved input
Use exactly:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_merged_review_v1.qs`

If this exact path is absent, stop and report available candidate paths rather than silently substituting another object.

## Task 004d design
Phase A:
- run scDblFinder independently for each `project_sample_id`;
- use raw counts;
- explicitly set `dbr.per1k=0.008`;
- use cluster-based mode for samples >=500 cells and random mode for smaller samples;
- record all scores/classes and neutrophil-retention QC;
- remove only cells called `scDblFinder.class == "doublet"`.

Phase B:
- build a new singlet-only merged QS object without overwriting the Task 004c input.

Phase C:
- replace the preliminary Task 004 per-cell score annotation with a cluster-based broad annotation;
- preprocess/cluster within each dataset without batch integration;
- use cluster markers + canonical lineage programs + Xue author labels as a reference anchor;
- retain ambiguous clusters as `Uncertain/Mixed`;
- do not call malignancy at this stage.

Approved task specification:
`tasks/task_004d.md`

## Implementation readiness
The canonical Task 004d implementation is now present in the control repository:
- `scripts/R/task004d_scdblfinder_filter.R`
- `scripts/R/task004e_broad_annotation.R`
- `scripts/R/task004d_finalize_report.R`
- `scripts/R/task004de_install_deps.R`
- `scripts/bash/task004de_doublet_and_broad_annotation.sh`
- `config/task004d_xue_author_to_broad.tsv`

Static delimiter/quote checks passed for all implementation scripts. Server-side execution has not yet been recorded; Task 004d remains incomplete until the runtime outputs validate and the final report is committed.

## Next Codex command
`Execute task_004d.`

## Hold point
Task 005 remains PAUSED until Task 004d doublet filtering and v2 broad annotation have been reviewed.
