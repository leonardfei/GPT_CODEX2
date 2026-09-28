# Task 004d report — per-sample doublet filtering and corrected broad annotation

## Status

**BLOCKED — read-only preflight completed; analysis phases not run.**

Task 005 was not started.

## Scope and input preflight

The exact approved input was found at:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_merged_review_v1.qs`

The preflight passed the task's structural and metadata checks:

- 1,490,852 cells and 38,025 unique HGNC feature symbols;
- 8 datasets, 194 unique `project_sample_id` captures, and 132 patients;
- one RNA counts layer stored as `task004b_chunked_counts`, with cell order matching the Seurat object;
- no duplicate cell or feature identifiers, missing required metadata, or zero-count/zero-feature metadata entries;
- 65,546 cells carry the preliminary Task 004 neutrophil label, for audit only.

The source file was 4,562,511,399 bytes. Its SHA-256 was `1ac8c7ba371879509fa44b1cb2626f0e85fae5da12a035dee27879a097eb6352`.

Capture size ranged from 179 to 20,399 cells. The only capture below 500 cells was `nature_xue::A105_HCC` (179 cells); the approved small-sample mode therefore remains `clusters=NULL` for this capture if/when scoring runs. No sample was merged with another capture.

Preflight tables are in `results/task004d_input_preflight.csv`, `results/task004d_input_metadata_missingness.csv`, `results/task004d_preflight_by_sample.csv`, and `results/task004d_preflight_by_dataset_tissue.csv`.

## Execution and blocker

The compute host uses R 4.3.3. At the last verified dependency check, `xgboost` 1.7.11.1 was available in the task-specific library, while `scDblFinder`, `SingleCellExperiment`, and `BiocParallel` were still missing. The installer process was still present but idle/waiting after approximately 18 minutes; the server-side installation log could not be retrieved during the subsequent access interruption. Its final state is therefore unknown.

The local installer is `scripts/R/task004d_install_deps.R`; it targets the R 4.3 / Bioconductor 3.18 compatibility lane, checks the `xgboost` compatibility guard, and smoke-tests both approved scDblFinder modes before analysis. No alternative doublet caller was substituted.

No scDblFinder call was completed. Consequently, no filtering was performed and the following were **not** created: the all-cell call table, doublet QC tables/figure, neutrophil retention audit, singlet QS object, cluster/marker outputs, corrected broad labels, or Phase C annotation figure. The canonical source QS was not modified. Server-side output state after the interrupted installer is unverified.

## Software and reproducibility record

Preflight runtime: R 4.3.3, qs 0.27.3, Seurat 5.3.0, SeuratObject 5.2.0. The preflight script is `scripts/R/task004d_preflight.R`; the dependency installer is `scripts/R/task004d_install_deps.R`. No analysis seed was consumed because Phase A did not start. Approved future scoring parameters remain those in `tasks/task_004d.md` and `SCIENTIFIC_DECISIONS.md` D015.

## Resume requirement

Re-establish authorized compute-host access with a safe authentication method (prefer an SSH key), inspect the existing installer process and log before rerunning anything, finish dependency verification/smoke tests, then resume the task from Phase A. Do not rerun preflight in place of the required analysis, and do not start Task 005.
