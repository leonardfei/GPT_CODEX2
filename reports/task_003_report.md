# Task 003 report — extension-cohort preparation

## Execution boundary

Task 003 was run after synchronizing the local checkout from `origin/main`. The task stopped after audit, subset, corrected QC, and object preparation; no Task 004 annotation or integration was performed.

Large raw inputs and prepared Seurat objects remain on the server under `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/`; only the small audit artifacts are Git-tracked locally.

## Input inventory

The recursive inventory contains 104 uploaded files across CRA002308, nature_xue, and in_house. Primary files used for processing have SHA-256 checksums; uploaded gzip files were tested with `gzip -t`.

## CRA002308

- 14 cell-called 10x Matrix Market samples were validated: N01–N07 and T01–T07, with 36,601 features and recoverable barcodes per sample.
- N01–N07 were harmonized to Adjacent and T01–T07 to Tumor; the numeric suffix supplies the matched patient key.
- The supplementary document states that live nucleated cells were flow-sorted after doublet exclusion from tumor and peri-tumor tissues. This is composition-altering selection, so `abundance_eligible=NO`; the cohort remains eligible for atlas/state analysis.
- Corrected Task 002 identity-independent thresholds were derived from all source cells per sample and applied without using neutrophil marker status to choose thresholds.

## nature_xue

The author object loaded as Seurat 5.3.0 with 20,002 features and 1,092,172 cells. RNA `counts` and `data` layers were present; the author metadata columns included `Sample`, `Cancer_type`, and `clusters`.
The final subset contains 675,539 cells and excludes non-HCC entities, including ICC-linked adjacent samples. It retains HCC tumor cells and AL cells only when the source sample identifier explicitly identifies an HCC sample; AL is harmonized to Adjacent.
Author `clusters` and all original metadata are preserved in namespaced fields, including `source_author_annotation`. The original object was read-only and unchanged. Because uploaded metadata do not document whole-tissue versus enriched/sorted sampling, `abundance_eligible=CONDITIONAL`; no Task 002 QC thresholds were reapplied.

## in_house

- 20 cell-called 10x Matrix Market samples were validated: YJCA01–YJCA10 and YJP01–YJP10, with unique numeric suffix pairing confirmed.
- YJCA was harmonized to Tumor and YJP to Adjacent; clinical fields not present in uploaded metadata remain `unknown`.
- The uploaded files do not establish whether sampling was whole-tissue, enriched, or sorted. Therefore `abundance_eligible=CONDITIONAL` pending provenance, while atlas inclusion is `YES`.

## QC/object validation

The 34 matrix samples produced 402,074 source cells and 392,457 corrected-QC cells. Every validated matrix object has a nonnegative RNA counts layer and globally unique cell IDs.

Prepared server objects:
- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/CRA002308/`
- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/in_house/`
- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/nature_xue/nature_xue_HCC_T_AL_prepared.rds`

## Tracked outputs

- `results/task003_extension_file_inventory.csv`
- `results/task003_extension_sample_manifest.csv`
- `results/task003_extension_cohort_summary.csv`
- `results/task003_extension_qc_audit.csv`
- `results/task003_nature_xue_subset_audit.csv`
- `results/task003_extension_object_validation.csv`

Task 003 is complete at this boundary. Proceed to Task 004 only after reviewing these extension-cohort audits.
