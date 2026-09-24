# Task 003 metadata correction — nature_xue pairing and CRA002308 abundance eligibility

## Status
COMPLETED

## Scope
This is a metadata-only correction after Task 003 review.

Do **not** redownload data.
Do **not** rerun count-matrix QC.
Do **not** change the retained cell set.
Do **not** execute Task 004.

## Correction 1 — nature_xue patient identity
The current prepared object incorrectly treats every Sample value as a distinct patient.

Use the base A-number as the patient identity:
- A013_HCC + A013_HCC_N -> nature_xue_A013
- A074_HCC + A074_HCC_IM -> nature_xue_A074
- A119_HCC + A119_HCC_IM1 + A119_HCC_IM2 + A119_HCC_N -> nature_xue_A119
- analogous rule for all HCC samples

Strip:
- _HCC
- _HCC_N
- _HCC_IM
- _HCC_IM<number>

Expected after correction:
- 92 sample records
- 79 unique patients
- 10 patients with both Tumor and Adjacent
- 69 Tumor-only patients
- no Adjacent-only patients
- multiregion Tumor patients include A074 and A119

For patients with both Tumor and Adjacent:
- paired_status = paired
- paired_id = patient_id

For Tumor-only patients:
- paired_status = unpaired
- paired_id = unknown

Do not merge multiple Tumor samples from A119 or A074 at the cell/object level. They remain separate sample_id values sharing one patient_id.

## Correction 2 — CRA002308 abundance eligibility
Change CRA002308 from abundance_eligible=NO to abundance_eligible=CONDITIONAL.

Rationale:
The protocol documents flow sorting of live nucleated cells after doublet exclusion. This can alter composition and therefore the cohort should not be pooled as an unbiased whole-tissue absolute abundance dataset. However, no lineage-specific CD45 or immune enrichment is documented in the prepared provenance, and Tumor/Adjacent samples underwent a matched protocol.

Allowed use:
- atlas/state analysis: YES
- paired Tumor–Adjacent abundance sensitivity analysis within CRA002308: YES
- pooled cross-cohort absolute cell-fraction analysis without cohort-level effect modeling: NO

## Implementation
Run:
`scripts/R/task003_correct_metadata.R`

Then update the matrix cohort manifest/QC audit so all CRA002308 rows have:
- abundance_eligible = CONDITIONAL
- the revised selection_strategy text

Regenerate:
- results/task003_extension_sample_manifest.csv
- results/task003_extension_cohort_summary.csv
- reports/task_003_report.md

## Required new output
- results/task003_nature_xue_pairing_audit.csv

## Validation
Verify:
1. nature_xue prepared object still contains exactly 675,539 cells.
2. Cell IDs and RNA counts/data layers are unchanged.
3. nature_xue manifest has 92 sample rows, 79 unique patients and 10 paired patients.
4. A119 HCC/IM1/IM2/N all share nature_xue_A119.
5. CRA002308 still contains the same 158,741 corrected-QC cells across 14 samples.
6. No count/QC thresholds or pass/fail cell lists change.
7. The final Task 003 report explicitly states the corrected abundance rule.

## Hold point
Stop after this correction and push the updated report/audits. Task 004 remains paused until Web GPT/user review.
