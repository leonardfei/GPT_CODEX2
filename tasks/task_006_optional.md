# Task 006 — Optional raw-data reprocessing

## Status
PENDING / CONDITIONAL

## Trigger
Run only if an important cohort cannot be prepared from the uploaded/public processed representation or shows implausibly poor neutrophil recovery and raw sequencing data can be processed reproducibly.

## Requirements
- obtain Web GPT/user approval for the cohort-specific preprocessing plan;
- preserve the existing processed representation as a comparator;
- record software/version/reference/chemistry/cell-calling parameters;
- do not tune solely to maximize neutrophil number;
- never overwrite source files or existing prepared objects.

## Examples
- CRA002308 contains only BAM/FASTQ and requires a count-generation plan;
- a count cohort requires raw-droplet reprocessing to test neutrophil recovery.

## Output
A cohort-specific task report and explicit recommendation on whether the reprocessed representation should replace or supplement the current atlas input.
