# Scientific Decisions

This file stores stable project decisions. Do not alter an existing decision without explicit instruction.

## D001 — Primary biological target
Mature neutrophils are a primary cell type of interest. QC must explicitly measure neutrophil retention rather than optimize only global atlas quality.

## D002 — Statistical unit for abundance
For Tumor–Adjacent abundance comparisons, the biological replicate is the patient/sample, not the individual cell. Prefer paired comparisons within each dataset. Cross-dataset evidence should be synthesized at the dataset/effect-size level rather than by pooling all cells.

## D003 — Composition-biased datasets
Datasets with immune enrichment, FACS gating, CD45 enrichment, artificial remixing, or another composition-altering sampling scheme are not automatically eligible for unbiased whole-tissue abundance estimates. They may still contribute to cell-state and atlas analyses.

## D004 — Raw counts
Raw RNA counts must be preserved whenever available. Integrated expression must never replace source counts. If an author-processed object lacks a raw counts layer, do not fabricate one; record `counts_available=FALSE` and preserve the author-provided expression layer.

## D005 — Integration
The default planned integration method is Seurat v5 RPCA, subject to Task 005 QC and resource feasibility. Biological tissue, etiology, and neutrophil state must not be regressed out merely to improve visual mixing.

## D006 — Multiple testing
Use Benjamini–Hochberg correction unless a task specifies otherwise.

## D007 — Reproducibility
Use fixed seeds for stochastic steps where possible and record them.

## D008 — Cohort inclusion
Current atlas candidates are:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340
- CRA002308
- nature_xue
- in_house

GSE202642 and GSE290298 remain on hold and excluded unless explicitly reactivated.

## D009 — QC threshold independence from neutrophil identity
Neutrophil/granulocyte marker status must not determine QC thresholds used to evaluate neutrophil retention.

For count-matrix cohorts requiring project secondary QC:
- sample-level minimum feature thresholds are bounded to 100–300 features;
- sample-level minimum count thresholds are bounded to 200–500 UMIs/counts;
- sample-level mitochondrial ceilings are bounded to 20–30%;
- thresholds are derived from all cells in the sample.

Neutrophil marker rules are audit tools until broad annotation.

## D010 — New-cohort preparation
CRA002308 and in_house must be audited before processing. If they contain cell-called count matrices, apply the corrected Task 002 QC logic. If only raw FASTQ/BAM or incompatible inputs are present, do not improvise a replacement workflow; report the exact blocker and required preprocessing.

nature_xue is an author-processed Seurat object. Do not discard or overwrite author annotations. Subset to human HCC Tumor and adjacent liver (AL) only, preserve author cell types/cluster labels as separate metadata, and harmonize project metadata before integration.

## D011 — Abundance eligibility is separate from atlas inclusion
A cohort can be included in the integrated atlas while being excluded from whole-tissue abundance meta-analysis. Task 003 must assign and document `abundance_eligible = YES / NO / CONDITIONAL` for each new cohort based on sampling/selection design and available provenance.

## D012 — nature_xue patient identity and CRA002308 abundance use
For nature_xue, patient identity is the base A-number. Sample suffixes `_HCC`, `_HCC_N`, and `_HCC_IM<number>` denote samples/regions from the same base patient rather than separate patients. A patient is paired only when both Tumor and Adjacent samples are present.

For CRA002308, live nucleated-cell flow sorting after doublet exclusion is treated as composition-altering but not as lineage-specific immune enrichment. Therefore:
- atlas/state analysis: allowed;
- paired within-cohort Tumor–Adjacent abundance sensitivity analysis: allowed;
- pooled cross-cohort absolute cell-fraction analysis as if unbiased whole tissue: not allowed.

CRA002308 `abundance_eligible` is therefore `CONDITIONAL`, not `NO`.

## Decision log
- 2026-09-23: D001–D007 initialized.
- 2026-09-24: original five-cohort phase-1 set defined; GSE202642/GSE290298 held.
- 2026-09-24: D009 added after corrected Task 002 review.
- 2026-09-24: D008, D010 and D011 updated/added to expand the atlas with CRA002308, nature_xue and in_house.
- 2026-09-24: D012 added after Task 003 review to correct nature_xue patient/pairing metadata and set CRA002308 abundance eligibility to CONDITIONAL.
