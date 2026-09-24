# Scientific Decisions

This file stores stable project decisions. Do not alter an existing decision without explicit instruction.

## D001 — Primary biological target
Mature neutrophils are a primary cell type of interest. QC must explicitly measure neutrophil retention rather than optimize only global atlas quality.

## D002 — Statistical unit for abundance
For Tumour–Adjacent abundance comparisons, the biological replicate is the patient/sample, not the individual cell. Prefer paired comparisons within each dataset. Cross-dataset evidence should be synthesized at the dataset/effect-size level rather than by pooling all cells.

## D003 — Composition-biased datasets
Datasets with immune enrichment, FACS gating, CD45 enrichment, or artificial remixing of cell fractions may be used for cell-state analyses but are not automatically eligible for unbiased whole-tissue cell-abundance comparisons.

## D004 — Raw counts
Raw RNA counts must be preserved in the integrated Seurat object. Integrated values must never replace source counts.

## D005 — Integration
The default planned integration method is Seurat v5 RPCA, subject to Task 004 QC. Biological tissue and etiology signals must not be intentionally regressed out without a separate scientific decision.

## D006 — Multiple testing
Use Benjamini–Hochberg correction unless a task specifies otherwise.

## D007 — Reproducibility
Use fixed seeds for stochastic steps where possible and record them.

## D008 — Phase-1 analysis cohort
For Tasks 002–004, include only the following five datasets:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Do not include GSE202642 or GSE290298 in preprocessing, annotation, abundance analysis, or Seurat integration at this stage. Retain their downloaded/metadata files unchanged for possible later supplementary analyses.

## D009 — QC threshold independence from neutrophil identity
Neutrophil or granulocyte marker status must not determine the QC thresholds used to evaluate neutrophil retention.

For the corrected Task 002 secondary QC:
- sample-level minimum feature thresholds are bounded to 100–300 features;
- sample-level minimum count thresholds are bounded to 200–500 UMIs/counts;
- sample-level mitochondrial ceilings are bounded to 20–30%;
- thresholds are derived from all cells in the sample, not from candidate-neutrophil subsets.

Neutrophil marker rules are audit tools only until Task 003 annotation.

The initial Task 002 QC run is preserved for comparison but is superseded for downstream analysis until the corrected run is approved.

## Decision log
- 2026-09-23: D001–D007 initialized for project start.
- 2026-09-24: D008 added by user decision; GSE202642 and GSE290298 excluded from phase-1 analysis/integration.
- 2026-09-24: D009 added after Web GPT review of initial Task 002; QC threshold selection must be independent of neutrophil-candidate status.
