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

## Decision log
- 2026-09-23: D001–D007 initialized for project start.
- 2026-09-24: D008 added by user decision; GSE202642 and GSE290298 excluded from phase-1 analysis/integration.
