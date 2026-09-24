# Project Context

## Project
HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Primary objective
Build a reproducible tumour–adjacent human HCC scRNA-seq atlas while maximizing technically defensible retention of mature neutrophils.

## Phase-1 approved cohort
The current analysis/integration cohort contains five datasets:
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Two audited datasets are retained on hold and are not included in Tasks 002–004:
- GSE290298
- GSE202642

## Workflow
1. Task 001: audit/download seven candidate datasets — completed.
2. Task 002: neutrophil-friendly QC on the five approved datasets.
3. Task 003: broad annotation and neutrophil confirmation.
4. Task 004: integrate the five approved datasets into a Seurat v5 object.
5. Task 005: raw-data reprocessing only when needed.

## Compute architecture
Control repository: `leonardfei/GPT_CODEX2`

Server project root:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/`

Expected server directories:
- raw_data/
- processed_data/
- objects/
- results/
- figures/
- logs/
- tmp/
- scripts/

The server is not assumed to have GitHub/general internet access.

## Final phase-1 object
Target path:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_5datasets_integrated_v1.rds`

The final object must retain dataset, sample, patient, tissue, pairing, etiology, MVI where available, platform/protocol, QC and neutrophil-confidence metadata.
