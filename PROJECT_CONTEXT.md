# Project Context

## Project
HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Primary objective
Build a reproducible human HCC Tumor–Adjacent/adjacent-liver scRNA-seq atlas while maximizing technically defensible retention of mature neutrophils.

## Current cohort structure

### Core five cohorts already QC-completed
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340

Corrected Task 002 objects are the authoritative inputs for these five cohorts:
`objects/task002_corrected_seurat/`

### Newly added cohorts
- CRA002308 — uploaded under `raw_data/CRA002308`; expected 7 HCC Tumor and 7 matched adjacent/normal liver samples.
- nature_xue — uploaded under `raw_data/nature_xue`; author-processed Seurat object from Xue et al., Nature 2022. Retain only human HCC Tumor and adjacent liver (AL) cells/samples.
- in_house — uploaded under `raw_data/in_house`; in-house HCC cohort. YJCA denotes Tumor and YJP denotes Adjacent.

### On hold / excluded
- GSE202642
- GSE290298

## Workflow
1. Task 001 — audit/download original public datasets — completed.
2. Task 002 — corrected neutrophil-preserving QC for the original five cohorts — completed and accepted for downstream use.
3. Task 003 — audit, subset and prepare CRA002308, nature_xue and in_house.
4. Task 004 — harmonized broad annotation and neutrophil confirmation across all prepared cohorts.
5. Task 005 — integrate the final approved eight cohorts into a Seurat v5 atlas.
6. Task 006 — optional raw-data reprocessing when required.

## Compute architecture
Control repository: `leonardfei/GPT_CODEX2`

Server root:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/`

The server is not assumed to have GitHub/general internet access.

## Final target object
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_8datasets_integrated_v1.rds`

The final object must preserve provenance fields including dataset, sample, patient, tissue, pairing, etiology, MVI where available, platform/protocol, selection strategy, QC provenance, abundance eligibility, author annotations where applicable, project broad annotation and neutrophil-confidence fields.
