# Project Context

## Project
HCC_Peritumoral_Neutrophil_scRNA_Atlas

## Primary objective
Integrate downloadable human hepatocellular carcinoma single-cell transcriptomic datasets containing tumour and adjacent/non-tumour liver tissue, while maximizing technically defensible retention of mature neutrophils.

## Initial phase
1. Audit and download seven public datasets.
2. Perform sample-level neutrophil-friendly preprocessing and QC.
3. Establish broad cell annotations and confirm neutrophils.
4. Integrate eligible datasets into a Seurat v5 object while retaining raw RNA counts and dataset/sample provenance.
5. Reprocess raw sequencing only when needed and feasible.

## Initial accessions
- GSE282701
- GSE242889
- GSE326201
- GSE149614
- GSE299340
- GSE290298
- GSE202642

All cohort/sample counts and pairing information are provisional until Task 001 validates accession-level metadata.

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

The server is not assumed to have GitHub/general internet access. Dataset download feasibility must be tested per source. If a source is blocked, Task 001 must return an explicit local-download/upload plan.

## Final phase-1 object
Target path:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/HCC_TA_7datasets_integrated_v1.rds`

The final object must retain dataset, sample, patient, tissue, pairing, etiology, MVI where available, platform/protocol, QC and neutrophil-confidence metadata.
