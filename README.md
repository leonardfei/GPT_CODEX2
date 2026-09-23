# GPT_CODEX2

GPT–Codex research workflow for the project **HCC_Peritumoral_Neutrophil_scRNA_Atlas**.

## Scientific objective
Build a reproducible human HCC tumour–adjacent scRNA-seq atlas that preserves mature neutrophils as far as technically possible, then integrate eligible public datasets into a Seurat object for downstream peritumoral-neutrophil analyses.

## Architecture
- **GitHub/local computer = control layer:** tasks, code, configuration, reports, QC summaries, small results and version history.
- **Compute server = data/compute layer:** downloaded public data, large intermediate objects, Seurat objects and heavy computation.
- The compute server must not require GitHub access.

Server project root:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/`

## Initial datasets
GSE282701, GSE242889, GSE326201, GSE149614, GSE299340, GSE290298, GSE202642.

Start with `tasks/task_001.md`.