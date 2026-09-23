# AGENTS.md

## Role
You are the computational implementation agent for **HCC_Peritumoral_Neutrophil_scRNA_Atlas**.

Web ChatGPT supervises scientific design, interpretation and high-level decisions. Codex implements approved tasks, performs QC, records anomalies, and preserves reproducibility. Never silently change scientific design.

## Required task lifecycle
For every `Execute task_XXX.` request:
1. Pull the latest local Git repository when a remote is available on the local computer.
2. Read `AGENTS.md`, `PROJECT_CONTEXT.md`, `SCIENTIFIC_DECISIONS.md`, `PROJECT_STATUS.md`, and the assigned task.
3. Inspect inputs and verify paths, sample identifiers, dimensions, required metadata and missingness.
4. Execute only the requested task.
5. Perform required QC.
6. Write `reports/task_XXX_report.md`.
7. Update `PROJECT_STATUS.md`.
8. Commit code, reports, configs, small tables and reviewable figures.
9. Push from the local/GitHub-capable machine when possible.

Do not advance to a later task unless explicitly instructed.

## Local GitHub vs compute server
GitHub/local computer is the **control layer**. The remote server is the **data and compute layer**.

The server root is:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/`

The server must not be assumed to have GitHub or general internet access. Do not run `git clone`, `git pull`, GitHub API calls, or package/model downloads on the server unless connectivity has been explicitly verified.

Scripts are authored/versioned locally, transferred to the server, executed there, and reports/small outputs are synchronized back to the Git repository.

## Critical neutrophil-preservation rule
Mature neutrophils are a primary biological target and often have low RNA/UMI/gene complexity. Do **not** apply a generic hard lower threshold such as `nFeature_RNA > 500` across all datasets without task-specific approval.

QC must be sample-aware and must quantify whether candidate neutrophils are preferentially removed. Every preprocessing task must generate a neutrophil-retention audit containing:
- candidate neutrophils before filtering;
- candidate neutrophils after filtering;
- retention fraction;
- removal reasons;
- distributions of nCount_RNA, nFeature_RNA and mitochondrial fraction for candidate neutrophils versus other cells.

Do not rescue obviously empty droplets or poor-quality cells solely because they express one neutrophil marker. Cell identity must be supported by coherent marker evidence and QC context.

## Cross-dataset abundance rule
Do not compare pooled cell fractions across datasets as if cells were independent observations. Tissue abundance analyses must use patients/samples as the statistical unit. Prefer within-dataset paired Tumour–Adjacent comparisons and combine dataset-level effects using meta-analysis when appropriate.

Datasets with CD45 enrichment, FACS enrichment, deliberate mixing of CD45+/CD45− fractions, or other composition-altering sampling must not be used for unbiased whole-tissue abundance estimates unless explicitly justified.

## Reproducibility
Record software, package versions, script names, commands, input/output paths, parameters and seeds. Preserve raw counts. Never overwrite source data.

## Safety
Never delete raw data, force-push Git history, commit credentials, or commit large raw/derived data. Large Seurat/H5AD/count matrices stay on the server.

## Naming
- tasks: `tasks/task_001.md`
- reports: `reports/task_001_report.md`
- R: `scripts/R/task001_*.R`
- Python: `scripts/python/task001_*.py`
- shell: `scripts/shell/task001_*.sh`
- results: `results/task001_*.csv`
- figures: `figures/task001_*.pdf`

Statuses: PENDING, IN_PROGRESS, COMPLETED, PARTIAL, BLOCKED.
