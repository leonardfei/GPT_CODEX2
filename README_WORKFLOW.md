# GPT–Codex Workflow

## Separation of responsibilities
**Web GPT:** scientific design, task specification, review, interpretation and next-step decisions.

**Codex:** code implementation, server execution, QC, numerical outputs, reports and reproducibility.

**GitHub/local computer:** shared control layer and version history.

**Compute server:** data and heavy computation only.

## Standard loop
`Web GPT → task_XXX.md → GitHub → Codex → server execution → QC/report → GitHub → Web GPT review`

Because the compute server may not reach GitHub, Codex must not assume the server can pull this repository. Transfer only the required scripts/configs to the server, run them there, and sync back small outputs/reports.

## Web GPT primarily reads
- PROJECT_CONTEXT.md
- SCIENTIFIC_DECISIONS.md
- PROJECT_STATUS.md
- reports/task_XXX_report.md

## Codex execution trigger
When instructed `Execute task_XXX.`, Codex must execute the full lifecycle for that task and stop before the next task unless instructed otherwise.
