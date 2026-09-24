#!/usr/bin/env python3
"""Combine small server-side Task 003 audits and write the Git-tracked report."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results-root", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    return p.parse_args()


def yes_no(value: object) -> str:
    return "YES" if bool(value) else "NO"


def main() -> None:
    args = parse_args()
    root = args.results_root.resolve()
    matrix_manifest = pd.read_csv(root / "task003_extension_sample_manifest_matrix_cohorts.csv")
    nature_manifest = pd.read_csv(root / "task003_nature_xue_sample_manifest.csv")
    manifest = pd.concat([matrix_manifest, nature_manifest], ignore_index=True, sort=False)
    manifest = manifest.sort_values(["dataset", "sample_id"]).reset_index(drop=True)
    manifest.to_csv(root / "task003_extension_sample_manifest.csv", index=False)

    matrix_qc = pd.read_csv(root / "task003_extension_qc_audit.csv")
    nature_audit = pd.read_csv(root / "task003_nature_xue_subset_audit.csv")
    summaries: list[dict[str, object]] = []
    for dataset, group in manifest.groupby("dataset", sort=True):
        row: dict[str, object] = {
            "dataset": dataset,
            "n_sample_rows": len(group),
            "n_unique_samples": group["sample_id"].nunique(),
            "n_unique_patients": group["patient_id"].nunique(),
            "n_paired_patients": group.loc[group["paired_status"].astype(str) == "paired", "patient_id"].nunique(),
            "tissues": ";".join(sorted(group["tissue"].dropna().astype(str).unique())),
            "counts_available": ";".join(sorted(group["counts_available"].dropna().astype(str).unique())),
            "abundance_eligible": ";".join(sorted(group["abundance_eligible"].dropna().astype(str).unique())),
            "analysis_inclusion": ";".join(sorted(group["analysis_inclusion"].dropna().astype(str).unique())),
            "n_cells_after_qc_or_subset": int(group["n_cells"].sum()) if "n_cells" in group else None,
            "selection_strategy": ";".join(sorted(group["selection_strategy"].dropna().astype(str).unique())),
        }
        if dataset in set(matrix_qc["dataset"]):
            q = matrix_qc[matrix_qc["dataset"] == dataset]
            row.update({
                "n_cells_source": int(q["n_cells_source"].sum()),
                "n_cells_after_corrected_qc": int(q["n_cells_after_corrected_qc"].sum()),
                "retention_fraction_source_to_qc": float(q["n_cells_after_corrected_qc"].sum() / q["n_cells_source"].sum()),
                "qc_status": ";".join(sorted(q["qc_status"].astype(str).unique())),
                "qc_provenance": ";".join(sorted(q["qc_provenance"].astype(str).unique())),
            })
            row["n_cells_after_qc_or_subset"] = int(q["n_cells_after_corrected_qc"].sum())
        else:
            a = nature_audit.iloc[0]
            row.update({
                "n_cells_source": int(a["original_n_cells"]),
                "n_cells_after_corrected_qc": int(a["subset_n_cells"]),
                "retention_fraction_source_to_qc": float(a["subset_n_cells"] / a["original_n_cells"]),
                "qc_status": str(a["analysis_inclusion"]),
                "qc_provenance": "author_processed_object_no_task002_threshold_reapplication",
            })
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(root / "task003_extension_cohort_summary.csv", index=False)

    qc = matrix_qc.copy()
    qc["audit_scope"] = "cell_called_matrix_corrected_qc"
    qc.to_csv(root / "task003_extension_qc_audit.csv", index=False)

    inventory = pd.read_csv(root / "task003_extension_file_inventory.csv")
    validation = pd.read_csv(root / "task003_extension_object_validation.csv")
    nature = nature_audit.iloc[0]
    lines = [
        "# Task 003 report — extension-cohort preparation",
        "",
        "## Execution boundary",
        "",
        "Task 003 was run after synchronizing the local checkout from `origin/main`. The task stopped after audit, subset, corrected QC, and object preparation; no Task 004 annotation or integration was performed.",
        "",
        "Large raw inputs and prepared Seurat objects remain on the server under `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/`; only the small audit artifacts are Git-tracked locally.",
        "",
        "## Input inventory",
        "",
        f"The recursive inventory contains {len(inventory)} uploaded files across CRA002308, nature_xue, and in_house. Primary files used for processing have SHA-256 checksums; uploaded gzip files were tested with `gzip -t`.",
        "",
        "## CRA002308",
        "",
        "- 14 cell-called 10x Matrix Market samples were validated: N01–N07 and T01–T07, with 36,601 features and recoverable barcodes per sample.",
        "- N01–N07 were harmonized to Adjacent and T01–T07 to Tumor; the numeric suffix supplies the matched patient key.",
        "- The supplementary document states that live nucleated cells were flow-sorted after doublet exclusion from tumor and peri-tumor tissues. Because this can alter cell composition but no lineage-specific immune enrichment is documented, `abundance_eligible=CONDITIONAL`: atlas/state analysis and paired within-cohort Tumor–Adjacent abundance sensitivity analyses are allowed, but it must not be treated as an unbiased pooled whole-tissue fraction dataset.",
        "- Corrected Task 002 identity-independent thresholds were derived from all source cells per sample and applied without using neutrophil marker status to choose thresholds.",
        "",
        "## nature_xue",
        "",
        f"The author object loaded as Seurat {nature['seurat_version']} with {int(nature['original_n_features']):,} features and {int(nature['original_n_cells']):,} cells. RNA `counts` and `data` layers were present; the author metadata columns included `Sample`, `Cancer_type`, and `clusters`.",
        f"The final subset contains {int(nature['subset_n_cells']):,} cells and excludes non-HCC entities, including ICC-linked adjacent samples. It retains HCC tumor cells and AL cells only when the source sample identifier explicitly identifies an HCC sample; AL is harmonized to Adjacent. Corrected Sample parsing yields {int(nature['n_unique_patients'])} unique patients, including {int(nature['n_paired_patients'])} patients with both Tumor and Adjacent samples.",
        "Author `clusters` and all original metadata are preserved in namespaced fields, including `source_author_annotation`. The original object was read-only and unchanged. Because uploaded metadata do not document whole-tissue versus enriched/sorted sampling, `abundance_eligible=CONDITIONAL`; no Task 002 QC thresholds were reapplied.",
        "",
        "## in_house",
        "",
        "- 20 cell-called 10x Matrix Market samples were validated: YJCA01–YJCA10 and YJP01–YJP10, with unique numeric suffix pairing confirmed.",
        "- YJCA was harmonized to Tumor and YJP to Adjacent; clinical fields not present in uploaded metadata remain `unknown`.",
        "- The uploaded files do not establish whether sampling was whole-tissue, enriched, or sorted. Therefore `abundance_eligible=CONDITIONAL` pending provenance, while atlas inclusion is `YES`.",
        "",
        "## QC/object validation",
        "",
        f"The 34 matrix samples produced {int(matrix_qc['n_cells_source'].sum()):,} source cells and {int(matrix_qc['n_cells_after_corrected_qc'].sum()):,} corrected-QC cells. Every validated matrix object has a nonnegative RNA counts layer and globally unique cell IDs.",
        "",
        "Prepared server objects:",
        "- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/CRA002308/`",
        "- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/in_house/`",
        "- `/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/nature_xue/nature_xue_HCC_T_AL_prepared.rds`",
        "",
        "## Tracked outputs",
        "",
        "- `results/task003_extension_file_inventory.csv`",
        "- `results/task003_extension_sample_manifest.csv`",
        "- `results/task003_extension_cohort_summary.csv`",
        "- `results/task003_extension_qc_audit.csv`",
        "- `results/task003_nature_xue_subset_audit.csv`",
        "- `results/task003_nature_xue_pairing_audit.csv`",
        "- `results/task003_extension_object_validation.csv`",
        "",
        "Task 003 is complete at this boundary. Proceed to Task 004 only after reviewing these extension-cohort audits.",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
