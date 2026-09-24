#!/usr/bin/env python3
"""Run the corrected, identity-independent Task 002 secondary QC.

The initial Task 002 objects and cell-level QC files are treated as historical
comparison artifacts. This run re-reads the existing server-side count inputs,
derives thresholds from all source cells only, and writes corrected staged
matrices under a separate output root.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from scipy.sparse import csr_matrix, vstack


SCRIPT_DIR = Path(__file__).resolve().parent
INITIAL_SCRIPT = SCRIPT_DIR / "task002_preprocess_qc.py"
spec = importlib.util.spec_from_file_location("task002_initial", INITIAL_SCRIPT)
if spec is None or spec.loader is None:
    raise ImportError(f"Cannot import initial Task 002 implementation: {INITIAL_SCRIPT}")
initial = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = initial
spec.loader.exec_module(initial)

APPROVED_DATASETS = initial.APPROVED_DATASETS
CORE_MARKERS = initial.CORE_MARKERS
SUPPORT_MARKERS = initial.SUPPORT_MARKERS
MatrixData = initial.MatrixData
SampleSpec = initial.SampleSpec


def parse_args() -> object:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    return parser.parse_args()


def corrected_thresholds(
    n_features: np.ndarray, n_counts: np.ndarray, percent_mt: np.ndarray
) -> dict[str, object]:
    """Derive thresholds from all source cells, with explicit bounds."""

    feature_q01 = float(np.quantile(n_features, 0.01))
    count_q01 = float(np.quantile(n_counts, 0.01))
    mt_q98 = float(np.quantile(percent_mt, 0.98))
    return {
        "min_nFeature_RNA": int(max(100, min(300, math.floor(feature_q01)))),
        "min_nCount_RNA": int(max(200, min(500, math.floor(count_q01)))),
        "max_percent_mt": float(min(30, max(20, math.ceil(mt_q98)))),
        "threshold_basis": "all_source_cells_q01_q01_q98_bounded_identity_independent",
        "all_feature_q01": feature_q01,
        "all_count_q01": count_q01,
        "all_percent_mt_q98": mt_q98,
    }


def summarize(values: np.ndarray, prefix: str) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            f"{prefix}_median": float("nan"),
            f"{prefix}_q1": float("nan"),
            f"{prefix}_q3": float("nan"),
            f"{prefix}_iqr": float("nan"),
        }
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return {
        f"{prefix}_median": float(median),
        f"{prefix}_q1": float(q1),
        f"{prefix}_q3": float(q3),
        f"{prefix}_iqr": float(q3 - q1),
    }


def fraction(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def marker_presence(matrix: csr_matrix, gene_map: dict[str, np.ndarray], marker: str) -> np.ndarray:
    rows = gene_map.get(marker)
    if rows is None:
        return np.zeros(matrix.shape[1], dtype=bool)
    return np.asarray(matrix[rows, :].sum(axis=0)).ravel() > 0


def load_initial_passes(project_root: Path, sample: SampleSpec, barcodes: list[str]) -> np.ndarray:
    path = project_root / "processed_data/task002/cell_qc" / sample.dataset / f"{sample.sample_id}.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Initial Task 002 cell QC is missing: {path}")
    cell_qc = pd.read_csv(path, sep="\t", compression="gzip", dtype={"barcode": str})
    if len(cell_qc) != len(barcodes) or set(cell_qc["barcode"]) != set(barcodes):
        raise ValueError(f"Initial QC barcode identity mismatch for {sample.dataset}/{sample.sample_id}")
    indexed = cell_qc.set_index("barcode").loc[barcodes]
    return indexed["passes_qc"].astype(bool).to_numpy()


def make_review_page(
    pdf: PdfPages,
    sample: SampleSpec,
    qc: pd.DataFrame,
    thresholds: dict[str, object],
    summary: dict[str, object],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    initial_pass = qc["passes_initial_qc"].to_numpy(dtype=bool)
    corrected_pass = qc["passes_qc"].to_numpy(dtype=bool)
    rescued = corrected_pass & ~initial_pass
    initial_only = initial_pass & ~corrected_pass
    colors = np.full(len(qc), "#bdbdbd", dtype=object)
    colors[initial_pass & corrected_pass] = "#2166ac"
    colors[rescued] = "#1b9e77"
    colors[initial_only] = "#d95f02"
    axes[0, 0].scatter(
        np.log10(qc["nFeature_RNA"].to_numpy() + 1),
        np.log10(qc["nCount_RNA"].to_numpy() + 1),
        c=colors,
        s=4,
        alpha=0.35,
        linewidths=0,
    )
    axes[0, 0].axvline(np.log10(float(thresholds["min_nFeature_RNA"]) + 1), color="black", ls="--", lw=0.8)
    axes[0, 0].axhline(np.log10(float(thresholds["min_nCount_RNA"]) + 1), color="black", ls="--", lw=0.8)
    axes[0, 0].set_xlabel("log10(nFeature_RNA + 1)")
    axes[0, 0].set_ylabel("log10(nCount_RNA + 1)")
    axes[0, 0].set_title("Corrected QC comparison")

    for mask, label, color in [
        (np.ones(len(qc), dtype=bool), "source", "#777777"),
        (initial_pass, "initial pass", "#d95f02"),
        (corrected_pass, "corrected pass", "#2166ac"),
        (rescued, "rescued", "#1b9e77"),
    ]:
        if mask.any():
            axes[0, 1].hist(
                np.log10(qc.loc[mask, "nFeature_RNA"] + 1),
                bins=40,
                histtype="step",
                lw=1.2,
                label=label,
                color=color,
            )
    axes[0, 1].axvline(np.log10(float(thresholds["min_nFeature_RNA"]) + 1), color="black", ls="--", lw=0.8)
    axes[0, 1].set_xlabel("log10(nFeature_RNA + 1)")
    axes[0, 1].set_ylabel("cells")
    axes[0, 1].set_title("Detected genes")
    axes[0, 1].legend(frameon=False, fontsize=7)

    box_data = []
    labels = []
    for metric in ["nFeature_RNA", "nCount_RNA", "percent.mt"]:
        values = qc.loc[rescued, metric].to_numpy(dtype=float)
        if values.size:
            box_data.append(values)
            labels.append(f"{metric}\nrescued")
        values = qc.loc[corrected_pass, metric].to_numpy(dtype=float)
        if values.size:
            box_data.append(values)
            labels.append(f"{metric}\ncorrected")
    if box_data:
        axes[1, 0].boxplot(box_data, labels=labels, showfliers=False)
    axes[1, 0].set_title("Rescued versus corrected-pass distributions")
    axes[1, 0].tick_params(axis="x", labelrotation=70, labelsize=7)
    axes[1, 0].set_ylabel("value; percent.mt is %")

    labels = ["source", "initial", "corrected", "HC initial", "HC corrected", "broad initial", "broad corrected"]
    values = [
        int(summary["n_cells_source"]),
        int(summary["n_cells_initial_qc"]),
        int(summary["n_cells_corrected_qc"]),
        int(summary["n_high_confidence_after_initial"]),
        int(summary["n_high_confidence_after_corrected"]),
        int(summary["n_broad_after_initial"]),
        int(summary["n_broad_after_corrected"]),
    ]
    axes[1, 1].bar(labels, values, color=["#777777", "#d95f02", "#2166ac", "#f4a582", "#67a9cf", "#fdbb84", "#74c476"])
    axes[1, 1].set_title("Cell counts and audit populations")
    axes[1, 1].tick_params(axis="x", labelrotation=75, labelsize=7)
    axes[1, 1].set_ylabel("cells")
    for index, value in enumerate(values):
        axes[1, 1].text(index, value, str(value), ha="center", va="bottom", fontsize=7)

    flag = "FLAGGED" if summary["flag_review"] else "no mandatory flag"
    fig.suptitle(
        f"Task 002 corrected QC | {sample.dataset} | {sample.sample_id} | {sample.tissue} | {sample.patient_id}\n"
        f"min genes={thresholds['min_nFeature_RNA']}; min counts={thresholds['min_nCount_RNA']}; "
        f"max mt={thresholds['max_percent_mt']:.0f}% | {flag}",
        fontsize=11,
    )
    pdf.savefig(fig)
    plt.close(fig)


def process_sample(
    project_root: Path,
    out_root: Path,
    sample: SampleSpec,
    data: MatrixData,
    structure_format: str,
    stage_source: str,
    audit_rows: list[dict[str, object]],
    threshold_rows: list[dict[str, object]],
    comparison_rows: list[dict[str, object]],
    rescued_rows: list[dict[str, object]],
    structure_rows: list[dict[str, object]],
    stage_rows: list[dict[str, object]],
    review_rows: list[tuple[SampleSpec, pd.DataFrame, dict[str, object], dict[str, object]]],
) -> None:
    matrix = data.matrix.tocsr()
    n_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    n_features = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(float)
    gene_map = initial.gene_mask_map(data.gene_names)
    mt_keys = [key for key in gene_map if key.startswith("MT-") or key.startswith("MT.")]
    mt_rows = np.concatenate([gene_map[key] for key in mt_keys]) if mt_keys else np.array([], dtype=np.int64)
    mt_counts = np.asarray(matrix[mt_rows, :].sum(axis=0)).ravel() if mt_rows.size else np.zeros(matrix.shape[1])
    percent_mt = np.divide(mt_counts * 100.0, n_counts, out=np.zeros_like(n_counts), where=n_counts > 0)

    core_hits, support_hits, missing_markers = initial.marker_hits(matrix, data.gene_names)
    s100a8 = marker_presence(matrix, gene_map, "S100A8")
    s100a9 = marker_presence(matrix, gene_map, "S100A9")
    high_confidence = (core_hits >= 2) | ((core_hits >= 1) & (support_hits >= 2))
    broad = high_confidence | ((support_hits >= 4) & s100a8 & s100a9)

    thresholds = corrected_thresholds(n_features, n_counts, percent_mt)
    corrected_pass = (
        (n_features >= thresholds["min_nFeature_RNA"])
        & (n_counts >= thresholds["min_nCount_RNA"])
        & (percent_mt <= thresholds["max_percent_mt"])
    )
    initial_pass = load_initial_passes(project_root, sample, data.barcodes)
    rescued = corrected_pass & ~initial_pass
    initial_only = initial_pass & ~corrected_pass
    count_q995 = np.quantile(n_counts, 0.995) if len(n_counts) >= 20 else np.max(n_counts)
    feature_q995 = np.quantile(n_features, 0.995) if len(n_features) >= 20 else np.max(n_features)
    doublet_proxy = (n_counts >= count_q995) & (n_features >= feature_q995)

    qc = pd.DataFrame(
        {
            "barcode": data.barcodes,
            "dataset": sample.dataset,
            "sample_id": sample.sample_id,
            "patient_id": sample.patient_id,
            "tissue": sample.tissue,
            "nCount_RNA": n_counts,
            "nFeature_RNA": n_features,
            "percent.mt": percent_mt,
            "core_marker_hits": core_hits,
            "support_marker_hits": support_hits,
            "candidate_neutrophil": high_confidence,
            "high_confidence_neutrophil": high_confidence,
            "broad_granulocyte_like": broad,
            "passes_initial_qc": initial_pass,
            "passes_corrected_qc": corrected_pass,
            "passes_qc": corrected_pass,
            "doublet_risk_proxy": doublet_proxy,
            "min_nFeature_RNA": thresholds["min_nFeature_RNA"],
            "min_nCount_RNA": thresholds["min_nCount_RNA"],
            "max_percent_mt": thresholds["max_percent_mt"],
        }
    )
    qc["removal_reason"] = np.where(
        corrected_pass,
        "pass",
        np.where(
            n_features < thresholds["min_nFeature_RNA"],
            "low_nFeature_RNA",
            np.where(n_counts < thresholds["min_nCount_RNA"], "low_nCount_RNA", "high_percent.mt"),
        ),
    )

    n_source = len(qc)
    n_initial = int(initial_pass.sum())
    n_corrected = int(corrected_pass.sum())
    n_rescued = int(rescued.sum())
    n_initial_only = int(initial_only.sum())
    rescued_broad = rescued & broad
    rescued_coreless_broad = rescued_broad & (core_hits == 0)
    flag_reasons: list[str] = []
    hc_retention = fraction(int((high_confidence & corrected_pass).sum()), int(high_confidence.sum()))
    broad_retention = fraction(int((broad & corrected_pass).sum()), int(broad.sum()))
    mt_rescued_median = float(np.median(percent_mt[rescued])) if n_rescued else float("nan")
    if np.isfinite(hc_retention) and hc_retention < 0.90:
        flag_reasons.append("corrected_high_confidence_retention_lt_90pct")
    if (n_corrected - n_initial) / n_source > 0.20:
        flag_reasons.append("corrected_total_retention_gain_gt_20_percentage_points")
    if np.isfinite(mt_rescued_median) and mt_rescued_median > 25:
        flag_reasons.append("rescued_median_percent_mt_gt_25")
    if n_rescued and int(rescued_broad.sum()) and rescued_coreless_broad.sum() / rescued_broad.sum() > 0.50:
        flag_reasons.append("more_than_half_rescued_broad_cells_without_core_marker")
    if n_corrected > n_source:
        flag_reasons.append("corrected_retained_exceeds_source")
    flag_review = bool(flag_reasons)

    summary: dict[str, object] = {
        "n_cells_source": n_source,
        "n_cells_initial_qc": n_initial,
        "n_cells_corrected_qc": n_corrected,
        "n_cells_rescued": n_rescued,
        "n_cells_initial_only": n_initial_only,
        "n_high_confidence_before": int(high_confidence.sum()),
        "n_high_confidence_after_initial": int((high_confidence & initial_pass).sum()),
        "n_high_confidence_after_corrected": int((high_confidence & corrected_pass).sum()),
        "n_broad_before": int(broad.sum()),
        "n_broad_after_initial": int((broad & initial_pass).sum()),
        "n_broad_after_corrected": int((broad & corrected_pass).sum()),
        "high_confidence_retention_corrected": hc_retention,
        "broad_retention_corrected": broad_retention,
        "flag_review": flag_review,
        "flag_reasons": ";".join(flag_reasons) if flag_reasons else "none",
    }

    cell_qc_path = out_root / "cell_qc" / sample.dataset / f"{sample.sample_id}.tsv.gz"
    cell_qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(cell_qc_path, sep="\t", index=False, compression="gzip")
    kept_indices = np.flatnonzero(corrected_pass)
    staged_data = MatrixData(
        matrix=matrix[:, kept_indices].tocsr(),
        gene_ids=list(data.gene_ids),
        gene_names=list(data.gene_names),
        feature_types=list(data.feature_types),
        barcodes=[data.barcodes[index] for index in kept_indices],
        source_member=data.source_member,
        structure_note=data.structure_note + "; corrected stage contains passes_corrected_qc cells",
    )
    stage_dir = out_root / "staged_matrices" / sample.dataset / sample.sample_id
    matrix_path, features_path, barcodes_path = initial.matrix_to_stage(staged_data, stage_dir)
    stage_rows.append(
        {
            **sample.__dict__,
            "matrix_path": str(matrix_path),
            "features_path": str(features_path),
            "barcodes_path": str(barcodes_path),
            "cell_qc_path": str(cell_qc_path),
            "n_features": matrix.shape[0],
            "n_cells_source": n_source,
            "n_cells_before_qc": n_source,
            "n_cells_staged_after_qc": n_corrected,
            "n_cells_after_corrected_qc": n_corrected,
            "source_member": stage_source,
        }
    )
    structure_rows.append(
        {
            **sample.__dict__,
            "structure_format": structure_format,
            "source_member": stage_source,
            "n_features": matrix.shape[0],
            "n_cells": n_source,
            "feature_ids_recovered": bool(len(data.gene_ids) == matrix.shape[0] and all(data.gene_ids)),
            "feature_names_recovered": bool(len(data.gene_names) == matrix.shape[0] and all(data.gene_names)),
            "barcodes_recovered": bool(len(data.barcodes) == matrix.shape[1] and all(data.barcodes)),
            "orientation": "features_by_barcodes",
            "structure_check": "PASS",
            "structure_note": data.structure_note,
            "missing_marker_genes": ";".join(missing_markers),
        }
    )

    threshold_rows.append({**sample.__dict__, **thresholds, **summary, "missing_marker_genes": ";".join(missing_markers)})
    comparison_rows.append({**sample.__dict__, **summary, "source_member": stage_source})
    rescued_row = {**sample.__dict__, "n_rescued_cells": n_rescued}
    for metric in ["nFeature_RNA", "nCount_RNA", "percent.mt"]:
        rescued_row.update(summarize(qc.loc[rescued, metric].to_numpy(), f"rescued_{metric}"))
    rescued_row.update(
        {
            "rescued_fraction_with_core_marker": fraction(int((rescued & (core_hits >= 1)).sum()), n_rescued),
            "rescued_fraction_high_confidence": fraction(int((rescued & high_confidence).sum()), n_rescued),
            "rescued_fraction_broad_granulocyte_like": fraction(int(rescued_broad.sum()), n_rescued),
            "rescued_fraction_s100a8_s100a9_only": fraction(
                int((rescued & s100a8 & s100a9 & (core_hits == 0)).sum()), n_rescued
            ),
            "rescued_fraction_broad_without_core": fraction(int(rescued_coreless_broad.sum()), int(rescued_broad.sum())),
        }
    )
    rescued_rows.append(rescued_row)
    audit_row = {**sample.__dict__, **summary}
    audit_row.update(
        {
            "n_high_confidence_neutrophils_before_qc": int(high_confidence.sum()),
            "n_high_confidence_neutrophils_after_initial_qc": int((high_confidence & initial_pass).sum()),
            "n_high_confidence_neutrophils_after_corrected_qc": int((high_confidence & corrected_pass).sum()),
            "corrected_high_confidence_neutrophil_retention_fraction": hc_retention,
            "n_broad_granulocyte_like_before_qc": int(broad.sum()),
            "n_broad_granulocyte_like_after_initial_qc": int((broad & initial_pass).sum()),
            "n_broad_granulocyte_like_after_corrected_qc": int((broad & corrected_pass).sum()),
            "corrected_broad_granulocyte_retention_fraction": broad_retention,
            "n_doublet_proxy_review": int(doublet_proxy.sum()),
            "missing_marker_genes": ";".join(missing_markers),
        }
    )
    audit_rows.append(audit_row)
    review_rows.append((sample, qc, thresholds, summary))


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    out_root = args.out_root.resolve()
    results_root = project_root / "results"
    figures_root = project_root / "figures"
    for path in [out_root / "staged_matrices", out_root / "cell_qc", results_root, figures_root]:
        path.mkdir(parents=True, exist_ok=True)
    samples = initial.read_manifest(args.manifest)
    by_dataset: dict[str, list[SampleSpec]] = defaultdict(list)
    for sample in samples:
        by_dataset[sample.dataset].append(sample)

    audit_rows: list[dict[str, object]] = []
    threshold_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    rescued_rows: list[dict[str, object]] = []
    structure_rows: list[dict[str, object]] = []
    stage_rows: list[dict[str, object]] = []
    review_rows: list[tuple[SampleSpec, pd.DataFrame, dict[str, object], dict[str, object]]] = []

    gse149_count_cells: list[str] | None = None
    gse149_sample_info: dict[str, dict[str, object]] | None = None
    gse149_parts: dict[str, list[csr_matrix]] = {}
    gse149_gene_names: list[str] = []
    if "GSE149614" in by_dataset:
        _, gse149_count_cells, gse149_sample_info = initial.load_gse149614_inputs(project_root)
        gse149_parts = {sample_id: [] for sample_id in gse149_sample_info}
        for gene_chunk, matrix_parts in initial.iterate_gse149614_chunks(
            project_root, gse149_count_cells, gse149_sample_info
        ):
            gse149_gene_names.extend(gene_chunk)
            for sample_id, matrix_part in matrix_parts.items():
                gse149_parts[sample_id].append(matrix_part)

    with tempfile.TemporaryDirectory(prefix="task002_corrected_h5_") as temp_name:
        temp_dir = Path(temp_name)
        for dataset, dataset_samples in by_dataset.items():
            if dataset == "GSE149614":
                assert gse149_count_cells is not None and gse149_sample_info is not None
                for sample in dataset_samples:
                    sample_key = sample.sample_title if sample.sample_title in gse149_sample_info else sample.sample_id
                    if sample_key not in gse149_sample_info:
                        raise KeyError(f"GSE149614 sample mapping missing for {sample.sample_id}")
                    info = gse149_sample_info[sample_key]
                    matrix = vstack(gse149_parts[sample_key], format="csr").astype(np.int32)
                    data = MatrixData(
                        matrix=matrix,
                        gene_ids=list(gse149_gene_names),
                        gene_names=list(gse149_gene_names),
                        feature_types=["Gene Expression"] * len(gse149_gene_names),
                        barcodes=list(info["barcodes"]),
                        source_member=(
                            "GSE149614_HCC.scRNAseq.S71915.count.txt.gz;"
                            "GSE149614_HCC.metadata.updated.txt.gz;"
                            f"site={sample.tissue}"
                        ),
                        structure_note="chunked combined gene-by-cell count text; primary Tumor/Adjacent cells only",
                    )
                    process_sample(
                        project_root, out_root, sample, data, "combined_count_text_plus_metadata", data.source_member,
                        audit_rows, threshold_rows, comparison_rows, rescued_rows, structure_rows, stage_rows, review_rows,
                    )
                    del matrix, data
                continue

            archive_path = project_root / "raw_data" / dataset / f"{dataset}_RAW.tar"
            with tarfile.open(archive_path, mode="r:") as archive:
                for sample in dataset_samples:
                    if dataset in {"GSE282701", "GSE299340"}:
                        data = initial.load_10x_matrix_from_tar(
                            archive, sample.sample_id, ("_matrix.mtx.gz", "_features.tsv.gz", "_barcodes.tsv.gz")
                        )
                        structure_format = "10x_matrix_market_gz"
                    elif dataset == "GSE242889":
                        data = initial.load_gse242889_sample(archive, sample.sample_id)
                        structure_format = "nested_tar_matrix_market"
                    elif dataset == "GSE326201":
                        data = initial.load_gse326201_h5(archive, sample.sample_id, temp_dir)
                        structure_format = "cellranger_filtered_h5"
                    else:
                        raise ValueError(f"Unexpected dataset: {dataset}")
                    process_sample(
                        project_root, out_root, sample, data, structure_format, data.source_member,
                        audit_rows, threshold_rows, comparison_rows, rescued_rows, structure_rows, stage_rows, review_rows,
                    )

    audit_df = pd.DataFrame(audit_rows).sort_values(["dataset", "sample_id"])
    threshold_df = pd.DataFrame(threshold_rows).sort_values(["dataset", "sample_id"])
    comparison_df = pd.DataFrame(comparison_rows).sort_values(["dataset", "sample_id"])
    rescued_df = pd.DataFrame(rescued_rows).sort_values(["dataset", "sample_id"])
    structure_df = pd.DataFrame(structure_rows).sort_values(["dataset", "sample_id"])
    stage_df = pd.DataFrame(stage_rows).sort_values(["dataset", "sample_id"])
    audit_df.to_csv(results_root / "task002_corrected_neutrophil_retention_audit.csv", index=False)
    threshold_df.to_csv(results_root / "task002_corrected_qc_thresholds_by_sample.csv", index=False)
    comparison_df.to_csv(results_root / "task002_initial_vs_corrected_qc.csv", index=False)
    rescued_df.to_csv(results_root / "task002_rescued_cells_summary.csv", index=False)
    structure_df.to_csv(results_root / "task002_corrected_input_structure_audit.csv", index=False)
    stage_df.to_csv(out_root / "task002_corrected_stage_manifest.csv", index=False)
    with (figures_root / "task002_corrected_qc_review.pdf").open("wb") as handle:
        with PdfPages(handle) as pdf:
            for sample, qc, thresholds, summary in review_rows:
                make_review_page(pdf, sample, qc, thresholds, summary)
    run_meta = {
        "task": "task_002_corrected",
        "script": "task002_corrected_preprocess_qc.py",
        "approved_datasets": sorted(APPROVED_DATASETS),
        "excluded_datasets": ["GSE202642", "GSE290298"],
        "n_samples": len(samples),
        "threshold_rule": "all-source-cell 1st percentiles for nFeature/nCount bounded to 100-300 and 200-500; all-source-cell 98th percentile percent.mt bounded to 20-30",
        "high_confidence_candidate_definition": "core_hits>=2 OR (core_hits>=1 AND support_hits>=2); audit only",
        "broad_granulocyte_definition": "high_confidence OR (support_hits>=4 AND S100A8 and S100A9 expressed); sensitivity audit only",
        "doublet_evaluation": "high nCount_RNA and nFeature_RNA 99.5th-percentile proxy retained as review flag; no hard removal",
        "initial_comparison_source": "processed_data/task002/cell_qc/*",
        "raw_counts_preserved": True,
        "corrected_objects_root": "objects/task002_corrected_seurat",
    }
    (results_root / "task002_corrected_run_metadata.json").write_text(json.dumps(run_meta, indent=2) + "\n", encoding="utf-8")
    print(f"Processed corrected QC for {len(samples)} samples across {len(by_dataset)} datasets")
    print(f"Audit: {results_root / 'task002_corrected_neutrophil_retention_audit.csv'}")
    print(f"Comparison: {results_root / 'task002_initial_vs_corrected_qc.csv'}")


if __name__ == "__main__":
    main()
