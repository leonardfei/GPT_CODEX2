#!/usr/bin/env python3
"""Prepare the CRA002308 and in_house 10x count-matrix cohorts for Task 003.

This script is intentionally limited to uploaded cell-called count matrices. It
does not infer cell identities for QC thresholds, and it never edits raw input
files. It writes server-side staged matrices, per-cell QC, and small audit
tables. The large staged matrices and Seurat objects remain on the server.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


SCRIPT_DIR = Path(__file__).resolve().parent
INITIAL_PATH = SCRIPT_DIR / "task002_preprocess_qc.py"
if not INITIAL_PATH.exists():
    INITIAL_PATH = SCRIPT_DIR.parent / "task002_preprocess_qc.py"
spec = importlib.util.spec_from_file_location("task002_initial", INITIAL_PATH)
if spec is None or spec.loader is None:
    raise ImportError(f"Unable to import {INITIAL_PATH}")
initial = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = initial
spec.loader.exec_module(initial)

CORE_MARKERS = initial.CORE_MARKERS
SUPPORT_MARKERS = initial.SUPPORT_MARKERS
MatrixData = initial.MatrixData


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-root", required=True, type=Path)
    p.add_argument("--results-root", required=True, type=Path)
    return p.parse_args()


def corrected_thresholds(n_features: np.ndarray, n_counts: np.ndarray, percent_mt: np.ndarray) -> dict[str, object]:
    q01_features = float(np.quantile(n_features, 0.01))
    q01_counts = float(np.quantile(n_counts, 0.01))
    q98_mt = float(np.quantile(percent_mt, 0.98))
    return {
        "min_nFeature_RNA": int(max(100, min(300, np.floor(q01_features)))),
        "min_nCount_RNA": int(max(200, min(500, np.floor(q01_counts)))),
        "max_percent_mt": float(min(30, max(20, np.ceil(q98_mt)))),
        "threshold_basis": "all_source_cells_q01_q01_q98_bounded_identity_independent",
        "all_feature_q01": q01_features,
        "all_count_q01": q01_counts,
        "all_percent_mt_q98": q98_mt,
    }


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                return digest.hexdigest()
            digest.update(chunk)


def gzip_ok(path: Path) -> str:
    try:
        result = subprocess.run(["gzip", "-t", str(path)], capture_output=True, text=True)
    except OSError as exc:
        return f"ERROR:{exc}"
    return "PASS" if result.returncode == 0 else f"FAIL:{(result.stderr or result.stdout).strip()[:200]}"


def classify_file(path: Path) -> str:
    name = path.name.lower()
    if name.endswith("matrix.mtx.gz"):
        return "10x Matrix Market count matrix"
    if name.endswith("features.tsv.gz"):
        return "10x feature annotation"
    if name.endswith("barcodes.tsv.gz"):
        return "10x cell barcode list"
    if name.endswith(".rds"):
        return "RDS author object"
    if name.endswith(".docx"):
        return "Word supplementary methods/document"
    return "other uploaded file"


def sample_from_path(cohort: str, path: Path) -> str:
    if cohort in {"CRA002308", "in_house"}:
        return path.parent.name
    return ""


def build_file_inventory(raw_root: Path, results_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort in ["CRA002308", "nature_xue", "in_house"]:
        cohort_root = raw_root / cohort
        for path in sorted(p for p in cohort_root.rglob("*") if p.is_file()):
            is_primary = (
                cohort in {"CRA002308", "in_house"}
                and path.name in {"matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"}
            ) or (cohort == "nature_xue" and path.name == "seu_A124.counts_anno.rds") or path.suffix.lower() == ".docx"
            rows.append(
                {
                    "dataset": cohort,
                    "relative_path": str(path.relative_to(raw_root)),
                    "absolute_path": str(path),
                    "filename": path.name,
                    "sample_id_from_path": sample_from_path(cohort, path),
                    "size_bytes": path.stat().st_size,
                    "extension": path.suffix.lower(),
                    "file_type": classify_file(path),
                    "primary_file_used": bool(is_primary),
                    "sha256": sha256(path) if is_primary else "",
                    "gzip_integrity": gzip_ok(path) if path.name.endswith(".gz") else "not_applicable",
                }
            )
    result = pd.DataFrame(rows).sort_values(["dataset", "relative_path"])
    results_root.mkdir(parents=True, exist_ok=True)
    result.to_csv(results_root / "task003_extension_file_inventory.csv", index=False)
    return result


def read_10x_dir(sample_dir: Path) -> tuple[MatrixData, str]:
    matrix_path = sample_dir / "matrix.mtx.gz"
    features_path = sample_dir / "features.tsv.gz"
    barcodes_path = sample_dir / "barcodes.tsv.gz"
    with matrix_path.open("rb") as handle:
        matrix = initial.read_mtx_stream(handle, compressed=True)
    with gzip.open(features_path, "rb") as handle:
        gene_ids, gene_names, feature_types = initial.parse_features(handle)
    with gzip.open(barcodes_path, "rb") as handle:
        barcodes = initial.parse_barcodes(handle)
    if matrix.shape != (len(gene_names), len(barcodes)):
        raise ValueError(f"Dimension mismatch in {sample_dir}: {matrix.shape}, {len(gene_names)}, {len(barcodes)}")
    return MatrixData(
        matrix=matrix.tocsr().astype(np.int32),
        gene_ids=gene_ids,
        gene_names=gene_names,
        feature_types=feature_types,
        barcodes=barcodes,
        source_member=str(sample_dir),
        structure_note="uploaded 10x Matrix Market cell-called count matrix; dimensions/features/barcodes validated",
    ), "10x_matrix_market_gz_cell_called"


def mt_percent(matrix: csr_matrix, gene_names: list[str]) -> np.ndarray:
    gene_map = initial.gene_mask_map(gene_names)
    mt_keys = [key for key in gene_map if key.startswith("MT-") or key.startswith("MT.")]
    rows = np.concatenate([gene_map[key] for key in mt_keys]) if mt_keys else np.array([], dtype=np.int64)
    mt_counts = np.asarray(matrix[rows, :].sum(axis=0)).ravel() if rows.size else np.zeros(matrix.shape[1])
    counts = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    return np.divide(mt_counts * 100.0, counts, out=np.zeros_like(counts), where=counts > 0)


def marker_presence(matrix: csr_matrix, gene_map: dict[str, np.ndarray], marker: str) -> np.ndarray:
    rows = gene_map.get(marker)
    if rows is None:
        return np.zeros(matrix.shape[1], dtype=bool)
    return np.asarray(matrix[rows, :].sum(axis=0)).ravel() > 0


def metadata(dataset: str, sample: str) -> dict[str, str]:
    if dataset == "CRA002308":
        match = re.fullmatch(r"([NT])(\d{2})", sample)
        if not match:
            raise ValueError(f"Unexpected CRA002308 sample name: {sample}")
        tissue_code, number = match.groups()
        tissue = "Tumor" if tissue_code == "T" else "Adjacent"
        patient = f"CRA002308_P{int(number):02d}"
        return {
            "dataset": dataset, "sample_id": sample, "patient_id": patient,
            "tissue": tissue, "paired_status": "paired", "paired_id": patient,
            "etiology": "unknown", "MVI": "unknown", "platform": "10x_scRNA_cell_called_matrix",
            "selection_strategy": "flow-sorted live nucleated cells after doublet exclusion; no lineage-specific immune enrichment documented; paired within-cohort abundance sensitivity analysis allowed",
            "matrix_type": "cell_called_count_matrix", "counts_available": "TRUE",
            "qc_provenance": "task003_corrected_task002_identity_independent",
            "abundance_eligible": "CONDITIONAL", "analysis_inclusion": "YES",
        }
    if dataset == "in_house":
        match = re.fullmatch(r"(YJCA|YJP)(\d{2})", sample)
        if not match:
            raise ValueError(f"Unexpected in_house sample name: {sample}")
        tissue_code, number = match.groups()
        tissue = "Tumor" if tissue_code == "YJCA" else "Adjacent"
        patient = f"in_house_P{int(number):02d}"
        return {
            "dataset": dataset, "sample_id": sample, "patient_id": patient,
            "tissue": tissue, "paired_status": "paired", "paired_id": patient,
            "etiology": "unknown", "MVI": "unknown", "platform": "10x_scRNA_cell_called_matrix",
            "selection_strategy": "unknown from uploaded files; provenance review required before abundance use",
            "matrix_type": "cell_called_count_matrix", "counts_available": "TRUE",
            "qc_provenance": "task003_corrected_task002_identity_independent",
            "abundance_eligible": "CONDITIONAL", "analysis_inclusion": "YES",
        }
    raise ValueError(dataset)


def process_sample(dataset: str, sample_dir: Path, out_root: Path) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    data, structure_format = read_10x_dir(sample_dir)
    meta = metadata(dataset, sample_dir.name)
    matrix = data.matrix
    n_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    n_features = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(float)
    percent_mt = mt_percent(matrix, data.gene_names)
    core_hits, support_hits, missing = initial.marker_hits(matrix, data.gene_names)
    gene_map = initial.gene_mask_map(data.gene_names)
    s100a8 = marker_presence(matrix, gene_map, "S100A8")
    s100a9 = marker_presence(matrix, gene_map, "S100A9")
    high_confidence = (core_hits >= 2) | ((core_hits >= 1) & (support_hits >= 2))
    broad = high_confidence | ((support_hits >= 4) & s100a8 & s100a9)
    thresholds = corrected_thresholds(n_features, n_counts, percent_mt)
    passes = (n_features >= thresholds["min_nFeature_RNA"]) & (n_counts >= thresholds["min_nCount_RNA"]) & (percent_mt <= thresholds["max_percent_mt"])
    count_q995 = np.quantile(n_counts, 0.995) if len(n_counts) >= 20 else np.max(n_counts)
    feature_q995 = np.quantile(n_features, 0.995) if len(n_features) >= 20 else np.max(n_features)
    doublet_proxy = (n_counts >= count_q995) & (n_features >= feature_q995)
    global_barcodes = [f"{dataset}__{sample_dir.name}__{barcode}" for barcode in data.barcodes]
    qc = pd.DataFrame({
        "barcode": global_barcodes, "source_barcode": data.barcodes,
        **{key: meta[key] for key in ["dataset", "sample_id", "patient_id", "tissue", "paired_status", "paired_id", "etiology", "MVI", "platform", "selection_strategy", "matrix_type", "counts_available", "qc_provenance", "abundance_eligible", "analysis_inclusion"]},
        "nCount_RNA": n_counts, "nFeature_RNA": n_features, "percent.mt": percent_mt,
        "core_marker_hits": core_hits, "support_marker_hits": support_hits,
        "candidate_neutrophil": high_confidence, "high_confidence_neutrophil": high_confidence,
        "broad_granulocyte_like": broad, "passes_initial_qc": "not_applicable",
        "passes_corrected_qc": passes, "passes_qc": passes, "qc_status": np.where(passes, "PASS_CORRECTED_QC", "FAIL_CORRECTED_QC"),
        "doublet_risk_proxy": doublet_proxy,
        "min_nFeature_RNA": thresholds["min_nFeature_RNA"], "min_nCount_RNA": thresholds["min_nCount_RNA"], "max_percent_mt": thresholds["max_percent_mt"],
    })
    qc["removal_reason"] = np.where(passes, "pass", np.where(n_features < thresholds["min_nFeature_RNA"], "low_nFeature_RNA", np.where(n_counts < thresholds["min_nCount_RNA"], "low_nCount_RNA", "high_percent.mt")))
    qc_path = out_root / "cell_qc" / dataset / f"{sample_dir.name}.tsv.gz"
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(qc_path, sep="\t", index=False, compression="gzip")
    kept = np.flatnonzero(passes)
    staged = MatrixData(matrix=matrix[:, kept].tocsr(), gene_ids=data.gene_ids, gene_names=data.gene_names, feature_types=data.feature_types, barcodes=[global_barcodes[i] for i in kept], source_member=str(sample_dir), structure_note=data.structure_note + "; staged corrected-QC pass cells with globally unique IDs")
    stage_dir = out_root / "staged_matrices" / dataset / sample_dir.name
    matrix_path, features_path, barcodes_path = initial.matrix_to_stage(staged, stage_dir)
    stage_row = {**meta, "matrix_path": str(matrix_path), "features_path": str(features_path), "barcodes_path": str(barcodes_path), "cell_qc_path": str(qc_path), "source_member": str(sample_dir), "structure_format": structure_format, "n_features": matrix.shape[0], "n_cells_source": len(qc), "n_cells_staged_after_qc": int(passes.sum()), "missing_marker_genes": ";".join(missing)}
    broad_retention = float((broad & passes).sum() / broad.sum()) if broad.sum() else float("nan")
    hc_retention = float((high_confidence & passes).sum() / high_confidence.sum()) if high_confidence.sum() else float("nan")
    audit = {**meta, "source_member": str(sample_dir), "n_features": matrix.shape[0], "n_cells_source": len(qc), "n_cells_after_corrected_qc": int(passes.sum()), "retention_fraction": float(passes.mean()), "n_high_confidence_neutrophils_before_qc": int(high_confidence.sum()), "n_high_confidence_neutrophils_after_corrected_qc": int((high_confidence & passes).sum()), "corrected_high_confidence_neutrophil_retention_fraction": hc_retention, "n_broad_granulocyte_like_before_qc": int(broad.sum()), "n_broad_granulocyte_like_after_corrected_qc": int((broad & passes).sum()), "corrected_broad_granulocyte_retention_fraction": broad_retention, "n_doublet_proxy_review": int(doublet_proxy.sum()), **thresholds, "flag_review": bool((np.isfinite(hc_retention) and hc_retention < 0.90) or (np.isfinite(broad_retention) and broad_retention < 0.90)), "qc_status": "PASS_WITH_REVIEW_FLAG" if ((np.isfinite(hc_retention) and hc_retention < 0.90) or (np.isfinite(broad_retention) and broad_retention < 0.90)) else "PASS", "missing_marker_genes": ";".join(missing)}
    manifest = {**meta, "source_member": str(sample_dir), "n_features": matrix.shape[0], "n_cells_source": len(qc), "n_cells_after_corrected_qc": int(passes.sum()), "qc_status": "PASS_WITH_REVIEW_FLAG" if audit["flag_review"] else "PASS", "prepared_stage": str(stage_dir)}
    del matrix, data, staged, qc
    return stage_row, audit, manifest


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    raw = project / "raw_data"
    out_root = args.out_root.resolve()
    results_root = args.results_root.resolve()
    inventory = build_file_inventory(raw, results_root)
    audit_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    stage_rows: list[dict[str, object]] = []
    for dataset in ["CRA002308", "in_house"]:
        sample_dirs = sorted(p for p in (raw / dataset).iterdir() if p.is_dir())
        expected = 14 if dataset == "CRA002308" else 20
        if len(sample_dirs) != expected:
            raise RuntimeError(f"Expected {expected} sample directories for {dataset}, found {len(sample_dirs)}")
        for sample_dir in sample_dirs:
            print(f"PROCESS {dataset}/{sample_dir.name}", flush=True)
            stage_row, audit, manifest = process_sample(dataset, sample_dir, out_root)
            stage_rows.append(stage_row)
            audit_rows.append(audit)
            manifest_rows.append(manifest)
    pd.DataFrame(audit_rows).sort_values(["dataset", "sample_id"]).to_csv(results_root / "task003_extension_qc_audit.csv", index=False)
    pd.DataFrame(manifest_rows).sort_values(["dataset", "sample_id"]).to_csv(results_root / "task003_extension_sample_manifest_matrix_cohorts.csv", index=False)
    pd.DataFrame(stage_rows).sort_values(["dataset", "sample_id"]).to_csv(out_root / "task003_extension_stage_manifest.csv", index=False)
    with (out_root / "task003_extension_run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump({"task": "task_003", "script": "task003_prepare_extensions.py", "datasets": ["CRA002308", "in_house"], "qc_rule": "corrected Task 002 all-source q01/q01/q98 bounded identity-independent", "raw_counts_preserved": True, "inventory_rows": len(inventory), "n_samples": len(manifest_rows)}, handle, indent=2)
        handle.write("\n")
    print(f"WROTE {len(manifest_rows)} matrix-cohort samples", flush=True)


if __name__ == "__main__":
    main()
