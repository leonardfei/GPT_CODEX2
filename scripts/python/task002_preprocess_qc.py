#!/usr/bin/env python3
"""Task 002: sample-aware, neutrophil-friendly scRNA-seq QC.

This script runs on the compute server. It never changes raw inputs. It
extracts a server-only Matrix Market staging representation, calculates
cell-level QC and a conservative coherent neutrophil candidate screen, chooses
sample-aware thresholds from the observed candidate distributions, writes the
required audit tables, and creates one multi-page review PDF.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import os
import re
import shutil
import tarfile
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator, Sequence

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_matrix, csr_matrix, issparse, vstack


APPROVED_DATASETS = {
    "GSE282701",
    "GSE242889",
    "GSE326201",
    "GSE149614",
    "GSE299340",
}

CORE_MARKERS = ["FCGR3B", "CSF3R", "CXCR2", "FPR1"]
SUPPORT_MARKERS = [
    "S100A8",
    "S100A9",
    "CTSG",
    "ELANE",
    "MPO",
    "LYZ",
    "MCEMP1",
    "FCAR",
    "FFAR2",
    "NAMPT",
    "CXCR4",
]


@dataclass
class SampleSpec:
    dataset: str
    sample_id: str
    sample_title: str
    patient_id: str
    tissue: str
    paired_status: str
    paired_id: str
    etiology: str
    mvi: str
    platform: str
    selection_strategy: str
    matrix_type: str
    data_representation: str
    phase1_role: str
    notes: str


@dataclass
class MatrixData:
    matrix: csr_matrix
    gene_ids: list[str]
    gene_names: list[str]
    feature_types: list[str]
    barcodes: list[str]
    source_member: str
    structure_note: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        required=True,
        type=Path,
        help="Server project root, e.g. /data/lf_data/...",
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    return parser.parse_args()


def read_manifest(path: Path) -> list[SampleSpec]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected: list[SampleSpec] = []
    for row in rows:
        dataset = row["dataset"]
        if dataset not in APPROVED_DATASETS:
            continue
        if dataset == "GSE149614" and row["tissue"] not in {"Tumor", "Adjacent"}:
            continue
        selected.append(
            SampleSpec(
                dataset=dataset,
                sample_id=row["sample_id"],
                sample_title=row.get("sample_title", ""),
                patient_id=row.get("patient_id", ""),
                tissue=row.get("tissue", ""),
                paired_status=row.get("paired_status", ""),
                paired_id=row.get("paired_id", ""),
                etiology=row.get("etiology", ""),
                mvi=row.get("MVI", ""),
                platform=row.get("platform", ""),
                selection_strategy=row.get("selection_strategy", ""),
                matrix_type=row.get("matrix_type", ""),
                data_representation=row.get("data_representation", ""),
                phase1_role=row.get("phase1_role", ""),
                notes=row.get("notes", ""),
            )
        )
    if not selected:
        raise RuntimeError("No approved samples found in manifest")
    return selected


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def unique_names(names: Sequence[str]) -> list[str]:
    seen: dict[str, int] = {}
    result: list[str] = []
    for name in names:
        base = name.strip() or "feature"
        count = seen.get(base, 0)
        result.append(base if count == 0 else f"{base}.{count}")
        seen[base] = count + 1
    return result


def read_text_lines(handle: BinaryIO) -> Iterator[str]:
    for raw in handle:
        if isinstance(raw, bytes):
            yield raw.decode("utf-8", errors="replace").rstrip("\r\n")
        else:
            yield str(raw).rstrip("\r\n")


def parse_features(handle: BinaryIO) -> tuple[list[str], list[str], list[str]]:
    ids: list[str] = []
    names: list[str] = []
    types: list[str] = []
    for line in read_text_lines(handle):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(f"Feature line has fewer than two columns: {line!r}")
        ids.append(fields[0])
        names.append(fields[1] or fields[0])
        types.append(fields[2] if len(fields) >= 3 else "Gene Expression")
    if not ids:
        raise ValueError("No features found")
    return ids, names, types


def parse_barcodes(handle: BinaryIO) -> list[str]:
    barcodes = [line for line in read_text_lines(handle) if line]
    if not barcodes:
        raise ValueError("No barcodes found")
    return barcodes


def read_mtx_stream(handle: BinaryIO, compressed: bool) -> csr_matrix:
    stream: BinaryIO | gzip.GzipFile = handle
    if compressed:
        stream = gzip.GzipFile(fileobj=handle, mode="rb")
    matrix = mmread(stream)
    if not issparse(matrix):
        matrix = csr_matrix(np.asarray(matrix))
    return matrix.tocsr().astype(np.int32)


def find_member(names: Iterable[str], patterns: Sequence[str]) -> str | None:
    names = list(names)
    for pattern in patterns:
        exact = [name for name in names if name == pattern]
        if exact:
            return exact[0]
    for pattern in patterns:
        matches = [name for name in names if name.endswith(pattern)]
        if matches:
            return matches[0]
    return None


def load_10x_matrix_from_tar(
    tar: tarfile.TarFile,
    sample_id: str,
    required_suffixes: tuple[str, str, str],
) -> MatrixData:
    names = tar.getnames()
    sample_members = [name for name in names if Path(name).name.startswith(sample_id + "_")]
    if not sample_members:
        raise FileNotFoundError(f"No archive members found for {sample_id}")
    matrix_member = next(
        (name for name in sample_members if name.endswith(required_suffixes[0])), None
    )
    feature_member = next(
        (name for name in sample_members if name.endswith(required_suffixes[1])), None
    )
    barcode_member = next(
        (name for name in sample_members if name.endswith(required_suffixes[2])), None
    )
    if not all([matrix_member, feature_member, barcode_member]):
        raise FileNotFoundError(
            f"Incomplete 10x members for {sample_id}: "
            f"matrix={matrix_member}, features={feature_member}, barcodes={barcode_member}"
        )
    matrix_file = tar.extractfile(matrix_member)
    feature_file = tar.extractfile(feature_member)
    barcode_file = tar.extractfile(barcode_member)
    if not matrix_file or not feature_file or not barcode_file:
        raise OSError(f"Unable to extract all members for {sample_id}")
    with matrix_file, feature_file, barcode_file:
        matrix = read_mtx_stream(matrix_file, matrix_member.endswith(".gz"))
        if feature_member.endswith(".gz"):
            with gzip.GzipFile(fileobj=feature_file, mode="rb") as features_stream:
                gene_ids, gene_names, feature_types = parse_features(features_stream)
        else:
            gene_ids, gene_names, feature_types = parse_features(feature_file)
        if barcode_member.endswith(".gz"):
            with gzip.GzipFile(fileobj=barcode_file, mode="rb") as barcodes_stream:
                barcodes = parse_barcodes(barcodes_stream)
        else:
            barcodes = parse_barcodes(barcode_file)
    if matrix.shape != (len(gene_names), len(barcodes)):
        raise ValueError(
            f"Dimension mismatch for {sample_id}: matrix={matrix.shape}, "
            f"features={len(gene_names)}, barcodes={len(barcodes)}"
        )
    return MatrixData(
        matrix=matrix,
        gene_ids=gene_ids,
        gene_names=gene_names,
        feature_types=feature_types,
        barcodes=barcodes,
        source_member=f"{matrix_member};{feature_member};{barcode_member}",
        structure_note="feature-by-barcode 10x Matrix Market; identifiers recovered",
    )


def load_gse242889_sample(outer: tarfile.TarFile, sample_id: str) -> MatrixData:
    outer_name = next(
        (name for name in outer.getnames() if Path(name).name.startswith(sample_id + "_")), None
    )
    if not outer_name:
        raise FileNotFoundError(f"No GSE242889 nested archive for {sample_id}")
    outer_file = outer.extractfile(outer_name)
    if not outer_file:
        raise OSError(f"Unable to extract nested archive for {sample_id}")
    with outer_file:
        archive_bytes = outer_file.read()
    with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:*") as inner:
        names = inner.getnames()
        matrix_member = find_member(names, ["matrix.mtx", "matrix.mtx.gz"])
        barcode_member = find_member(names, ["barcodes.tsv", "barcodes.tsv.gz"])
        feature_member = find_member(names, ["genes.tsv", "features.tsv", "genes.tsv.gz", "features.tsv.gz"])
        if not all([matrix_member, barcode_member, feature_member]):
            raise FileNotFoundError(
                f"GSE242889 identifiers not recoverable for {sample_id}: "
                f"matrix={matrix_member}, barcodes={barcode_member}, features={feature_member}"
            )
        matrix_file = inner.extractfile(matrix_member)
        barcode_file = inner.extractfile(barcode_member)
        feature_file = inner.extractfile(feature_member)
        if not matrix_file or not barcode_file or not feature_file:
            raise OSError(f"Unable to extract GSE242889 members for {sample_id}")
        with matrix_file, barcode_file, feature_file:
            matrix = read_mtx_stream(matrix_file, matrix_member.endswith(".gz"))
            if feature_member.endswith(".gz"):
                with gzip.GzipFile(fileobj=feature_file, mode="rb") as features_stream:
                    gene_ids, gene_names, feature_types = parse_features(features_stream)
            else:
                gene_ids, gene_names, feature_types = parse_features(feature_file)
            if barcode_member.endswith(".gz"):
                with gzip.GzipFile(fileobj=barcode_file, mode="rb") as barcodes_stream:
                    barcodes = parse_barcodes(barcodes_stream)
            else:
                barcodes = parse_barcodes(barcode_file)
    if matrix.shape != (len(gene_names), len(barcodes)):
        raise ValueError(
            f"Dimension mismatch for GSE242889 {sample_id}: matrix={matrix.shape}, "
            f"features={len(gene_names)}, barcodes={len(barcodes)}"
        )
    return MatrixData(
        matrix=matrix,
        gene_ids=gene_ids,
        gene_names=gene_names,
        feature_types=feature_types,
        barcodes=barcodes,
        source_member=f"{outer_name}!{matrix_member};{feature_member};{barcode_member}",
        structure_note="nested POSIX tar; matrix.mtx/genes.tsv/barcodes.tsv; identifiers recovered",
    )


def h5_values(value: object) -> list[str]:
    array = np.asarray(value)
    return [decode(item) for item in array.tolist()]


def load_gse326201_h5(outer: tarfile.TarFile, sample_id: str, temp_dir: Path) -> MatrixData:
    member = next(
        (name for name in outer.getnames() if Path(name).name.startswith(sample_id + "_")), None
    )
    if not member:
        raise FileNotFoundError(f"No GSE326201 H5 member for {sample_id}")
    extracted = temp_dir / f"{sample_id}.h5"
    with outer.extractfile(member) as source, extracted.open("wb") as target:
        if source is None:
            raise OSError(f"Unable to extract {member}")
        shutil.copyfileobj(source, target, length=1024 * 1024)
    with h5py.File(extracted, "r") as h5:
        root = h5["matrix"]
        data = np.asarray(root["data"][:])
        indices = np.asarray(root["indices"][:], dtype=np.int64)
        indptr = np.asarray(root["indptr"][:], dtype=np.int64)
        shape = tuple(int(x) for x in root["shape"][:])
        matrix = csc_matrix((data, indices, indptr), shape=shape).tocsr().astype(np.int32)
        barcodes = h5_values(root["barcodes"][:])
        features = root["features"]
        if "name" in features:
            gene_names = h5_values(features["name"][:])
        elif "gene_names" in features:
            gene_names = h5_values(features["gene_names"][:])
        else:
            raise KeyError(f"No feature names in {member}")
        if "id" in features:
            gene_ids = h5_values(features["id"][:])
        else:
            gene_ids = list(gene_names)
        if "feature_type" in features:
            feature_types = h5_values(features["feature_type"][:])
        else:
            feature_types = ["Gene Expression"] * len(gene_names)
    extracted.unlink(missing_ok=True)
    if matrix.shape != (len(gene_names), len(barcodes)):
        raise ValueError(
            f"Dimension mismatch for GSE326201 {sample_id}: matrix={matrix.shape}, "
            f"features={len(gene_names)}, barcodes={len(barcodes)}"
        )
    return MatrixData(
        matrix=matrix,
        gene_ids=gene_ids,
        gene_names=gene_names,
        feature_types=feature_types,
        barcodes=barcodes,
        source_member=member,
        structure_note="Cell Ranger filtered feature-barcode HDF5; identifiers recovered",
    )


def load_gse149614_inputs(project_root: Path) -> tuple[pd.DataFrame, list[str], dict[str, dict[str, object]]]:
    raw = project_root / "raw_data" / "GSE149614"
    metadata_path = raw / "GSE149614_HCC.metadata.updated.txt.gz"
    counts_path = raw / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str, compression="gzip")
    # GEO metadata use both "Adjacent" and "Normal" for non-tumour tissue;
    # the phase-1 manifest maps the latter to the standardized Adjacent label.
    metadata = metadata[metadata["site"].isin(["Tumor", "Adjacent", "Normal", "NTL"])].copy()
    with gzip.open(counts_path, "rt", encoding="utf-8", errors="replace") as handle:
        header_line = handle.readline().rstrip("\r\n")
    count_cells = header_line.split("\t")
    if len(count_cells) != len(set(count_cells)):
        raise ValueError("GSE149614 count header contains duplicate cell identifiers")
    metadata_cells = metadata["Cell"].astype(str).tolist()
    missing = sorted(set(metadata_cells) - set(count_cells))
    if missing:
        raise ValueError(f"GSE149614 metadata cells absent from counts: {missing[:3]}")
    cell_index = {cell: index for index, cell in enumerate(count_cells)}
    sample_info: dict[str, dict[str, object]] = {}
    for sample_id, group in metadata.groupby("sample", sort=False):
        sample_cells = group["Cell"].astype(str).tolist()
        sample_info[str(sample_id)] = {
            "barcodes": sample_cells,
            "indices": [cell_index[cell] for cell in sample_cells],
        }
    return metadata, count_cells, sample_info


def iterate_gse149614_chunks(
    project_root: Path,
    count_cells: list[str],
    sample_info: dict[str, dict[str, object]],
    chunksize: int = 256,
) -> Iterator[tuple[list[str], dict[str, csr_matrix]]]:
    counts_path = project_root / "raw_data/GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    dtype_map = {cell: np.int32 for cell in count_cells}
    reader = pd.read_csv(
        counts_path,
        sep="\t",
        compression="gzip",
        header=None,
        skiprows=1,
        names=["gene"] + count_cells,
        index_col=0,
        dtype=dtype_map,
        chunksize=chunksize,
    )
    for chunk in reader:
        gene_names = [str(value) for value in chunk.index.tolist()]
        parts: dict[str, csr_matrix] = {}
        for sample_id, info in sample_info.items():
            indices = info["indices"]
            dense = chunk.iloc[:, indices].to_numpy(dtype=np.int32, copy=True)
            parts[sample_id] = csr_matrix(dense)
            del dense
        yield gene_names, parts


def matrix_to_stage(data: MatrixData, stage_dir: Path) -> tuple[Path, Path, Path]:
    stage_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = stage_dir / "matrix.mtx"
    features_path = stage_dir / "features.tsv"
    barcodes_path = stage_dir / "barcodes.tsv"
    mmwrite(matrix_path, data.matrix.tocoo())
    with features_path.open("w", encoding="utf-8") as handle:
        for gene_id, gene_name, feature_type in zip(
            data.gene_ids, unique_names(data.gene_names), data.feature_types
        ):
            handle.write(f"{gene_id}\t{gene_name}\t{feature_type}\n")
    with barcodes_path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(data.barcodes) + "\n")
    return matrix_path, features_path, barcodes_path


def robust_iqr(values: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return (float(median), float(q1), float(q3))


def mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    med = np.median(values)
    return float(np.median(np.abs(values - med)))


def choose_thresholds(
    n_features: np.ndarray,
    n_counts: np.ndarray,
    percent_mt: np.ndarray,
    candidate: np.ndarray,
) -> dict[str, object]:
    candidate_count = int(candidate.sum())
    all_feature_floor = max(50, int(math.floor(np.quantile(n_features, 0.01))))
    all_count_floor = max(100, int(math.floor(np.quantile(n_counts, 0.01))))
    all_mt_ceiling = float(max(20.0, np.quantile(percent_mt, 0.98)))
    if candidate_count >= 10:
        candidate_features = n_features[candidate]
        candidate_counts = n_counts[candidate]
        candidate_mt = percent_mt[candidate]
        min_features = max(50, int(math.floor(np.quantile(candidate_features, 0.05))))
        min_counts = max(100, int(math.floor(np.quantile(candidate_counts, 0.05))))
        max_mt = float(max(20.0, math.ceil(np.quantile(candidate_mt, 0.99))))
        basis = "candidate_5th_percentile_for_lower_bounds_and_candidate_99th_percentile_mt"
    else:
        global_median_features = float(np.median(n_features))
        global_median_counts = float(np.median(n_counts))
        min_features = max(
            50,
            int(math.floor(global_median_features - 3.0 * max(mad(n_features), 1.0))),
            all_feature_floor,
        )
        min_counts = max(
            100,
            int(math.floor(global_median_counts - 3.0 * max(mad(n_counts), 1.0))),
            all_count_floor,
        )
        max_mt = all_mt_ceiling
        basis = "fallback_global_robust_floor_and_global_98th_percentile_mt_no_candidate_support"
    return {
        "candidate_count": candidate_count,
        "min_nFeature_RNA": int(min_features),
        "min_nCount_RNA": int(min_counts),
        "max_percent_mt": float(max_mt),
        "threshold_basis": basis,
        "all_feature_q01": float(np.quantile(n_features, 0.01)),
        "all_count_q01": float(np.quantile(n_counts, 0.01)),
        "all_percent_mt_q98": float(np.quantile(percent_mt, 0.98)),
        "candidate_feature_q05": float(np.quantile(n_features[candidate], 0.05)) if candidate_count else float("nan"),
        "candidate_count_q05": float(np.quantile(n_counts[candidate], 0.05)) if candidate_count else float("nan"),
        "candidate_percent_mt_q99": float(np.quantile(percent_mt[candidate], 0.99)) if candidate_count else float("nan"),
    }


def gene_mask_map(gene_names: Sequence[str]) -> dict[str, np.ndarray]:
    normalized: dict[str, list[int]] = defaultdict(list)
    for index, gene in enumerate(gene_names):
        key = re.sub(r"\..*$", "", str(gene).upper())
        normalized[key].append(index)
    return {key: np.asarray(indices, dtype=np.int64) for key, indices in normalized.items()}


def marker_hits(matrix: csr_matrix, gene_names: Sequence[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    gene_map = gene_mask_map(gene_names)
    core_hits = np.zeros(matrix.shape[1], dtype=np.int16)
    support_hits = np.zeros(matrix.shape[1], dtype=np.int16)
    missing: list[str] = []
    for marker in CORE_MARKERS:
        rows = gene_map.get(marker)
        if rows is None:
            missing.append(marker)
            continue
        core_hits += np.asarray(matrix[rows, :].sum(axis=0)).ravel() > 0
    for marker in SUPPORT_MARKERS:
        rows = gene_map.get(marker)
        if rows is None:
            missing.append(marker)
            continue
        support_hits += np.asarray(matrix[rows, :].sum(axis=0)).ravel() > 0
    # Require coherent evidence from at least three markers, with a core signal
    # or the paired S100A8/S100A9 context. This deliberately avoids a one-marker rescue.
    s100a8 = np.zeros(matrix.shape[1], dtype=bool)
    s100a9 = np.zeros(matrix.shape[1], dtype=bool)
    for marker, target in [("S100A8", s100a8), ("S100A9", s100a9)]:
        rows = gene_map.get(marker)
        if rows is not None:
            target[:] = np.asarray(matrix[rows, :].sum(axis=0)).ravel() > 0
    candidate = ((core_hits >= 1) & (support_hits >= 2)) | (
        (support_hits >= 3) & s100a8 & s100a9
    )
    return core_hits, support_hits, sorted(set(missing))


def qc_reason(row: pd.Series) -> str:
    reasons: list[str] = []
    if row["nFeature_RNA"] < row["min_nFeature_RNA"]:
        reasons.append("low_nFeature_RNA")
    if row["nCount_RNA"] < row["min_nCount_RNA"]:
        reasons.append("low_nCount_RNA")
    if row["percent.mt"] > row["max_percent_mt"]:
        reasons.append("high_percent.mt")
    return "pass" if not reasons else ";".join(reasons)


def make_review_page(
    pdf: PdfPages,
    sample: SampleSpec,
    qc: pd.DataFrame,
    thresholds: dict[str, object],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    candidate = qc["candidate_neutrophil"].astype(bool).to_numpy()
    passed = qc["passes_qc"].astype(bool).to_numpy()
    colors = np.where(candidate, "#b2182b", "#969696")
    axes[0, 0].scatter(
        np.log10(qc["nFeature_RNA"].to_numpy() + 1),
        np.log10(qc["nCount_RNA"].to_numpy() + 1),
        c=colors,
        s=4,
        alpha=0.35,
        linewidths=0,
    )
    axes[0, 0].set_xlabel("log10(nFeature_RNA + 1)")
    axes[0, 0].set_ylabel("log10(nCount_RNA + 1)")
    axes[0, 0].set_title("All cells; red = coherent neutrophil candidate")
    axes[0, 0].axvline(np.log10(float(thresholds["min_nFeature_RNA"]) + 1), color="black", ls="--", lw=0.8)
    axes[0, 0].axhline(np.log10(float(thresholds["min_nCount_RNA"]) + 1), color="black", ls="--", lw=0.8)

    for mask, label, color in [
        (np.ones(len(qc), dtype=bool), "all", "#777777"),
        (candidate, "candidate", "#b2182b"),
    ]:
        if mask.any():
            axes[0, 1].hist(
                np.log10(qc.loc[mask, "nFeature_RNA"] + 1),
                bins=40,
                histtype="step",
                lw=1.3,
                label=label,
                color=color,
            )
    axes[0, 1].axvline(np.log10(float(thresholds["min_nFeature_RNA"]) + 1), color="black", ls="--", lw=0.8)
    axes[0, 1].set_xlabel("log10(nFeature_RNA + 1)")
    axes[0, 1].set_ylabel("cells")
    axes[0, 1].set_title("Detected genes")
    axes[0, 1].legend(frameon=False)

    metric_data = []
    metric_labels = []
    for metric in ["nFeature_RNA", "nCount_RNA", "percent.mt"]:
        for label, mask in [("candidate", candidate), ("other", ~candidate)]:
            values = qc.loc[mask, metric].to_numpy(dtype=float)
            if values.size:
                metric_data.append(values)
                metric_labels.append(f"{metric}\n{label}")
    if metric_data:
        axes[1, 0].boxplot(metric_data, labels=metric_labels, showfliers=False)
    axes[1, 0].set_title("Candidate vs other distributions")
    axes[1, 0].tick_params(axis="x", labelrotation=70, labelsize=7)
    axes[1, 0].set_ylabel("value; percent.mt is %")

    labels = ["all cells", "candidate pre-QC", "candidate post-QC", "doublet proxy"]
    values = [
        len(qc),
        int(candidate.sum()),
        int((candidate & passed).sum()),
        int(qc["doublet_risk_proxy"].sum()),
    ]
    axes[1, 1].bar(labels, values, color=["#777777", "#b2182b", "#2166ac", "#f4a582"])
    axes[1, 1].set_title("Retention and review flags")
    axes[1, 1].tick_params(axis="x", labelrotation=35)
    axes[1, 1].set_ylabel("cells")
    for index, value in enumerate(values):
        axes[1, 1].text(index, value, str(value), ha="center", va="bottom", fontsize=8)

    fig.suptitle(
        f"Task 002 QC | {sample.dataset} | {sample.sample_id} | {sample.tissue} | {sample.patient_id}\n"
        f"min genes={thresholds['min_nFeature_RNA']}; min UMIs={thresholds['min_nCount_RNA']}; "
        f"max mt={thresholds['max_percent_mt']:.2f}%",
        fontsize=11,
    )
    pdf.savefig(fig)
    plt.close(fig)


def summarize_metric(values: np.ndarray, prefix: str) -> dict[str, float]:
    median, q1, q3 = robust_iqr(values)
    return {f"{prefix}_median": median, f"{prefix}_q1": q1, f"{prefix}_q3": q3, f"{prefix}_iqr": q3 - q1}


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    out_root = args.out_root.resolve()
    raw_root = project_root / "raw_data"
    stage_root = out_root / "staged_matrices"
    cell_qc_root = out_root / "cell_qc"
    results_root = project_root / "results"
    figures_root = project_root / "figures"
    logs_root = project_root / "logs"
    for path in [stage_root, cell_qc_root, results_root, figures_root, logs_root]:
        path.mkdir(parents=True, exist_ok=True)
    samples = read_manifest(args.manifest)
    by_dataset: dict[str, list[SampleSpec]] = defaultdict(list)
    for sample in samples:
        by_dataset[sample.dataset].append(sample)

    audit_rows: list[dict[str, object]] = []
    threshold_rows: list[dict[str, object]] = []
    structure_rows: list[dict[str, object]] = []
    stage_rows: list[dict[str, object]] = []
    all_qc_for_pdf: list[tuple[SampleSpec, pd.DataFrame, dict[str, object]]] = []

    gse149_count_cells: list[str] | None = None
    gse149_sample_info: dict[str, dict[str, object]] | None = None
    gse149_parts: dict[str, list[csr_matrix]] = {}
    gse149_gene_names: list[str] = []
    if "GSE149614" in by_dataset:
        _, gse149_count_cells, gse149_sample_info = load_gse149614_inputs(project_root)
        gse149_parts = {sample_id: [] for sample_id in gse149_sample_info}
        for gene_chunk, matrix_parts in iterate_gse149614_chunks(
            project_root, gse149_count_cells, gse149_sample_info
        ):
            gse149_gene_names.extend(gene_chunk)
            for sample_id, matrix_part in matrix_parts.items():
                gse149_parts[sample_id].append(matrix_part)

    with tempfile.TemporaryDirectory(prefix="task002_h5_") as temp_name:
        temp_dir = Path(temp_name)
        for dataset, dataset_samples in by_dataset.items():
            if dataset == "GSE149614":
                assert gse149_count_cells is not None and gse149_sample_info is not None
                for sample in dataset_samples:
                    sample_key = sample.sample_title if sample.sample_title in gse149_sample_info else sample.sample_id
                    if sample_key not in gse149_sample_info:
                        raise KeyError(
                            f"GSE149614 sample mapping missing for GEO {sample.sample_id}; "
                            f"tried sample title {sample.sample_title!r}"
                        )
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
                    structure_format = "combined_count_text_plus_metadata"
                    stage_source = data.source_member
                    process_sample(
                        sample,
                        data,
                        structure_format,
                        stage_source,
                        stage_root,
                        cell_qc_root,
                        audit_rows,
                        threshold_rows,
                        structure_rows,
                        stage_rows,
                        all_qc_for_pdf,
                    )
                    del matrix, data
                continue

            archive_path = raw_root / dataset / f"{dataset}_RAW.tar"
            with tarfile.open(archive_path, mode="r:") as archive:
                for sample in dataset_samples:
                    if dataset == "GSE282701":
                        data = load_10x_matrix_from_tar(
                            archive,
                            sample.sample_id,
                            ("_matrix.mtx.gz", "_features.tsv.gz", "_barcodes.tsv.gz"),
                        )
                        structure_format = "10x_matrix_market_gz"
                    elif dataset == "GSE299340":
                        data = load_10x_matrix_from_tar(
                            archive,
                            sample.sample_id,
                            ("_matrix.mtx.gz", "_features.tsv.gz", "_barcodes.tsv.gz"),
                        )
                        structure_format = "10x_matrix_market_gz"
                    elif dataset == "GSE242889":
                        data = load_gse242889_sample(archive, sample.sample_id)
                        structure_format = "nested_tar_matrix_market"
                    elif dataset == "GSE326201":
                        data = load_gse326201_h5(archive, sample.sample_id, temp_dir)
                        structure_format = "cellranger_filtered_h5"
                    else:
                        raise ValueError(f"Unexpected dataset: {dataset}")
                    process_sample(
                        sample,
                        data,
                        structure_format,
                        data.source_member,
                        stage_root,
                        cell_qc_root,
                        audit_rows,
                        threshold_rows,
                        structure_rows,
                        stage_rows,
                        all_qc_for_pdf,
                    )

    audit_df = pd.DataFrame(audit_rows).sort_values(["dataset", "sample_id"])
    threshold_df = pd.DataFrame(threshold_rows).sort_values(["dataset", "sample_id"])
    structure_df = pd.DataFrame(structure_rows).sort_values(["dataset", "sample_id"])
    stage_df = pd.DataFrame(stage_rows).sort_values(["dataset", "sample_id"])
    audit_df.to_csv(results_root / "task002_neutrophil_retention_audit.csv", index=False)
    threshold_df.to_csv(results_root / "task002_qc_thresholds_by_sample.csv", index=False)
    structure_df.to_csv(results_root / "task002_input_structure_audit.csv", index=False)
    stage_df.to_csv(out_root / "task002_stage_manifest.csv", index=False)
    with (figures_root / "task002_qc_review.pdf").open("wb") as handle:
        pdf = PdfPages(handle)
        for sample, qc, thresholds in all_qc_for_pdf:
            make_review_page(pdf, sample, qc, thresholds)
        pdf.close()
    run_meta = {
        "task": "task_002",
        "script": "task002_preprocess_qc.py",
        "approved_datasets": sorted(APPROVED_DATASETS),
        "excluded_datasets": ["GSE202642", "GSE290298"],
        "n_samples": len(samples),
        "candidate_definition": "(core_hits>=1 and support_hits>=2) or (support_hits>=3 and S100A8 and S100A9 expressed)",
        "qc_filter": "nFeature_RNA >= sample_min_nFeature_RNA and nCount_RNA >= sample_min_nCount_RNA and percent.mt <= sample_max_percent_mt",
        "doublet_evaluation": "high nCount_RNA and nFeature_RNA 99.5th-percentile proxy retained as review flag; no hard removal",
        "raw_counts_preserved": True,
    }
    (results_root / "task002_run_metadata.json").write_text(json.dumps(run_meta, indent=2) + "\n", encoding="utf-8")
    print(f"Processed {len(samples)} samples across {len(by_dataset)} datasets")
    print(f"Audit: {results_root / 'task002_neutrophil_retention_audit.csv'}")
    print(f"Thresholds: {results_root / 'task002_qc_thresholds_by_sample.csv'}")


def process_sample(
    sample: SampleSpec,
    data: MatrixData,
    structure_format: str,
    stage_source: str,
    stage_root: Path,
    cell_qc_root: Path,
    audit_rows: list[dict[str, object]],
    threshold_rows: list[dict[str, object]],
    structure_rows: list[dict[str, object]],
    stage_rows: list[dict[str, object]],
    all_qc_for_pdf: list[tuple[SampleSpec, pd.DataFrame, dict[str, object]]],
) -> None:
    matrix = data.matrix.tocsr()
    n_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    n_features = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(float)
    gene_map = gene_mask_map(data.gene_names)
    mt_rows = np.concatenate(
        [gene_map[key] for key in gene_map if key.startswith("MT-") or key.startswith("MT.")]
    ) if any(key.startswith("MT-") or key.startswith("MT.") for key in gene_map) else np.array([], dtype=np.int64)
    mt_counts = np.asarray(matrix[mt_rows, :].sum(axis=0)).ravel() if mt_rows.size else np.zeros(matrix.shape[1])
    percent_mt = np.divide(mt_counts * 100.0, n_counts, out=np.zeros_like(n_counts), where=n_counts > 0)
    core_hits, support_hits, missing_markers = marker_hits(matrix, data.gene_names)
    candidate = ((core_hits >= 1) & (support_hits >= 2)) | (
        (support_hits >= 3)
        & (np.asarray(matrix[gene_map["S100A8"], :].sum(axis=0)).ravel() > 0 if "S100A8" in gene_map else False)
        & (np.asarray(matrix[gene_map["S100A9"], :].sum(axis=0)).ravel() > 0 if "S100A9" in gene_map else False)
    )
    thresholds = choose_thresholds(n_features, n_counts, percent_mt, candidate)
    passes = (
        (n_features >= thresholds["min_nFeature_RNA"])
        & (n_counts >= thresholds["min_nCount_RNA"])
        & (percent_mt <= thresholds["max_percent_mt"])
    )
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
            "candidate_neutrophil": candidate,
            "passes_qc": passes,
            "doublet_risk_proxy": doublet_proxy,
            "min_nFeature_RNA": thresholds["min_nFeature_RNA"],
            "min_nCount_RNA": thresholds["min_nCount_RNA"],
            "max_percent_mt": thresholds["max_percent_mt"],
        }
    )
    qc["removal_reason"] = qc.apply(qc_reason, axis=1)

    # Stage only cells that pass the reviewed sample-specific QC. The source
    # count matrix remains untouched; the full pre-QC cell metrics are retained
    # in cell_qc so the retention audit remains reproducible.
    kept_indices = np.flatnonzero(passes)
    staged_data = MatrixData(
        matrix=matrix[:, kept_indices].tocsr(),
        gene_ids=list(data.gene_ids),
        gene_names=list(data.gene_names),
        feature_types=list(data.feature_types),
        barcodes=[data.barcodes[index] for index in kept_indices],
        source_member=data.source_member,
        structure_note=data.structure_note + "; staged matrix contains passes_qc cells",
    )
    stage_dir = stage_root / sample.dataset / sample.sample_id
    matrix_path, features_path, barcodes_path = matrix_to_stage(staged_data, stage_dir)
    cell_qc_path = cell_qc_root / sample.dataset / f"{sample.sample_id}.tsv.gz"
    cell_qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(cell_qc_path, sep="\t", index=False, compression="gzip")
    stage_rows.append(
        {
            **sample.__dict__,
            "matrix_path": str(matrix_path),
            "features_path": str(features_path),
            "barcodes_path": str(barcodes_path),
            "cell_qc_path": str(cell_qc_path),
            "n_features": matrix.shape[0],
            "n_cells_before_qc": matrix.shape[1],
            "n_cells_staged_after_qc": staged_data.matrix.shape[1],
            "source_member": stage_source,
        }
    )
    structure_rows.append(
        {
            **sample.__dict__,
            "structure_format": structure_format,
            "source_member": stage_source,
            "n_features": matrix.shape[0],
            "n_cells": matrix.shape[1],
            "feature_ids_recovered": bool(len(data.gene_ids) == matrix.shape[0] and all(data.gene_ids)),
            "feature_names_recovered": bool(len(data.gene_names) == matrix.shape[0] and all(data.gene_names)),
            "barcodes_recovered": bool(len(data.barcodes) == matrix.shape[1] and all(data.barcodes)),
            "orientation": "features_by_barcodes",
            "structure_check": "PASS",
            "structure_note": data.structure_note,
            "missing_marker_genes": ";".join(missing_markers),
        }
    )
    threshold_rows.append(
        {
            **sample.__dict__,
            **thresholds,
            "n_cells_before_qc": int(len(qc)),
            "n_cells_after_qc": int(passes.sum()),
            "n_candidate_neutrophils_before_qc": int(candidate.sum()),
            "n_candidate_neutrophils_after_qc": int((candidate & passes).sum()),
            "n_doublet_proxy_review": int(doublet_proxy.sum()),
            "missing_marker_genes": ";".join(missing_markers),
        }
    )
    candidate_metrics: dict[str, float] = {}
    noncandidate_metrics: dict[str, float] = {}
    for metric in ["nFeature_RNA", "nCount_RNA", "percent.mt"]:
        candidate_metrics.update(summarize_metric(qc.loc[candidate, metric].to_numpy(), f"candidate_{metric}"))
        noncandidate_metrics.update(summarize_metric(qc.loc[~candidate, metric].to_numpy(), f"noncandidate_{metric}"))
    candidate_before = int(candidate.sum())
    candidate_after = int((candidate & passes).sum())
    removal_counts = qc.loc[candidate & ~passes, "removal_reason"].value_counts().to_dict()
    audit_rows.append(
        {
            **sample.__dict__,
            "n_cells_before_qc": int(len(qc)),
            "n_cells_after_qc": int(passes.sum()),
            "n_candidate_neutrophils_before_qc": candidate_before,
            "n_candidate_neutrophils_after_qc": candidate_after,
            "neutrophil_retention_fraction": (candidate_after / candidate_before if candidate_before else float("nan")),
            "candidate_removed_low_nFeature_RNA": int((candidate & (n_features < thresholds["min_nFeature_RNA"])).sum()),
            "candidate_removed_low_nCount_RNA": int((candidate & (n_counts < thresholds["min_nCount_RNA"])).sum()),
            "candidate_removed_high_percent_mt": int((candidate & (percent_mt > thresholds["max_percent_mt"])).sum()),
            "candidate_major_removal_reasons": ";".join(f"{key}:{value}" for key, value in sorted(removal_counts.items())) or "none",
            "flag_suspiciously_low_retention": bool(candidate_before >= 10 and candidate_after / candidate_before < 0.8),
            "n_doublet_proxy_review": int(doublet_proxy.sum()),
            **candidate_metrics,
            **noncandidate_metrics,
            "missing_marker_genes": ";".join(missing_markers),
            "structure_note": data.structure_note,
        }
    )
    all_qc_for_pdf.append((sample, qc, thresholds))


if __name__ == "__main__":
    main()
