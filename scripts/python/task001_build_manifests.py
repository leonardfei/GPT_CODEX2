#!/usr/bin/env python3
"""Build Task 001 sample and download manifests from GEO SOFT metadata.

The script intentionally does metadata normalization only. It does not download
matrices, filter cells, merge objects, or infer missing patient relationships.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


ACCESSIONS = [
    "GSE282701",
    "GSE242889",
    "GSE326201",
    "GSE149614",
    "GSE299340",
    "GSE290298",
    "GSE202642",
]


def ftp_series(accession: str, filename: str) -> str:
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:6]}nnn/"
        f"{accession}/suppl/{filename}"
    )


DATASET = {
    "GSE282701": {
        "matrix_type": "filtered_counts",
        "data_representation": "12 per-sample 10x-style Matrix Market matrices inside GSE282701_RAW.tar",
        "platform": "SeekOne Digital Droplet 3-prime / Illumina NovaSeq 6000; Seeksoultools 1.2.2; GRCh38-2020-A",
        "selection": "Unselected (no enrichment stated)",
        "matrix_url": ftp_series("GSE282701", "GSE282701_RAW.tar"),
        "matrix_size_bytes": "856524800",
        "fastq_status": "No SRA hit in accession-term search; GEO matrix archive available",
        "abundance": "YES, paired unselected tumour/adjacent design; MVI is not explicitly sample-annotated",
        "ambiguity": "Series states no vascular invasion/metastatic features, but a dedicated MVI field is not supplied",
        "source_cells": "121848 cells (Series overall design; source-reported)",
        "raw_available": "NO separate unfiltered raw matrix observed",
        "filtered_available": "YES",
        "normalized_available": "NO",
        "role": "phase1_preferred",
    },
    "GSE242889": {
        "matrix_type": "raw_counts",
        "data_representation": "10 per-sample tar.gz archives; representative archive contains Matrix Market matrix.mtx",
        "platform": "Singleron GEXSCOPE / Illumina NovaSeq 6000; CeleScope 1.9.0; GRCh38 Ensembl 92",
        "selection": "Unselected (RBC lysis only; no cell enrichment stated)",
        "matrix_url": ftp_series("GSE242889", "GSE242889_RAW.tar"),
        "matrix_size_bytes": "1376245760",
        "fastq_status": "No SRA hit in accession-term search; GEO matrix archives available",
        "abundance": "YES, five paired tumour/adjacent pairs; analyse at patient level and retain MVI strata",
        "ambiguity": "Adjacent sample records do not independently repeat MVI status; representative tar had matrix.mtx without observed feature/barcode companions",
        "source_cells": "46789 cells (Series summary; source-reported)",
        "raw_available": "YES, author-described raw cell count matrices",
        "filtered_available": "UNKNOWN; cell calling/filtering details need server-side inspection",
        "normalized_available": "NO",
        "role": "phase1_preferred_pending_format_check",
    },
    "GSE326201": {
        "matrix_type": "filtered_counts",
        "data_representation": "18 per-sample Cell Ranger filtered feature-barcode HDF5 files inside GSE326201_RAW.tar",
        "platform": "10x Chromium Next GEM 5-prime v2 / Illumina NovaSeq 6000; Cell Ranger 6.0; human hg38",
        "selection": "Unselected (no enrichment stated)",
        "matrix_url": ftp_series("GSE326201", "GSE326201_RAW.tar"),
        "matrix_size_bytes": "244039680",
        "fastq_status": "YES: 18 public SRA runs (SRR37818017-SRR37818034) found",
        "abundance": "YES for the eight confirmed paired patients; two patients have tumour-only libraries",
        "ambiguity": "No MVI field; source provides filtered matrices rather than unfiltered droplets",
        "source_cells": "Not stated at Series level; per-library filtered H5 files are supplied",
        "raw_available": "NO separate unfiltered raw matrix observed",
        "filtered_available": "YES",
        "normalized_available": "NO",
        "role": "phase1_preferred",
    },
    "GSE149614": {
        "matrix_type": "raw_counts",
        "data_representation": "Combined tab-delimited raw count matrix (165349783 bytes compressed) plus sample/cell metadata; normalized matrix is optional",
        "platform": "10x Chromium 3-prime v2 / Illumina NovaSeq 6000; Cell Ranger 2.2.0; hg38",
        "selection": "Unselected (no enrichment stated)",
        "matrix_url": ftp_series("GSE149614", "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"),
        "matrix_size_bytes": "165349783",
        "fastq_status": "EGA raw-data access cited as EGAS00001004468; no SRA hit in accession-term search",
        "abundance": "CONDITIONAL: eight paired primary tumour/normal pairs; exclude two PVTT and one lymph-node libraries",
        "ambiguity": "Public matrix includes non-target tissues and is a processed cell set; etiology comes from updated metadata",
        "source_cells": "71915 rows in updated cell metadata (source-reported public cell set)",
        "raw_available": "YES (count.txt.gz)",
        "filtered_available": "UNKNOWN; matrix is cell-level processed data",
        "normalized_available": "YES (normalized.txt.gz, not preferred)",
        "role": "phase1_preferred_after_nonTA_exclusion",
    },
    "GSE299340": {
        "matrix_type": "raw_counts",
        "data_representation": "10 per-sample integer Cell Ranger Matrix Market matrices with features/barcodes inside GSE299340_RAW.tar",
        "platform": "10x Chromium / Illumina HiSeq 2000 (SOFT; series notes HiSeq or NovaSeq); Cell Ranger 7.0.0 matrix metadata",
        "selection": "Unselected (no enrichment stated)",
        "matrix_url": ftp_series("GSE299340", "GSE299340_RAW.tar"),
        "matrix_size_bytes": "843857920",
        "fastq_status": "No SRA hit in accession-term search; GEO matrix archive available",
        "abundance": "YES, five paired tumour/noncancerous pairs; source matrices reflect prior cell-level QC",
        "ambiguity": "Source cell type is labelled hepatocyte in GEO; public matrix is already cell-called/processed and MVI is patient-level",
        "source_cells": "73707 cells retained after source QC (Series overall design; not project pre-QC)",
        "raw_available": "YES integer matrix; unfiltered droplet matrix not observed",
        "filtered_available": "LIKELY cell-called; not separately labelled in the Series record",
        "normalized_available": "NO separate normalized file observed",
        "role": "phase1_preferred_with_sourceQC_caveat",
    },
    "GSE290298": {
        "matrix_type": "normalized",
        "data_representation": "One processed CSV (923544964 bytes compressed) with normalized floating-point expression; eight gene-expression libraries plus two scTCR-only records",
        "platform": "10x Chromium / Illumina NovaSeq 6000; Cell Ranger 6.1.2; GRCh38",
        "selection": "Unselected (no enrichment stated)",
        "matrix_url": ftp_series("GSE290298", "GSE290298_All_processed_data.csv.gz"),
        "matrix_size_bytes": "923544964",
        "fastq_status": "No SRA hit in accession-term search; processed CSV available",
        "abundance": "CONDITIONAL/NOT PREFERRED: four paired gene-expression pairs, but no raw or filtered count matrix was observed",
        "ambiguity": "Two GEO samples are scTCR-only; normalized processed CSV is unsuitable as the default raw-count phase-1 input",
        "source_cells": "Not stated in GEO Series record",
        "raw_available": "NO",
        "filtered_available": "NO separately identified file",
        "normalized_available": "YES",
        "role": "optional_normalized_only",
    },
    "GSE202642": {
        "matrix_type": "raw_counts",
        "data_representation": "Combined Matrix Market integer matrix, 36601 x 115732 with 180199647 nonzeros, plus features/barcodes",
        "platform": "10x Chromium 3-prime v3 / Illumina NovaSeq; Cell Ranger 4.0.0; human reference not explicitly repeated in file header",
        "selection": "FACS",
        "matrix_url": ftp_series("GSE202642", "GSE202642_matrix.mtx.gz"),
        "matrix_size_bytes": "698853827",
        "fastq_status": "No SRA hit in accession-term search; GEO combined matrix available",
        "abundance": "NO for unbiased whole-tissue abundance: GEO protocol explicitly says FACS-sorted cells",
        "ambiguity": "Seven tumour and four adjacent tissues are listed, but patient-level pairing is not verified from GEO sample records",
        "source_cells": "115732 matrix barcodes; Matrix Market dimensions are source-observed",
        "raw_available": "YES integer matrix; source protocol says Cell Ranger-filtered barcodes",
        "filtered_available": "YES (Cell Ranger-filtered barcodes stated in protocol)",
        "normalized_available": "NO",
        "role": "cell_state_only_FACS_biased",
    },
}


DOWNLOAD_FILES = {
    "GSE282701": [("GSE282701_RAW.tar", ftp_series("GSE282701", "GSE282701_RAW.tar"), "856524800", "filtered_counts")],
    "GSE242889": [("GSE242889_RAW.tar", ftp_series("GSE242889", "GSE242889_RAW.tar"), "1376245760", "raw_counts")],
    "GSE326201": [("GSE326201_RAW.tar", ftp_series("GSE326201", "GSE326201_RAW.tar"), "244039680", "filtered_counts")],
    "GSE149614": [
        ("GSE149614_HCC.scRNAseq.S71915.count.txt.gz", ftp_series("GSE149614", "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"), "165349783", "raw_counts"),
        ("GSE149614_HCC.metadata.updated.txt.gz", ftp_series("GSE149614", "GSE149614_HCC.metadata.updated.txt.gz"), "489882", "metadata"),
    ],
    "GSE299340": [("GSE299340_RAW.tar", ftp_series("GSE299340", "GSE299340_RAW.tar"), "843857920", "raw_counts")],
    "GSE290298": [("GSE290298_All_processed_data.csv.gz", ftp_series("GSE290298", "GSE290298_All_processed_data.csv.gz"), "923544964", "normalized")],
    "GSE202642": [
        ("GSE202642_matrix.mtx.gz", ftp_series("GSE202642", "GSE202642_matrix.mtx.gz"), "698853827", "raw_counts"),
        ("GSE202642_features.tsv.gz", ftp_series("GSE202642", "GSE202642_features.tsv.gz"), "333437", "metadata"),
        ("GSE202642_barcodes.tsv.gz", ftp_series("GSE202642", "GSE202642_barcodes.tsv.gz"), "583228", "metadata"),
    ],
}


def parse_soft(path: Path) -> list[dict[str, object]]:
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"sample_id": line.split("=", 1)[1].strip()}
                continue
            if current is None or not line.startswith("!Sample_") or " = " not in line:
                continue
            key, value = line[8:].split(" = ", 1)
            current.setdefault(key, []).append(value)
    if current is not None:
        samples.append(current)
    return samples


def values(row: dict[str, object], key: str) -> list[str]:
    value = row.get(key, [])
    return [str(x) for x in value] if isinstance(value, list) else [str(value)]


def first(row: dict[str, object], key: str) -> str:
    return values(row, key)[0] if values(row, key) else ""


def characteristics(row: dict[str, object]) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in values(row, "characteristics_ch1"):
        if ":" in item:
            key, value = item.split(":", 1)
            result[key.strip().lower()] = value.strip()
    return result


def parse_gse149_metadata(path: Path | None) -> dict[str, int]:
    if path is None or not path.exists():
        return {}
    counts: Counter[str] = Counter()
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[row["sample"]] += 1
    return dict(counts)


def etiology_from_hbv_hcv(raw: str) -> str:
    parts = [x.strip() for x in raw.split("/")]
    if len(parts) != 2:
        return "unknown"
    hbv, hcv = parts
    if hbv == "+":
        return "HBV"
    if hcv == "+":
        return "HCV"
    if hbv == "-" and hcv == "-":
        return "NBNC"
    return "other"


def build_rows(soft_dir: Path, gse149_counts: dict[str, int]) -> list[dict[str, str]]:
    parsed: dict[str, list[dict[str, object]]] = {}
    for accession in ACCESSIONS:
        parsed[accession] = parse_soft(soft_dir / f"{accession}_family.soft.gz")

    # Resolve MVI from the tumor records in GSE242889 before assigning it to pairs.
    gse242_mvi: dict[str, str] = {}
    for row in parsed["GSE242889"]:
        title = first(row, "title")
        match = re.match(r"(\d+)T$", title)
        if match:
            source = first(row, "source_name_ch1").lower()
            gse242_mvi[match.group(1)] = (
                "positive" if "mvi present" in source else "negative" if "mvi absent" in source else "unknown"
            )

    rows: list[dict[str, str]] = []
    for accession in ACCESSIONS:
        info = DATASET[accession]
        raw_rows = parsed[accession]
        patient_counts: Counter[str] = Counter()
        tissue_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
        preliminary: list[dict[str, str]] = []

        for row in raw_rows:
            sample_id = first(row, "sample_id")
            title = first(row, "title")
            chars = characteristics(row)
            tissue = "unknown"
            patient_id = "unknown"
            paired_id = "unknown"
            etiology = "unknown"
            mvi = "unknown"
            sample_role = "gene_expression"
            phase1_role = info["role"]
            notes: list[str] = []

            if accession == "GSE282701":
                match = re.match(r"P(\d+),(tumor|adjacent),(.+)$", title, re.I)
                if match:
                    patient_id = f"P{match.group(1)}"
                    paired_id = patient_id
                    tissue = "Tumor" if match.group(2).lower() == "tumor" else "Adjacent"
                    etiology = "HBV" if match.group(3).lower() == "hbv-po" else "NBNC"
                else:
                    notes.append("title did not match expected P#,tumor/adjacent,HBV-* pattern")
            elif accession == "GSE242889":
                match = re.match(r"(\d+)(NT|T)$", title)
                if match:
                    case = match.group(1)
                    patient_id = f"case_{case}"
                    paired_id = patient_id
                    tissue = "Adjacent" if match.group(2) == "NT" else "Tumor"
                    mvi = gse242_mvi.get(case, "unknown")
                notes.append("patient identity is represented by the paired 1NT/1T-style sample name")
            elif accession == "GSE326201":
                match = re.match(r"HCC patient (\d+), ([^,]+), (non-tumour|tumour) tissue", title, re.I)
                if match:
                    patient_id = f"patient_{int(match.group(1)):03d}"
                    paired_id = patient_id
                    tissue = "Adjacent" if match.group(3).lower().startswith("non") else "Tumor"
                    etiology = {"HBV": "HBV", "HCV": "HCV", "NBNC": "NBNC"}.get(match.group(2), "unknown")
                else:
                    notes.append("title did not match expected patient/etiology/tissue pattern")
            elif accession == "GSE149614":
                match = re.search(r"patient:\s*(HCC\d+)", " || ".join(values(row, "characteristics_ch1")), re.I)
                if match:
                    patient_id = match.group(1)
                    paired_id = patient_id
                tissue_text = chars.get("tissue", "")
                if "primary tumor" in tissue_text.lower():
                    tissue = "Tumor"
                elif "adjacent non-tumor" in tissue_text.lower():
                    tissue = "Adjacent"
                else:
                    tissue = "Other"
                    phase1_role = "exclude_non_TA"
                etiology = {"HCC01": "HBV", "HCC02": "HBV", "HCC03": "HCV", "HCC04": "HCV", "HCC05": "NBNC", "HCC06": "HBV", "HCC07": "NBNC", "HCC08": "NBNC", "HCC09": "HBV", "HCC10": "HBV"}.get(patient_id, "unknown")
                if sample_id in gse149_counts:
                    notes.append(f"public cell metadata rows={gse149_counts[sample_id]}")
            elif accession == "GSE299340":
                patient_raw = chars.get("patient id", "")
                patient_match = re.search(r"patient\s+(\d+)", patient_raw, re.I)
                if patient_match:
                    patient_id = f"patient_{int(patient_match.group(1))}"
                    paired_id = patient_id
                tissue = "Tumor" if "primary hepatocellular carcinoma" in chars.get("tissue", "").lower() else "Adjacent" if "noncancerous" in chars.get("tissue", "").lower() else "unknown"
                etiology = etiology_from_hbv_hcv(chars.get("hbv/hcv infection", ""))
                if "microvascular invasion" in chars:
                    mvi = "negative" if chars["microvascular invasion"] == "-" else "positive" if chars["microvascular invasion"] == "+" else "unknown"
            elif accession == "GSE290298":
                match = re.match(r"(N|T)(\d+),", title)
                if match:
                    patient_id = f"patient_{match.group(2)}"
                    paired_id = patient_id
                    tissue = "Adjacent" if match.group(1) == "N" else "Tumor"
                if "scTCR-seq" in title:
                    sample_role = "TCR_only"
                    phase1_role = "exclude_TCR_only"
                    notes.append("GEO title explicitly identifies scTCR-seq; do not use as gene-expression matrix")
            elif accession == "GSE202642":
                tissue = "Tumor" if "hepatocellular carcinoma tissues" in title.lower() else "Adjacent" if "adjacent liver tissues" in title.lower() else "unknown"
                code_match = re.search(r"\[([^]]+)\]", title)
                if code_match:
                    patient_id = f"specimen_{code_match.group(1)}"
                notes.append("patient-level pairing not verified; sample title/specimen code retained verbatim")

            patient_counts[patient_id] += 1
            tissue_counts[patient_id][tissue] += 1
            preliminary.append(
                {
                    "dataset": accession,
                    "sample_id": sample_id,
                    "sample_title": title,
                    "patient_id": patient_id,
                    "tissue": tissue,
                    "paired_status": "unknown",
                    "paired_id": paired_id,
                    "etiology": etiology,
                    "MVI": mvi,
                    "platform": info["platform"],
                    "selection_strategy": info["selection"],
                    "matrix_type": info["matrix_type"] if sample_role != "TCR_only" else "unknown",
                    "data_representation": info["data_representation"] if sample_role != "TCR_only" else "TCR-only record; no phase-1 expression matrix",
                    "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sample_id}",
                    "matrix_url": info["matrix_url"] if sample_role != "TCR_only" else "",
                    "download_status": "NOT_REQUIRED" if sample_role == "TCR_only" else "DOWNLOAD_BLOCKED",
                    "fastq_status": info["fastq_status"],
                    "n_cells_pre_qc": str(gse149_counts.get(sample_id, "")) if accession == "GSE149614" else "",
                    "sample_role": sample_role,
                    "phase1_role": phase1_role,
                    "notes": "; ".join(notes),
                }
            )

        for item in preliminary:
            if item["sample_role"] == "TCR_only":
                item["paired_status"] = "paired" if tissue_counts[item["patient_id"]].get("Tumor", 0) and tissue_counts[item["patient_id"]].get("Adjacent", 0) else "unknown"
            elif item["dataset"] == "GSE202642":
                item["paired_status"] = "unknown"
            else:
                has_t = tissue_counts[item["patient_id"]].get("Tumor", 0) > 0
                has_a = tissue_counts[item["patient_id"]].get("Adjacent", 0) > 0
                item["paired_status"] = "paired" if has_t and has_a else "unpaired"
            if item["dataset"] == "GSE149614" and item["tissue"] == "Other":
                item["paired_status"] = "unknown"
        rows.extend(preliminary)
    return rows


def inventory(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for accession in ACCESSIONS:
        info = DATASET[accession]
        group = [r for r in rows if r["dataset"] == accession]
        target = [r for r in group if r["sample_role"] == "gene_expression" and r["tissue"] in {"Tumor", "Adjacent"}]
        patients: dict[str, set[str]] = defaultdict(set)
        for r in target:
            patients[r["patient_id"]].add(r["tissue"])
        paired = sum(1 for tissues in patients.values() if tissues == {"Tumor", "Adjacent"})
        n_patients = len(patients) if accession != "GSE202642" else "unknown"
        output.append(
            {
                "dataset": accession,
                "n_patients": str(n_patients),
                "n_libraries_total": str(len(group)),
                "n_expression_libraries": str(len(target)),
                "n_tumour_samples": str(sum(r["tissue"] == "Tumor" for r in target)),
                "n_adjacent_samples": str(sum(r["tissue"] == "Adjacent" for r in target)),
                "n_other_samples": str(sum(r["tissue"] == "Other" for r in group)),
                "n_confirmed_paired_patients": str(paired),
                "matrix_format": info["data_representation"],
                "phase1_matrix_type": info["matrix_type"],
                "matrix_size_approx_bytes": info["matrix_size_bytes"],
                "source_cells_or_barcodes": info["source_cells"],
                "raw_counts_available": info["raw_available"],
                "filtered_counts_available": info["filtered_available"],
                "normalized_available": info["normalized_available"],
                "fastq_sra_availability": info["fastq_status"],
                "enrichment_or_selection_bias": info["selection"],
                "whole_tissue_abundance_suitability": info["abundance"],
                "metadata_ambiguities": info["ambiguity"],
                "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "download_status": "DOWNLOAD_BLOCKED" if info["role"] != "optional_normalized_only" else "DOWNLOAD_BLOCKED",
                "server_download_path": f"/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/raw_data/{accession}/",
                "notes": "Server connection unavailable during Task 001; use resumable script after access is restored.",
            }
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft-dir", type=Path, required=True)
    parser.add_argument("--gse149-metadata", type=Path)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--download-manifest", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    args = parser.parse_args()

    rows = build_rows(args.soft_dir, parse_gse149_metadata(args.gse149_metadata))
    sample_columns = [
        "dataset", "sample_id", "sample_title", "patient_id", "tissue", "paired_status", "paired_id", "etiology", "MVI", "platform", "selection_strategy", "matrix_type", "data_representation", "source_url", "matrix_url", "download_status", "fastq_status", "n_cells_pre_qc", "sample_role", "phase1_role", "notes",
    ]
    write_tsv(args.sample_manifest, rows, sample_columns)

    download_rows: list[dict[str, str]] = []
    for accession in ACCESSIONS:
        info = DATASET[accession]
        for filename, url, size, file_type in DOWNLOAD_FILES[accession]:
            download_rows.append({
                "dataset": accession,
                "phase1_role": info["role"],
                "filename": filename,
                "remote_url": url,
                "remote_size_bytes": size,
                "file_type": file_type,
                "data_representation": info["data_representation"],
                "target_server_path": f"/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/raw_data/{accession}/{filename}",
                "download_status": "DOWNLOAD_BLOCKED",
                "checksum_path": "/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/checksums/task001_download_checksums.sha256",
                "log_path": f"/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/logs/task001_{accession}.log",
                "verify_plan": "curl --continue-at -; sha256sum; gzip -t or tar -tf where applicable; inspect matrix dimensions",
                "notes": "No server-side file downloaded in this run because configured SSH hosts were unavailable.",
            })
    download_columns = [
        "dataset", "phase1_role", "filename", "remote_url", "remote_size_bytes", "file_type", "data_representation", "target_server_path", "download_status", "checksum_path", "log_path", "verify_plan", "notes",
    ]
    write_tsv(args.download_manifest, download_rows, download_columns)

    inventory_columns = list(inventory(rows)[0])
    write_tsv(args.inventory, inventory(rows), inventory_columns, delimiter=",")


if __name__ == "__main__":
    main()
