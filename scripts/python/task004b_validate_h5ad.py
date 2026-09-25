#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", required=True)
parser.add_argument("--expected-cells", type=int, default=1490852)
parser.add_argument("--expected-datasets", type=int, default=8)
parser.add_argument("--out", required=True)
args = parser.parse_args()

path = Path(args.h5ad)
if not path.exists():
    raise FileNotFoundError(path)

adata = ad.read_h5ad(path, backed="r")
required_obs = [
    "dataset", "sample_id", "patient_id", "tissue", "abundance_eligible",
    "project_broad_celltype", "neutrophil_confidence",
    "source_author_annotation", "qc_status", "annotation_status"
]
missing = [x for x in required_obs if x not in adata.obs.columns]
if missing:
    raise RuntimeError(f"Missing required AnnData obs fields: {missing}")

datasets = pd.Index(adata.obs["dataset"].astype(str).unique())
if adata.n_obs != args.expected_cells:
    raise RuntimeError(f"Cell count mismatch: {adata.n_obs} != {args.expected_cells}")
if len(datasets) != args.expected_datasets:
    raise RuntimeError(f"Dataset count mismatch: {len(datasets)} != {args.expected_datasets}")
if not adata.obs_names.is_unique:
    raise RuntimeError("AnnData obs_names are not globally unique")

first_obs = str(adata.obs_names[0])
last_obs = str(adata.obs_names[-1])
sha256 = hashlib.sha256(path.read_bytes()).hexdigest()

result = {
    "h5ad_path": str(path),
    "h5ad_size_bytes": path.stat().st_size,
    "sha256": sha256,
    "n_obs": int(adata.n_obs),
    "n_vars": int(adata.n_vars),
    "n_datasets": int(len(datasets)),
    "datasets": sorted(datasets.tolist()),
    "obs_names_unique": bool(adata.obs_names.is_unique),
    "first_obs_name": first_obs,
    "last_obs_name": last_obs,
    "required_obs_fields_present": True,
    "X_backed_type": str(type(adata.X)),
}
Path(args.out).write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
adata.file.close()
