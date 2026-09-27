#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--final-path", required=True)
    p.add_argument("--mapping", required=True)
    p.add_argument("--var", required=True)
    p.add_argument("--qs-validation", required=True)
    p.add_argument("--validation", required=True)
    p.add_argument("--block-cells", type=int, default=10000)
    return p.parse_args()


def decode(arr):
    arr = np.asarray(arr)
    if arr.dtype.kind in {"S", "O"}:
        return np.array([
            x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x)
            for x in arr
        ], dtype=object)
    return arr


def h5_index(group):
    key = group.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return decode(group[key][...])


args = parse_args()
source_path = Path(args.source)
out_path = Path(args.out)
final_path = Path(args.final_path)
mapping_path = Path(args.mapping)
var_path = Path(args.var)
validation_path = Path(args.validation)
qs_validation_path = Path(args.qs_validation)

if not source_path.exists():
    raise FileNotFoundError(source_path)
for p in (mapping_path, var_path, qs_validation_path):
    if not p.exists():
        raise FileNotFoundError(p)

out_path.parent.mkdir(parents=True, exist_ok=True)
validation_path.parent.mkdir(parents=True, exist_ok=True)
tmp_path = Path(str(out_path) + ".partial.h5ad")
shell_path = Path(str(out_path) + ".varshell.h5ad")
for p in (tmp_path, shell_path):
    if p.exists():
        p.unlink()

mapping = pd.read_csv(mapping_path)
var_meta = pd.read_csv(var_path).sort_values("canonical_index")
qs_val = pd.read_csv(qs_validation_path).iloc[0]

if mapping["raw_feature"].duplicated().any():
    raise RuntimeError("Feature mapping contains duplicated raw_feature entries")
if var_meta["canonical_symbol"].duplicated().any():
    raise RuntimeError("Canonical HGNC table contains duplicated symbols")
if len(mapping) != 68394:
    raise RuntimeError(f"Expected 68,394 source features, found {len(mapping)}")

new_symbols = var_meta["canonical_symbol"].astype(str).to_numpy()
n_new = len(new_symbols)
symbol_to_new = {g: i for i, g in enumerate(new_symbols)}

src = ad.read_h5ad(source_path, backed="r")
n_obs, n_old = src.shape
if (n_obs, n_old) != (1_490_852, 68_394):
    raise RuntimeError(f"Unexpected source H5AD shape: {src.shape}")

source_var_names = np.asarray(src.var_names.astype(str))
map_by_raw = mapping.set_index("raw_feature")
missing_map = [g for g in source_var_names if g not in map_by_raw.index]
if missing_map:
    raise RuntimeError(f"{len(missing_map)} source H5AD var names are absent from mapping; examples: {missing_map[:10]}")

old_to_new = np.full(n_old, -1, dtype=np.int64)
for i, gene in enumerate(source_var_names):
    row = map_by_raw.loc[gene]
    if bool(row["mapped_hgnc"]):
        sym = str(row["canonical_symbol"])
        old_to_new[i] = symbol_to_new[sym]

mapped_old = np.where(old_to_new >= 0)[0]
mapped_new = old_to_new[mapped_old]
if len(mapped_old) == 0:
    raise RuntimeError("No H5AD features mapped to HGNC")

with h5py.File(source_path, "r") as h5src:
    if "X" not in h5src or not isinstance(h5src["X"], h5py.Group):
        raise RuntimeError("Source H5AD X is not sparse-group encoded")
    source_dtype = h5src["X"]["data"].dtype
    source_x_encoding = h5src["X"].attrs.get("encoding-type", "")
    if isinstance(source_x_encoding, bytes):
        source_x_encoding = source_x_encoding.decode()
    if source_x_encoding != "csr_matrix":
        raise RuntimeError(f"Expected source CSR X, found {source_x_encoding}")

map_matrix = sp.csr_matrix(
    (
        np.ones(len(mapped_old), dtype=source_dtype),
        (mapped_old, mapped_new),
    ),
    shape=(n_old, n_new),
)

# Build a small AnnData shell solely to obtain spec-compliant var/varm/varp/layers encodings.
var_df = var_meta.set_index("canonical_symbol")[["hgnc_id", "locus_group", "locus_type"]].copy()
var_df.index = var_df.index.astype(str)
shell = ad.AnnData(
    X=sp.csr_matrix((0, n_new), dtype=source_dtype),
    obs=pd.DataFrame(index=pd.Index([], dtype=str)),
    var=var_df,
)
shell.write_h5ad(shell_path, compression="gzip")
del shell

source_total_counts = 0.0
hgnc_total_counts = 0.0
source_nnz = 0
hgnc_nnz = 0
global_offset = 0

with h5py.File(source_path, "r") as h5src, \
     h5py.File(shell_path, "r") as h5shell, \
     h5py.File(tmp_path, "w") as h5out:

    for k, v in h5src.attrs.items():
        h5out.attrs[k] = v

    # Copy observation-side structures unchanged.
    for key in ("obs", "obsm", "obsp", "uns"):
        if key in h5src:
            h5src.copy(key, h5out)

    # Copy new feature-side structures from the HGNC shell.
    for key in ("var", "varm", "varp", "layers"):
        if key in h5shell:
            h5shell.copy(key, h5out)

    gx = h5out.create_group("X")
    gx.attrs["encoding-type"] = "csr_matrix"
    gx.attrs["encoding-version"] = "0.1.0"
    gx.attrs["shape"] = np.asarray([n_obs, n_new], dtype=np.int64)

    element_chunk = 1_000_000
    data_ds = gx.create_dataset(
        "data",
        shape=(0,),
        maxshape=(None,),
        dtype=source_dtype,
        chunks=(element_chunk,),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
    )
    indices_ds = gx.create_dataset(
        "indices",
        shape=(0,),
        maxshape=(None,),
        dtype=np.int32,
        chunks=(element_chunk,),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
    )
    indptr_ds = gx.create_dataset(
        "indptr",
        shape=(n_obs + 1,),
        dtype=np.int64,
        chunks=(min(1_000_000, n_obs + 1),),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
    )
    indptr_ds[0] = 0

    for start in range(0, n_obs, args.block_cells):
        end = min(n_obs, start + args.block_cells)
        print(f"Harmonising H5AD rows {start:,}-{end:,} / {n_obs:,}", flush=True)

        block = src.X[start:end]
        if not sp.isspmatrix_csr(block):
            block = block.tocsr()
        block.sum_duplicates()
        block.sort_indices()

        source_total_counts += float(block.data.sum(dtype=np.float64))
        source_nnz += int(block.nnz)

        new_block = block @ map_matrix
        if not sp.isspmatrix_csr(new_block):
            new_block = new_block.tocsr()
        new_block.sum_duplicates()
        new_block.eliminate_zeros()
        new_block.sort_indices()

        hgnc_total_counts += float(new_block.data.sum(dtype=np.float64))
        hgnc_nnz += int(new_block.nnz)

        next_offset = global_offset + new_block.nnz
        data_ds.resize((next_offset,))
        indices_ds.resize((next_offset,))
        data_ds[global_offset:next_offset] = new_block.data.astype(source_dtype, copy=False)
        indices_ds[global_offset:next_offset] = new_block.indices.astype(np.int32, copy=False)
        indptr_ds[start + 1:end + 1] = global_offset + new_block.indptr[1:].astype(np.int64, copy=False)
        global_offset = next_offset

        del block, new_block

if global_offset != hgnc_nnz:
    raise RuntimeError("Final H5AD CSR offset mismatch")

try:
    src.file.close()
except Exception:
    pass
del src, map_matrix

# Structural validation before atomically promoting the file.
check = ad.read_h5ad(tmp_path, backed="r")
if check.shape != (1_490_852, n_new):
    raise RuntimeError(f"Harmonised H5AD shape mismatch: {check.shape}")
if not check.obs_names.is_unique:
    raise RuntimeError("Harmonised H5AD obs_names are not unique")
if not check.var_names.is_unique:
    raise RuntimeError("Harmonised H5AD var_names are not unique")
if not np.array_equal(np.asarray(check.var_names.astype(str)), new_symbols):
    raise RuntimeError("Harmonised H5AD var order differs from canonical HGNC order")
try:
    check.file.close()
except Exception:
    pass

with h5py.File(tmp_path, "r") as h5:
    if h5["X"].attrs.get("encoding-type", "") not in ("csr_matrix", b"csr_matrix"):
        raise RuntimeError("Harmonised H5AD X is not CSR")
    final_nnz = int(h5["X"]["data"].shape[0])
    if final_nnz != hgnc_nnz:
        raise RuntimeError("Harmonised H5AD nnz mismatch")
    obs_names = h5_index(h5["obs"])
    var_names = h5_index(h5["var"])
    if len(obs_names) != 1_490_852 or len(var_names) != n_new:
        raise RuntimeError("HDF5 axis lengths mismatch")

qs_source_total = float(qs_val["old_total_counts"])
qs_hgnc_total = float(qs_val["hgnc_total_counts"])
tol_old = max(1.0, abs(qs_source_total) * 1e-10)
tol_new = max(1.0, abs(qs_hgnc_total) * 1e-10)
if abs(source_total_counts - qs_source_total) > tol_old:
    raise RuntimeError(f"Source count sum disagrees with QS: H5AD={source_total_counts}, QS={qs_source_total}")
if abs(hgnc_total_counts - qs_hgnc_total) > tol_new:
    raise RuntimeError(f"HGNC count sum disagrees with QS: H5AD={hgnc_total_counts}, QS={qs_hgnc_total}")

os.replace(tmp_path, out_path)
if shell_path.exists():
    shell_path.unlink()

result = {
    "status": "VALIDATED",
    "source_h5ad": str(source_path),
    "harmonised_h5ad_temp": str(out_path),
    "harmonised_h5ad_final": str(final_path),
    "h5ad_size_bytes": out_path.stat().st_size,
    "n_obs": 1_490_852,
    "source_features": 68_394,
    "hgnc_features": n_new,
    "mapped_source_feature_rows": int(len(mapped_old)),
    "source_nnz": int(source_nnz),
    "hgnc_nnz": int(hgnc_nnz),
    "source_total_counts": source_total_counts,
    "hgnc_total_counts": hgnc_total_counts,
    "count_retention_fraction": hgnc_total_counts / source_total_counts,
    "obs_names_unique": True,
    "var_names_unique": True,
    "X_encoding": "csr_matrix",
    "block_cells": args.block_cells,
    "mapping_file": str(mapping_path),
    "var_file": str(var_path),
    "qs_cross_validation": "PASSED",
}
validation_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
