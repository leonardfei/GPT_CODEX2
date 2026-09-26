#!/usr/bin/env python3
import argparse, csv, json, os
from collections import Counter
from pathlib import Path
import anndata as ad
import h5py
import numpy as np

p=argparse.ArgumentParser()
p.add_argument("--manifest", required=True); p.add_argument("--out", required=True)
p.add_argument("--validation", required=True); p.add_argument("--max-loaded-elems", type=int, default=50_000_000)
a=p.parse_args()
rows=list(csv.DictReader(open(a.manifest, encoding="utf-8")))
order=["GSE282701","GSE242889","GSE326201","GSE149614","GSE299340","CRA002308","nature_xue","in_house"]
by={r["dataset"]:r for r in rows}
if set(by)!=set(order): raise RuntimeError("Parts manifest dataset mismatch")
files={d:by[d]["h5ad_path"] for d in order}
for d,f in files.items():
    if not Path(f).exists(): raise FileNotFoundError(f"{d}: {f}")
expected_nnz=sum(int(by[d]["nnz"]) for d in order)
out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
tmp=Path(str(out)+".partial")
for x in (tmp,out):
    if x.exists(): x.unlink()
ad.experimental.concat_on_disk(files,tmp,axis="obs",join="outer",merge="first",uns_merge=None,label=None,index_unique=None,fill_value=0,pairwise=False,max_loaded_elems=a.max_loaded_elems)
os.replace(tmp,out)
lazy=ad.experimental.read_lazy(out)
if tuple(lazy.shape)!=(1490852,68394): raise RuntimeError(f"shape mismatch: {lazy.shape}")
try:
    if getattr(lazy,"file",None) is not None: lazy.file.close()
except Exception: pass

def dec(arr):
    arr=np.asarray(arr)
    if arr.dtype.kind in {"S","O"}:
        return np.array([v.decode() if isinstance(v,(bytes,bytearray)) else str(v) for v in arr],dtype=object)
    return arr
def idx(g):
    k=g.attrs.get("_index","_index"); k=k.decode() if isinstance(k,bytes) else k
    return dec(g[k][...])
def col(g,n):
    z=g[n]
    if isinstance(z,h5py.Dataset): return dec(z[...])
    enc=z.attrs.get("encoding-type",""); enc=enc.decode() if isinstance(enc,bytes) else enc
    if enc=="categorical":
        codes=np.asarray(z["codes"][...]); cats=dec(z["categories"][...]); out=np.empty(len(codes),dtype=object)
        miss=codes<0; out[miss]=None; good=~miss; out[good]=cats[codes[good]]; return out
    raise RuntimeError(f"Unsupported obs encoding {n}: {enc}")

with h5py.File(out,"r") as f:
    on=idx(f["obs"]); vn=idx(f["var"])
    if len(on)!=1490852 or len(vn)!=68394: raise RuntimeError("axis length mismatch")
    if len(set(on.tolist()))!=len(on): raise RuntimeError("duplicate obs_names")
    if len(set(vn.tolist()))!=len(vn): raise RuntimeError("duplicate var_names")
    req=["dataset","project_sample_id","project_patient_id","tissue","project_broad_celltype","neutrophil_confidence","source_author_annotation","annotation_status"]
    miss=[x for x in req if x not in f["obs"]]
    if miss: raise RuntimeError(f"missing obs: {miss}")
    ds=col(f["obs"],"dataset"); ss=col(f["obs"],"project_sample_id"); pp=col(f["obs"],"project_patient_id")
    tt=col(f["obs"],"tissue"); aa=col(f["obs"],"annotation_status")
    if len(set(x for x in ds if x is not None))!=8: raise RuntimeError("dataset count mismatch")
    if len(set(x for x in ss if x is not None))!=194: raise RuntimeError("sample count mismatch")
    if len(set(x for x in pp if x is not None))!=132: raise RuntimeError("patient count mismatch")
    tc=Counter(x for x in tt if x is not None)
    if tc.get("Tumor",0)!=1039293 or tc.get("Adjacent",0)!=451559: raise RuntimeError(f"tissue counts mismatch: {tc}")
    if set(x for x in aa if x is not None)!={"preliminary_unvalidated_task004"}: raise RuntimeError("annotation status mismatch")
    x=f["X"]
    if not isinstance(x,h5py.Group) or "data" not in x: raise RuntimeError("X is not sparse")
    enc=x.attrs.get("encoding-type",""); enc=enc.decode() if isinstance(enc,bytes) else enc
    nnz=int(x["data"].shape[0])
    if nnz!=expected_nnz: raise RuntimeError(f"nnz mismatch {nnz} != {expected_nnz}")

res={"status":"VALIDATED","h5ad_path":str(out),"h5ad_size_bytes":out.stat().st_size,"anndata_version":ad.__version__,"n_obs":1490852,"n_vars":68394,"n_datasets":8,"n_samples":194,"n_patients":132,"tumor_cells":1039293,"adjacent_cells":451559,"obs_names_unique":True,"var_names_unique":True,"X_encoding":enc,"nnz":nnz,"expected_nnz":expected_nnz,"join":"outer","fill_value":0}
Path(a.validation).write_text(json.dumps(res,indent=2),encoding="utf-8")
print(json.dumps(res,indent=2))
