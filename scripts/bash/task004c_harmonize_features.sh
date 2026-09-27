#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
export R_LIBS_USER="${ROOT}/.task004b_Rlib"
PYTHON="${ROOT}/.task004b_pyenv/bin/python"
cd "${ROOT}"

SOURCE_QS="${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.qs"
SOURCE_H5AD="${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.h5ad"
OUT_QS="${ROOT}/objects/HCC_TA_8datasets_merged_review_HGNC_v1.qs"
OUT_H5AD="${ROOT}/objects/HCC_TA_8datasets_merged_review_HGNC_v1.h5ad"

[[ -s "${SOURCE_QS}" ]] || { echo "ERROR: missing source QS" >&2; exit 2; }
[[ -s "${SOURCE_H5AD}" ]] || { echo "ERROR: missing source H5AD" >&2; exit 2; }
[[ -x "${PYTHON}" ]] || { echo "ERROR: Task 004b Python environment is missing" >&2; exit 2; }

FREE_KB=$(df -Pk "${ROOT}" | awk 'NR==2 {print $4}')
FREE_GB=$((FREE_KB / 1024 / 1024))
echo "Free disk: ${FREE_GB} GB"
if (( FREE_GB < 50 )); then
  echo "ERROR: require at least 50 GB free for the new HGNC QS/H5AD outputs." >&2
  exit 3
fi

echo "=== Harmonise merged QS to HGNC-approved symbols ==="
Rscript scripts/R/task004c_harmonize_features.R \
  --project-root "${ROOT}" \
  --r-lib "${ROOT}/.task004b_Rlib" \
  --source-qs "${SOURCE_QS}" \
  --out-qs "${OUT_QS}"

echo "=== Harmonise merged H5AD directly on disk ==="
"${PYTHON}" scripts/python/task004c_harmonize_merged_h5ad.py \
  --source "${SOURCE_H5AD}" \
  --out "${OUT_H5AD}" \
  --mapping "${ROOT}/results/task004c_merged_feature_mapping.csv" \
  --var "${ROOT}/results/task004c_hgnc_var.csv" \
  --qs-validation "${ROOT}/results/task004c_qs_validation.csv" \
  --validation "${ROOT}/results/task004c_h5ad_validation.json" \
  --block-cells 10000

test -s "${OUT_QS}"
test -s "${OUT_H5AD}"
grep -q VALIDATED "${ROOT}/results/task004c_qs_validation.csv"
grep -q '"status": "VALIDATED"' "${ROOT}/results/task004c_h5ad_validation.json"
echo "Task 004c direct merged-object harmonisation completed."
echo "Original v1 QS/H5AD retained unchanged."
