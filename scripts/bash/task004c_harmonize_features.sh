#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
PROMOTE_ONLY="${2:-}"
export R_LIBS_USER="${ROOT}/.task004b_Rlib"
export PATH="${ROOT}/tmp/r_env/bin:${PATH}"
PYTHON="${ROOT}/.task004b_pyenv/bin/python"
cd "${ROOT}"
mkdir -p "${ROOT}/logs"
LOG="${ROOT}/logs/task004c_harmonize_features.log"
exec > >(tee -a "${LOG}") 2>&1
echo "=== Task 004c HGNC feature harmonisation ==="
date

QS="${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.qs"
H5AD="${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.h5ad"
TMP_QS="${QS}.task004c_tmp"
TMP_H5AD="${H5AD}.task004c_tmp.h5ad"
QS_BACKUP="${QS}.task004c_backup"
H5AD_BACKUP="${H5AD}.task004c_backup"

[[ -s "${QS}" ]] || { echo "ERROR: missing source QS" >&2; exit 2; }
[[ -s "${H5AD}" ]] || { echo "ERROR: missing source H5AD" >&2; exit 2; }
[[ -x "${ROOT}/tmp/r_env/bin/Rscript" ]] || { echo "ERROR: Task 004b R runtime is missing" >&2; exit 2; }
[[ -x "${PYTHON}" ]] || { echo "ERROR: Task 004b Python environment is missing" >&2; exit 2; }
[[ ! -e "${QS_BACKUP}" && ! -e "${H5AD_BACKUP}" ]] || {
  echo "ERROR: Task 004c backup path exists; inspect it before retrying." >&2
  exit 2
}
if [[ "${PROMOTE_ONLY}" == "--promote-only" ]]; then
  echo "=== Reusing previously validated temporary QS/H5AD files ==="
  test -s "${TMP_QS}"
  test -s "${TMP_H5AD}"
  grep -q VALIDATED "${ROOT}/results/task004c_qs_validation.csv"
  grep -q '"status": "VALIDATED"' "${ROOT}/results/task004c_h5ad_validation.json"
else
  rm -f "${TMP_QS}" "${TMP_H5AD}"

  FREE_KB=$(df -Pk "${ROOT}" | awk 'NR==2 {print $4}')
  FREE_GB=$((FREE_KB / 1024 / 1024))
  echo "Free disk: ${FREE_GB} GB"
  if (( FREE_GB < 50 )); then
    echo "ERROR: require at least 50 GB free for validated temporary HGNC replacements." >&2
    exit 3
  fi

  echo "=== Build and validate temporary HGNC QS ==="
  Rscript scripts/R/task004c_harmonize_features.R \
    --project-root "${ROOT}" \
    --r-lib "${ROOT}/.task004b_Rlib" \
    --source-qs "${QS}" \
    --out-qs "${TMP_QS}" \
    --final-qs "${QS}"

  echo "=== Build and validate temporary HGNC H5AD ==="
  "${PYTHON}" scripts/python/task004c_harmonize_merged_h5ad.py \
    --source "${H5AD}" \
    --out "${TMP_H5AD}" \
    --final-path "${H5AD}" \
    --mapping "${ROOT}/results/task004c_merged_feature_mapping.csv" \
    --var "${ROOT}/results/task004c_hgnc_var.csv" \
    --qs-validation "${ROOT}/results/task004c_qs_validation.csv" \
    --validation "${ROOT}/results/task004c_h5ad_validation.json" \
    --block-cells 10000
fi

test -s "${TMP_QS}"
test -s "${TMP_H5AD}"
grep -q VALIDATED "${ROOT}/results/task004c_qs_validation.csv"
grep -q '"status": "VALIDATED"' "${ROOT}/results/task004c_h5ad_validation.json"

echo "=== Both temporary files validated; atomically overwrite merged v1 files ==="
rollback_replacement() {
  status=$?
  if (( status != 0 )); then
    echo "ERROR: replacement interrupted; restoring any staged original files." >&2
    if [[ -e "${QS_BACKUP}" ]]; then
      mv -f "${QS_BACKUP}" "${QS}" || echo "ERROR: could not restore QS backup ${QS_BACKUP}" >&2
    fi
    if [[ -e "${H5AD_BACKUP}" ]]; then
      mv -f "${H5AD_BACKUP}" "${H5AD}" || echo "ERROR: could not restore H5AD backup ${H5AD_BACKUP}" >&2
    fi
  fi
  return "${status}"
}
trap rollback_replacement EXIT

mv "${QS}" "${QS_BACKUP}"
mv "${H5AD}" "${H5AD_BACKUP}"
mv "${TMP_QS}" "${QS}"
mv "${TMP_H5AD}" "${H5AD}"

test -s "${QS}"
test -s "${H5AD}"
trap - EXIT
rm -f "${QS_BACKUP}" "${H5AD_BACKUP}"
echo "Task 004c completed: original merged v1 QS/H5AD paths now contain HGNC-standardised data."
echo "No separate HGNC-named duplicate files were retained."
