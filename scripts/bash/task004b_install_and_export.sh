#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
R_LIB="${ROOT}/.task004b_Rlib"
PY_ENV="${ROOT}/.task004b_pyenv"
PARTS_DIR="${ROOT}/tmp/task004b_h5ad_parts"
LOG="${ROOT}/logs/task004b_export_formats.log"
mkdir -p "${R_LIB}" "${PARTS_DIR}" "${ROOT}/logs" "${ROOT}/results"
exec > >(tee -a "${LOG}") 2>&1
echo "=== Task 004b format export ==="; date; echo "ROOT=${ROOT}"
CHECKPOINT="${ROOT}/objects/task004_merge_review/checkpoint/HCC_TA_8datasets_merged_review_v1_checkpoint.rds"
[[ -s "${CHECKPOINT}" ]] || { echo "ERROR: missing checkpoint ${CHECKPOINT}" >&2; exit 2; }
FREE_KB=$(df -Pk "${ROOT}" | awk 'NR==2 {print $4}')
FREE_GB=$((FREE_KB / 1024 / 1024))
echo "Free disk space: ${FREE_GB} GB"
(( FREE_GB >= 60 )) || { echo "ERROR: need at least 60 GB free disk space" >&2; exit 3; }
export R_LIBS_USER="${R_LIB}"
Rscript "${ROOT}/scripts/R/task004b_install_export_deps.R" --project-root "${ROOT}" --lib "${R_LIB}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PY_VER=$("${PYTHON_BIN}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
echo "System Python: ${PY_VER}"
[[ -x "${PY_ENV}/bin/python" ]] || "${PYTHON_BIN}" -m venv "${PY_ENV}"
"${PY_ENV}/bin/python" -m pip install --upgrade pip setuptools wheel
PY_MAJOR=$("${PY_ENV}/bin/python" -c 'import sys; print(sys.version_info.major)')
PY_MINOR=$("${PY_ENV}/bin/python" -c 'import sys; print(sys.version_info.minor)')
if (( PY_MAJOR == 3 && PY_MINOR >= 12 )); then
  "${PY_ENV}/bin/pip" install "anndata==0.13.4" "h5py>=3.11" "scipy>=1.14" "pandas>=2.3"
elif (( PY_MAJOR == 3 && PY_MINOR == 11 )); then
  "${PY_ENV}/bin/pip" install "anndata>=0.12.19,<0.13" "h5py>=3.11" "scipy>=1.12" "pandas>=2.2"
else
  echo "ERROR: Python >=3.11 required; found ${PY_MAJOR}.${PY_MINOR}" >&2; exit 4
fi
echo "--- Export QS ---"
Rscript "${ROOT}/scripts/R/task004b_export_qs.R" --project-root "${ROOT}" --lib "${R_LIB}"
echo "--- Export cohort H5AD parts ---"
Rscript "${ROOT}/scripts/R/task004b_export_h5ad_parts.R" --project-root "${ROOT}" --lib "${R_LIB}" --out-dir "${PARTS_DIR}"
echo "--- On-disk H5AD concat ---"
"${PY_ENV}/bin/python" "${ROOT}/scripts/python/task004b_concat_h5ad.py" --manifest "${ROOT}/results/task004b_h5ad_parts_manifest.csv" --out "${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.h5ad" --validation "${ROOT}/results/task004b_h5ad_validation.json"
test -s "${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.qs"
test -s "${ROOT}/objects/HCC_TA_8datasets_merged_review_v1.h5ad"
grep -q VALIDATED "${ROOT}/results/task004b_qs_validation.csv"
grep -q '"status": "VALIDATED"' "${ROOT}/results/task004b_h5ad_validation.json"
echo "Both requested outputs validated."
rm -f "${PARTS_DIR}"/*.h5ad
echo "Checkpoint retained for rollback: ${CHECKPOINT}"
date
