#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
R_LIB="${ROOT}/.task004b_Rlib"
PY_ENV="${ROOT}/.task004b_pyenv"
PY_WHEEL_DIR="${ROOT}/tmp/task004b_py_wheels"
PARTS_DIR="${ROOT}/tmp/task004b_h5ad_parts"
LOG="${ROOT}/logs/task004b_export_formats.log"
mkdir -p "${R_LIB}" "${PARTS_DIR}" "${ROOT}/logs" "${ROOT}/results"
export PATH="${PY_ENV}:${PY_ENV}/bin:${ROOT}/tmp/r_env/bin:${PATH}"
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
SYS_PY_MAJOR=$("${PYTHON_BIN}" -c 'import sys; print(sys.version_info.major)')
SYS_PY_MINOR=$("${PYTHON_BIN}" -c 'import sys; print(sys.version_info.minor)')
echo "System Python: ${SYS_PY_MAJOR}.${SYS_PY_MINOR}"

if (( SYS_PY_MAJOR == 3 && SYS_PY_MINOR >= 11 )); then
  if [[ ! -x "${PY_ENV}/bin/python" ]]; then
    if ! "${PYTHON_BIN}" -m venv "${PY_ENV}"; then
      rm -rf "${PY_ENV}"
      CONDA_BIN="$(command -v conda || true)"
      if [[ -z "${CONDA_BIN}" ]]; then
        for candidate in "/data/lf_data/miniconda3/bin/conda" "/data/miniconda3/bin/conda" "/opt/conda/bin/conda"; do
          if [[ -x "${candidate}" ]]; then CONDA_BIN="${candidate}"; break; fi
        done
      fi
      [[ -n "${CONDA_BIN}" ]] || { echo "ERROR: venv failed and conda was not found." >&2; exit 4; }
      "${CONDA_BIN}" create -y -p "${PY_ENV}" -c conda-forge --override-channels python=3.12 pip cmake
    fi
  fi
  if [[ ! -x "${PY_ENV}/bin/cmake" ]]; then
    "${PY_ENV}/bin/pip" install cmake
  fi
else
  CONDA_BIN="$(command -v conda || true)"
  if [[ -z "${CONDA_BIN}" ]]; then
    for candidate in "${HOME}/miniconda3/bin/conda" "/data/lf_data/miniconda3/bin/conda" "/opt/conda/bin/conda"; do
      if [[ -x "${candidate}" ]]; then CONDA_BIN="${candidate}"; break; fi
    done
  fi
  if [[ -z "${CONDA_BIN}" ]]; then
    echo "ERROR: system Python is <3.11 and conda was not found." >&2
    exit 4
  fi
  if [[ -d "${PY_ENV}" && ! -x "${PY_ENV}/bin/python" ]]; then rm -rf "${PY_ENV}"; fi
  if [[ ! -x "${PY_ENV}/bin/python" ]] || ! "${PY_ENV}/bin/python" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3,11) else 1)'; then
    rm -rf "${PY_ENV}"
    "${CONDA_BIN}" create -y -p "${PY_ENV}" python=3.12 pip
  fi
fi

"${PY_ENV}/bin/python" -m pip --version
PY_MAJOR=$("${PY_ENV}/bin/python" -c 'import sys; print(sys.version_info.major)')
PY_MINOR=$("${PY_ENV}/bin/python" -c 'import sys; print(sys.version_info.minor)')
if ! "${PY_ENV}/bin/python" -c 'import anndata, h5py, scipy, pandas, pydantic; assert anndata.__version__ == "0.13.4"; assert tuple(map(int, h5py.__version__.split(".")[:2])) >= (3, 11); assert tuple(map(int, scipy.__version__.split(".")[:2])) >= (1, 14); assert tuple(map(int, pandas.__version__.split(".")[:2])) >= (2, 3); assert tuple(map(int, pydantic.__version__.split(".")[:2])) >= (2, 13)' 2>/dev/null; then
  if [[ -d "${PY_WHEEL_DIR}" ]] && compgen -G "${PY_WHEEL_DIR}/*.whl" >/dev/null; then
    echo "Installing Python export dependencies from local wheelhouse: ${PY_WHEEL_DIR}"
    if (( PY_MAJOR == 3 && PY_MINOR >= 12 )); then
      "${PY_ENV}/bin/pip" install --no-index --find-links "${PY_WHEEL_DIR}" "anndata==0.13.4" "h5py>=3.11" "scipy>=1.14" "pandas>=2.3" "pydantic>=2.13.5"
    elif (( PY_MAJOR == 3 && PY_MINOR == 11 )); then
      "${PY_ENV}/bin/pip" install --no-index --find-links "${PY_WHEEL_DIR}" "anndata>=0.12.19,<0.13" "h5py>=3.11" "scipy>=1.12" "pandas>=2.2"
    else
      echo "ERROR: isolated Python >=3.11 could not be prepared; found ${PY_MAJOR}.${PY_MINOR}" >&2
      exit 4
    fi
  elif (( PY_MAJOR == 3 && PY_MINOR >= 12 )); then
    "${PY_ENV}/bin/pip" install "anndata==0.13.4" "h5py>=3.11" "scipy>=1.14" "pandas>=2.3" "pydantic>=2.13.5"
  elif (( PY_MAJOR == 3 && PY_MINOR == 11 )); then
    "${PY_ENV}/bin/pip" install "anndata>=0.12.19,<0.13" "h5py>=3.11" "scipy>=1.12" "pandas>=2.2"
  else
    echo "ERROR: isolated Python >=3.11 could not be prepared; found ${PY_MAJOR}.${PY_MINOR}" >&2
    exit 4
  fi
else
  echo "Python export dependencies already satisfy the task requirements."
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
