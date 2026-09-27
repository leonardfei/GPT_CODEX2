#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
export R_LIBS_USER="${ROOT}/.task004b_Rlib"
cd "${ROOT}"
FREE_KB=$(df -Pk "${ROOT}" | awk 'NR==2 {print $4}')
FREE_GB=$((FREE_KB / 1024 / 1024))
echo "Free disk: ${FREE_GB} GB"
if (( FREE_GB < 20 )); then
  echo "ERROR: require at least 20 GB free disk for harmonised cohort objects and audit outputs." >&2
  exit 3
fi
Rscript scripts/R/task004c_harmonize_features.R \
  --project-root "${ROOT}" \
  --r-lib "${ROOT}/.task004b_Rlib"
