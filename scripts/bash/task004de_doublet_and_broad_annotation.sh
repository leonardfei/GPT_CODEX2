#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
cd "${ROOT}"

SOURCE="${ROOT}/objects/merge/HCC_TA_8datasets_merged_review_v1.qs"
SINGLETS="${ROOT}/objects/merge/HCC_TA_8datasets_singlets_v1.qs"
FINAL="${ROOT}/objects/merge/HCC_TA_8datasets_singlets_broad_v1.qs"

[[ -s "${SOURCE}" ]] || {
  echo "ERROR: exact approved source QS is missing: ${SOURCE}" >&2
  echo "Available merge candidates:" >&2
  find "${ROOT}/objects/merge" -maxdepth 1 -type f \( -name '*.qs' -o -name '*.h5ad' \) -print >&2 || true
  exit 2
}

export R_LIBS_USER="${ROOT}/.task004d_Rlib:${ROOT}/.task004b_Rlib"

echo "=== Task 004d dependency check/install ==="
Rscript scripts/R/task004de_install_deps.R "${ROOT}"

echo "=== Task 004d Phase A/B: per-sample scDblFinder and singlet object ==="
Rscript scripts/R/task004d_scdblfinder_filter.R \
  --project-root "${ROOT}" \
  --source-qs "${SOURCE}" \
  --out-qs "${SINGLETS}" \
  --seed 44000

grep -q VALIDATED "${ROOT}/results/task004d_scdblfinder_validation.csv"
[[ -s "${SINGLETS}" ]]

echo "=== Task 004d Phase C: corrected broad annotation ==="
Rscript scripts/R/task004e_broad_annotation.R \
  --project-root "${ROOT}" \
  --source-qs "${SINGLETS}" \
  --out-qs "${FINAL}" \
  --shared-features "${ROOT}/results/task004c_shared_hgnc_features_8of8.txt" \
  --xue-map "${ROOT}/config/task004d_xue_author_to_broad.tsv" \
  --sketch-cells 50000 \
  --projection-block 2000 \
  --seed 40500

grep -q VALIDATED "${ROOT}/results/task004d_broad_annotation_validation.csv"
[[ -s "${FINAL}" ]]

echo "=== Task 004d final validation report ==="
Rscript scripts/R/task004d_finalize_report.R "${ROOT}"

[[ -s "${ROOT}/reports/task_004d_report.md" ]]
[[ -s "${ROOT}/results/task004d_object_checksums.csv" ]]

echo "Task 004d completed and validated."
echo "Singlets: ${SINGLETS}"
echo "Final broad-annotated singlets: ${FINAL}"
echo "Task 005 remains paused."
