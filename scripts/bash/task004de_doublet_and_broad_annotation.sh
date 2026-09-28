#!/usr/bin/env bash
set -euo pipefail
ROOT="$1"
if [[ -z "$ROOT" ]]; then
  ROOT="/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas"
fi
cd "$ROOT"

SOURCE="$ROOT/objects/merge/HCC_TA_8datasets_merged_review_v1.qs"
FILTERED="$ROOT/objects/merge/HCC_TA_8datasets_merged_scdblfinder_filtered_v1.qs"
ANNOTATED="$ROOT/objects/merge/HCC_TA_8datasets_merged_scdblfinder_filtered_broad_v1.qs"

[[ -s "$SOURCE" ]] || { echo "ERROR: missing source QS: $SOURCE" >&2; exit 2; }

export R_LIBS_USER="$ROOT/.task004de_Rlib:$ROOT/.task004b_Rlib"

Rscript scripts/R/task004de_install_deps.R "$ROOT"

echo "=== Task 004d: per-sample scDblFinder ==="
Rscript scripts/R/task004d_scdblfinder_filter.R \
  --project-root "$ROOT" \
  --source-qs "$SOURCE" \
  --out-qs "$FILTERED"

grep -q VALIDATED "$ROOT/results/task004d_scdblfinder_validation.csv"
[[ -s "$FILTERED" ]]

echo "=== Task 004e: broad annotation ==="
Rscript scripts/R/task004e_broad_annotation.R \
  --project-root "$ROOT" \
  --source-qs "$FILTERED" \
  --out-qs "$ANNOTATED" \
  --shared-features "$ROOT/results/task004c_shared_hgnc_features_8of8.txt" \
  --sketch-cells 50000 \
  --projection-block 2000 \
  --resolution 0.6 \
  --seed 40500

grep -q VALIDATED "$ROOT/results/task004e_broad_annotation_validation.csv"
[[ -s "$ANNOTATED" ]]

echo "Task 004d/e completed."
echo "Filtered:  $FILTERED"
echo "Annotated: $ANNOTATED"
