#!/usr/bin/env python3
"""Apply the Task 003 metadata-only correction to small audit tables.

The prepared objects are corrected by scripts/R/task003_correct_metadata.R.
This companion step keeps the Git-tracked matrix manifest and QC audit in
agreement with those object-level metadata changes, without touching counts or
cell pass/fail decisions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CRA_SELECTION = (
    "flow-sorted live nucleated cells after doublet exclusion; "
    "no lineage-specific immune enrichment documented; paired within-cohort "
    "abundance sensitivity analysis allowed"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.results_root.resolve()

    manifest_path = root / "task003_extension_sample_manifest_matrix_cohorts.csv"
    manifest = pd.read_csv(manifest_path)
    cra = manifest["dataset"].eq("CRA002308")
    manifest.loc[cra, "abundance_eligible"] = "CONDITIONAL"
    manifest.loc[cra, "selection_strategy"] = CRA_SELECTION
    manifest.to_csv(manifest_path, index=False)

    qc_path = root / "task003_extension_qc_audit.csv"
    qc = pd.read_csv(qc_path)
    cra_qc = qc["dataset"].eq("CRA002308")
    qc.loc[cra_qc, "abundance_eligible"] = "CONDITIONAL"
    qc.loc[cra_qc, "selection_strategy"] = CRA_SELECTION
    qc.to_csv(qc_path, index=False)

    print(f"Updated CRA002308 metadata rows: manifest={int(cra.sum())}, qc={int(cra_qc.sum())}")


if __name__ == "__main__":
    main()
