# Project Status

## Current task
Task 004c — HGNC feature harmonisation — READY TO EXECUTE

Task 004b is COMPLETED:
- 1,490,852 cells
- 68,394 raw union features
- 8 datasets
- 194 samples
- 132 patients
- QS and H5AD validated

## Why Task 004c
The 68,394-feature union is retained for archival/raw storage but should not be used directly as the cross-cohort integration feature universe.

Task 004c will:
- preserve all Task 004/004b source objects unchanged;
- download and checksum the current HGNC complete set;
- standardise exact symbols, version-stripped Ensembl IDs, Entrez IDs, previous symbols and unique aliases to current HGNC-approved symbols;
- refuse ambiguous mappings;
- sum raw counts when multiple source rows map to one approved HGNC symbol;
- create eight HGNC-harmonised counts+metadata cohort objects;
- quantify the strict 8/8 shared HGNC gene set and >=7/8 set.

Unmapped/custom/HBV features remain in source/raw objects but are excluded from the shared human integration universe.

## Next command
Execute task_004c.

Equivalent:
bash scripts/bash/task004c_harmonize_features.sh /data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas

Task 005 remains paused pending review of Task 004c.