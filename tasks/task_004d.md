# Task 004d — Per-sample scDblFinder filtering and corrected broad cell-type annotation

## Status
PENDING — APPROVED FOR CODEX EXECUTION

## Goal
Starting from the HGNC-harmonised eight-cohort merged QS object, identify and remove predicted doublets with scDblFinder on a per-capture/sample basis, then rebuild a singlet-only merged object and perform a corrected broad cell-type annotation using standard transcriptomic preprocessing, clustering, cluster markers, canonical lineage programs, and the Xue author annotation as a reference anchor.

This task replaces the preliminary Task 004 broad labels as the project working annotation. Task 004 labels remain available only for audit/comparison.

Do not execute Task 005 integration in this task.

## Canonical input
Use exactly:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_merged_review_v1.qs`

Preflight expectation:
- 1,490,852 cells;
- 38,025 unique HGNC-approved symbols;
- 8 datasets;
- 194 globally unique `project_sample_id` values;
- 132 `project_patient_id` values;
- one RNA counts layer;
- globally unique cell IDs.

If this exact path is absent, do not silently substitute another object. Report the available candidate path(s) and stop.

## Phase A — scDblFinder doublet calling

### Biological unit
Run doublet detection independently for each `project_sample_id`.

Treat `project_sample_id` as the independent 10x/library capture unit unless the input metadata explicitly proves otherwise. Do not run one scDblFinder model on all 1.49 million cells together, and do not use patient ID as the capture unit.

### Input expression
Use raw RNA counts only.

The merged QS uses the custom `task004b_chunked_counts` representation. Implement/restore the necessary `dim`, `dimnames`, and subset methods and extract one sample at a time into a standard sparse matrix. Never materialise the full 1.49M × 38,025 matrix as a single `dgCMatrix`.

### scDblFinder version
Use a current compatible Bioconductor release of:
- `scDblFinder`;
- `SingleCellExperiment`;
- `BiocParallel`.

Record exact versions.

Require `scDblFinder >= 1.24.8` if xgboost >=3.1 is installed. Prefer the current stable release compatible with the server R/Bioconductor version.

Explicitly set `dbr.per1k = 0.008`; do not rely on an old package default.

If packages are missing, first check whether the server has working package access. The server must not be assumed to have internet. If installation is impossible, stop with an exact dependency blocker; do not substitute DoubletFinder/Scrublet without approval.

### scDblFinder parameters
Default per sample:

```r
set.seed(44000 + sample_index)
sce <- scDblFinder::scDblFinder(
    sce,
    clusters = TRUE,
    dbr = NULL,
    dbr.per1k = 0.008,
    dbr.sd = NULL,
    iter = 2,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
)
```

Small-sample rule:
- if a sample has <500 retained cells, use `clusters = NULL` to avoid unstable fast clustering;
- do not merge a small sample with another sample for doublet detection.

Do not manually tune the final threshold sample-by-sample unless the algorithm fails. Any deviation must be documented.

### Required per-cell metadata
Collect for every original cell:
- `cell_id`;
- dataset;
- sample/patient/tissue identifiers;
- `project_sample_id`;
- `scDblFinder.score`;
- `scDblFinder.class`;
- any additional scDblFinder origin/cluster fields returned by the installed version;
- preliminary Task 004 broad label for audit only.

### Filtering
A cell is removed from the singlet object only when:

`scDblFinder.class == "doublet"`

Do not add a secondary UMI/nFeature doublet cutoff.

Do not rescue a predicted doublet solely because it was previously labelled neutrophil.

### Required doublet QC
Report overall, per dataset, and per sample:
- source cells;
- singlets;
- predicted doublets;
- doublet fraction;
- median/IQR `scDblFinder.score`;
- nCount/nFeature distributions in singlets vs doublets;
- Tumor vs Adjacent doublet fractions;
- number/fraction of preliminary Task 004 neutrophils called doublet;
- neutrophil retention fraction after filtering;
- preliminary broad-cell-type composition among predicted doublets.

Add review flags for unusually high/low called-doublet fractions, but do not silently change thresholds to make rates look typical.

For CRA002308, record that experimental/source preprocessing already reported doublet exclusion, but still run the same computational audit because no per-cell known-doublet labels are available.

For nature_xue, preserve all author labels and explicitly state that scDblFinder is an additional project-level computational filter on an author-processed object.

### Doublet outputs
Write:
- `results/task004d_scdblfinder_cell_calls.csv.gz`
- `results/task004d_scdblfinder_by_sample.csv`
- `results/task004d_scdblfinder_by_dataset.csv`
- `results/task004d_neutrophil_doublet_retention.csv`
- `figures/task004d_scdblfinder_qc.pdf`

The all-cell calls file is the audit trail for removed cells.

## Phase B — build the singlet-only merged object

Create a new object; do not overwrite the Task 004c input.

Primary intermediate:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_singlets_v1.qs`

Requirements:
- all and only scDblFinder singlets;
- 38,025 HGNC feature space retained unless a technical reason requires otherwise;
- original cell IDs and metadata preserved;
- source author annotations preserved;
- doublet-detection provenance added to `misc`;
- raw counts preserved;
- no normalization/batch correction stored as a replacement for raw counts.

Validate exact singlet count against the cell-call table.

## Phase C — corrected broad cell-type annotation

### Core principle
Do not reuse the Task 004 per-cell marker-score labels as the final annotation.

Annotation must be cluster-based and supported by:
1. standard transcriptomic preprocessing and clustering;
2. cluster-level marker evidence;
3. canonical lineage marker programs;
4. Xue author labels as an independent reference anchor where available.

No cross-dataset batch integration is performed here.

### Analysis feature universe
Use:
`results/task004c_shared_hgnc_features_8of8.txt`

as the candidate shared feature universe for cross-cohort-comparable PCA/clustering features.

Within each dataset:
- subset singlets by dataset;
- LogNormalize (`scale.factor=10000`);
- select ~3,000 HVGs;
- restrict integration/comparable PCA features to the Task 004c 8/8 HGNC shared set;
- exclude mitochondrial and ribosomal genes from PCA/HVG features if they dominate the selected set;
- ScaleData on selected features only;
- PCA with 30 PCs initially;
- inspect elbow/variance but keep 30 PCs as the reproducible default unless clear failure;
- FindNeighbors on PCs 1:30;
- FindClusters with Leiden/Louvain resolution 0.4 and 0.8 for sensitivity;
- use resolution 0.8 as the initial annotation partition unless it fragments coherent broad lineages.

Do not regress out tissue, patient, etiology, neutrophil state, cell cycle, or mitochondrial fraction merely to improve mixing.

### Scalability
For datasets >250,000 singlets, full-cell normalization/HVG calculation is allowed, but if full-cell graph clustering is infeasible under available RAM, use an explicit Seurat v5 sketch workflow:
- leverage-score sketch up to 50,000 cells;
- cluster/annotate the sketch;
- project cluster labels back to every full-data cell;
- retain all cells in the final object.

This is a documented computational sketch, not silent biological downsampling.

### Broad annotation classes
Use these final broad labels:

- `Hepatocyte/Epithelial`
- `T_cell`
- `NK_cell`
- `B_cell`
- `Plasma_cell`
- `Monocyte/Macrophage`
- `Neutrophil`
- `Dendritic_cell`
- `Mast_cell`
- `Endothelial`
- `Fibroblast/Mesenchymal`
- `Uncertain/Mixed`

Do not call cells malignant at this stage. Malignancy requires a later CNV/tumour-cell analysis.

Do not use `Cycling` as a lineage label; store cycling as a state/flag while retaining the inferred parent lineage when possible.

### Canonical marker programs
Use coherent multi-marker evidence, including at minimum:

- Hepatocyte/Epithelial: `ALB, APOA1, APOA2, TTR, ASGR1, KRT8, KRT18, EPCAM, KRT19, KRT7`
- T cell: `CD3D, CD3E, TRAC, CD2, LTB`
- NK cell: `NKG7, GNLY, KLRD1, PRF1, CTSW`, with absent/low CD3 program
- B cell: `CD79A, CD79B, MS4A1, CD37, CD74, CD22, CD19`
- Plasma cell: `MZB1, JCHAIN, XBP1, DERL3, SDC1, IGKC`
- Monocyte/Macrophage: `LST1, TYROBP, FCER1G, CTSS, LILRB1, C1QC, APOC1, SPP1`
- Neutrophil: `FCGR3B, CSF3R, CXCR2, FPR1, FCAR, S100A8, S100A9, NAMPT, MCEMP1, ELANE, MPO`
- Dendritic cell: `FCER1A, CD1C, CLEC10A, CLEC9A, XCR1, GZMB, TCF4`
- Mast cell: `TPSAB1, TPSB2, KIT, MS4A2, HDC, CPA3`
- Endothelial: `PECAM1, VWF, EMCN, KDR, RAMP2, ENG, ESM1`
- Fibroblast/Mesenchymal: `COL1A1, COL1A2, COL3A1, DCN, LUM, COL6A1, COL6A2, PDGFRA, FAP, THY1`

A single marker is never sufficient for a cluster label.

### Cluster marker evidence
For every cluster, produce:
- top positive markers;
- average expression/detection of each canonical program;
- top and second-best lineage program;
- margin between them;
- fraction of cells expressing >=2 markers from the proposed lineage.

Use pseudobulk/average-expression summaries when full FindAllMarkers is computationally prohibitive. Broad annotation does not require cell-level p-values for every gene.

### Xue reference anchor
For `nature_xue`:
- preserve `source_author_annotation`;
- tabulate every distinct author label;
- create and commit `config/task004d_xue_author_to_broad.tsv`;
- map fine author labels to the broad classes above;
- leave unclear labels as `Uncertain/Mixed` rather than guessing.

For each Xue-derived cluster, report the author-label consensus and compatibility with cluster marker evidence.

For non-Xue datasets, use the Xue broad labels as a reference anchor at the cluster level (e.g. centroid correlation/reference transfer), but the final call must remain compatible with canonical markers. Do not copy a reference label when lineage markers contradict it.

### Annotation decision rule
Assign a broad label only if cluster-level evidence is coherent.

At minimum, require one of:
- dominant canonical program plus >=2 compatible lineage markers and no strong incompatible lineage program; or
- >=70% Xue/reference broad-label consensus plus compatible canonical markers.

If two incompatible lineages are both strong, or reference and canonical evidence conflict materially, label `Uncertain/Mixed` and flag for review.

For `Neutrophil`, require coherent granulocyte evidence including at least one core marker among `FCGR3B, CSF3R, FPR1, FCAR, ELANE, MPO` at cluster level; S100A8/S100A9 alone are insufficient.

### Final per-cell annotation metadata
Add:
- `broad_celltype_v2`
- `broad_annotation_confidence` = high / medium / low / uncertain
- `broad_annotation_basis`
- `broad_cluster_id`
- `broad_cluster_resolution`
- `broad_reference_label`
- `broad_reference_score`
- `task004_preliminary_label` copied from the old label for audit
- `annotation_status = "task004d_broad_annotation_v2"`

### Broad annotation QC
Produce:
- overall and per-dataset broad cell counts/fractions;
- Tumor/Adjacent counts by broad lineage;
- Xue author-label vs v2 crosswalk/confusion table;
- Task 004 preliminary vs v2 crosswalk;
- cluster marker table;
- uncertain/mixed cluster list;
- neutrophil marker coherence audit;
- UMAP/PCA review plots by dataset, tissue, cluster, v2 broad label, and source author label for Xue.

Do not interpret pooled cell fractions as biological abundance comparisons.

## Final objects

Create:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_singlets_broad_v1.qs`

and, if the existing Task 004c H5AD export machinery can be reused safely without changing scientific results:

`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/merge/HCC_TA_8datasets_singlets_broad_v1.h5ad`

The H5AD is requested for Scanpy interoperability and must contain raw counts in `X` or a clearly documented raw-count layer plus the v2 annotation in `obs`.

Do not overwrite the Task 004c merged input.

## Required tracked outputs
- `results/task004d_scdblfinder_cell_calls.csv.gz`
- `results/task004d_scdblfinder_by_sample.csv`
- `results/task004d_scdblfinder_by_dataset.csv`
- `results/task004d_neutrophil_doublet_retention.csv`
- `results/task004d_broad_celltype_counts.csv`
- `results/task004d_cluster_annotation.csv`
- `results/task004d_cluster_markers.csv.gz`
- `results/task004d_xue_author_vs_v2.csv`
- `results/task004d_task004_vs_v2.csv`
- `results/task004d_uncertain_clusters.csv`
- `config/task004d_xue_author_to_broad.tsv`
- `figures/task004d_scdblfinder_qc.pdf`
- `figures/task004d_broad_annotation_qc.pdf`
- `reports/task_004d_report.md`

## Validation
Final report must state:
- exact scDblFinder/package versions;
- exact parameters and seeds;
- input/output checksums and sizes;
- cells before doublet filtering;
- predicted doublets and fraction;
- singlets retained;
- per-dataset/sample doublet rates;
- preliminary-neutrophil retention;
- final v2 broad cell counts;
- number/fraction `Uncertain/Mixed`;
- Xue author-label concordance at broad level;
- whether sketching was used for any dataset;
- all deviations/errors and how they were resolved.

## Stop rule
Stop after Task 004d final objects and QC validate.

Do not execute Task 005 integration.

