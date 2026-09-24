# Task 003 — Prepare CRA002308, nature_xue and in_house for atlas inclusion

## Status
PENDING

## Goal
Audit and prepare the three newly uploaded cohorts so they can enter the same downstream annotation/integration framework as the five corrected-QC cohorts.

## Server location
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/raw_data/`

Expected cohort directories/names:
- CRA002308
- nature_xue
- in_house

Do not redownload these cohorts from the internet.

## Part A — inventory before processing
Recursively inventory each cohort directory before opening large files.

Record:
- path and filename;
- size;
- extension/file type;
- probable matrix/object type;
- sample identifiers encoded in filenames;
- presence of matrices/features/barcodes/H5/H5AD/RDS/RData/FASTQ/BAM;
- checksums for the primary files actually used.

Do not assume the uploaded representation matches the public repository representation.

## Part B — CRA002308

### Expected biological design
Seven HCC Tumor samples and seven matched adjacent/normal liver samples.

### Required audit
1. Identify all 14 expected biological samples and recover patient pairing.
2. Normalize tissue labels to:
   - T/Tumor -> `Tumor`
   - N/P/Normal/Adjacent/AL as supported by source metadata -> `Adjacent`
3. Verify the pair key from sample IDs. Do not invent pairing if filenames/metadata conflict.
4. Determine whether uploaded files are cell-called count matrices, Seurat/H5AD objects, BAM/FASTQ, or another representation.
5. Determine/record whether the study used immune-cell enrichment or any composition-altering selection. Assign `abundance_eligible` only after this is established.

### Processing
If cell-called count matrices are available:
- preserve source counts;
- calculate nFeature_RNA, nCount_RNA and percent.mt;
- apply the corrected Task 002 identity-independent secondary QC:
  - min_nFeature_RNA from all-source-cell q01 bounded to 100–300;
  - min_nCount_RNA from all-source-cell q01 bounded to 200–500;
  - max_percent_mt from all-source-cell q98 bounded to 20–30%;
- retain a `qc_status` field;
- perform the same high-confidence/broad-granulocyte retention audit used in corrected Task 002.

If only raw BAM/FASTQ or another non-count representation is present:
- do not fabricate count matrices;
- do not silently use a different quantification pipeline;
- stop CRA002308 processing at the audit stage and report the exact preprocessing required.

## Part C — nature_xue author-processed Seurat object

### Source
Xue R et al. Liver tumour immune microenvironment subtypes and neutrophil heterogeneity. Nature. 2022;612:141–147.

### Required handling
1. Identify and load the author-provided Seurat object.
2. Record Seurat object class/version compatibility, assays, layers/slots, reductions, cell count, feature count and metadata columns.
3. Preserve the original object unchanged.
4. Identify exact author metadata fields encoding:
   - species;
   - cancer type;
   - tissue/source;
   - sample/patient;
   - author broad/fine cell type or cluster.
5. Subset **human HCC only**.
6. Within HCC retain only:
   - Tumor
   - AL / adjacent liver
7. Exclude PB, mouse, ICC, CHC and other cancer/metastatic/non-HCC entities.
8. Preserve author annotations in `source_author_annotation` and additional namespaced metadata as needed.
9. Harmonize tissue to `Tumor` / `Adjacent`, with AL -> Adjacent.
10. Do not reapply the Task 002 QC thresholds blindly to this author-processed object.
11. If raw counts exist, preserve them. If not, explicitly set `counts_available=FALSE`; do not invent counts.
12. Record the object's selection/sampling design and assign `abundance_eligible` separately from atlas inclusion.

Create a subset/harmonized server object without modifying the original author object.

## Part D — in_house

### Known naming
- YJCA = HCC Tumor
- YJP = Adjacent

### Required audit
1. Inventory all uploaded sample files.
2. Recover patient/sample pairing from the naming scheme only after confirming a unique mapping.
3. Normalize:
   - YJCA -> Tumor
   - YJP -> Adjacent
4. Verify platform/matrix format and whether source matrices are raw, filtered/cell-called, or already processed.
5. Preserve in-house sample identifiers and create stable project patient IDs.
6. If count matrices are present, apply the same corrected Task 002 identity-independent QC and neutrophil-retention audit.
7. Record HBV/HCV/NBNC, MVI and other clinical fields only if available in uploaded metadata; use unknown rather than guessing.

## Harmonized metadata required for all three
At minimum:
- dataset
- sample_id
- patient_id
- tissue
- paired_status
- paired_id
- etiology
- MVI
- platform
- selection_strategy
- matrix_type
- counts_available
- qc_provenance
- qc_status
- abundance_eligible
- analysis_inclusion
- source_author_annotation when available

Cell IDs must be globally unique and should include dataset + sample provenance.

## Required outputs

### Git-tracked
- `results/task003_extension_file_inventory.csv`
- `results/task003_extension_sample_manifest.csv`
- `results/task003_extension_cohort_summary.csv`
- `results/task003_extension_qc_audit.csv`
- `results/task003_nature_xue_subset_audit.csv`
- `reports/task_003_report.md`

### Server-side
Prepared objects under:
`/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas/objects/task003_extension/`

Recommended names:
- `CRA002308_prepared.rds`
- `nature_xue_HCC_T_AL_prepared.rds`
- `in_house_prepared.rds`

Per-sample objects are acceptable for CRA002308/in_house if more appropriate.

## Completion criteria
- nature_xue human HCC Tumor/AL subset is prepared and author annotations preserved;
- CRA002308 and in_house are either prepared from compatible count/object inputs or clearly marked BLOCKED with exact input/preprocessing requirements;
- abundance eligibility is explicitly assigned and justified for each new cohort;
- no Task 004 annotation/integration is performed yet.

Stop after Task 003 and report.
