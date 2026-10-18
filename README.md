# A neuron-to-astrocyte interleukin-6 relay drives schizophrenia-like behaviors

This repository contains the custom analysis code and processed inputs used in
the associated manuscript. It covers human dorsolateral prefrontal cortex
(DLPFC) transcriptomics, human peripheral-blood expression, mouse medial
prefrontal cortex (mPFC) expression, human and mouse single-nucleus RNA
sequencing (snRNA-seq), and improved in situ sequencing (IISS).

## Release status

This archive is ready for submission-stage review. Public-dataset identifiers
are fixed where available. Accession information for the newly generated mouse
mPFC snRNA-seq and human small-RNA/cohort VII mRNA-seq datasets is available
from the corresponding author upon request during peer review.

## Repository layout

This archive is distributed in two parts. The analysis code (`codes/`) is hosted on
GitHub, and the processed inputs (`source_data/`) are deposited at Mendeley Data
(DOI: [insert on deposit]). All paths in this README are relative to a reconstructed
archive root, that is a directory containing both `codes/` and `source_data/`. To run
the scripts, download `source_data/` from Mendeley Data and place it next to `codes/`.

```text
.
├── codes/
│   ├── 01_bulk_transcriptomics/
│   ├── 02_human_blood_expression/
│   ├── 03_mouse_expression/
│   ├── 04_snRNAseq/
│   │   ├── human_batiuk/
│   │   │   └── environment/sessionInfo_step5.txt
│   │   ├── human_hbcc/
│   │   │   └── environment/sessionInfo_step5.txt
│   │   └── mouse_mPFC/
│   │       └── environment/sessionInfo_step5.txt
│   ├── 05_IISS/
│   ├── 06_CAMK2A_AUCell_validation/
│   │   └── environment/run_info_step6.txt
│   ├── README.md
│   ├── environment.yml
│   ├── requirements.txt
│   └── R_PACKAGE_VERSIONS.tsv

├── source_data/
    ├── IISS/
    │   ├── processed_spatial/
    │   └── README.md
    ├── snRNAseq/
    │   ├── Batiuk/
    │   │   └── step6_CAMK2A_validation/
    │   ├── HBCC/
    │   ├── Mouse_mPFC/
    │   ├── batiuk_sample_manifest.csv
    │   ├── hbcc_donor_filter_manifest.csv
    │   └── mouse_mPFC_snRNAseq_sample_manifest.csv
    ├── bulk_DLPFC_meta_analysis_input.csv
    ├── cohort_VII_IL6_let7b_expression_data.xlsx
    ├── cohorts_I_IV_miRNA_volcano_plot_data.xlsx
    └── mouse_mPFC_Il6_let7b_qPCR_data.tsv
```

## Data access

| Dataset                                        | Identifier and input                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ---------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Batiuk neuron-enriched DLPFC snRNA-seq         | Processed data: Zenodo record [6921620](https://doi.org/10.5281/zenodo.6921620), version 2. Download `snRNA-seq_and_spatial_transcriptomics.zip` and use the enclosed `snRNA-seq_raw_countmatrices.RDS` and `annotations_final.RDS`. Raw data: EGA study `EGAS00001006495`. Original code: Zenodo record [6583332](https://doi.org/10.5281/zenodo.6583332).                                                                                                                                                                |
| PsychAD HBCC DLPFC snRNA-seq                   | CZ CELLxGENE collection `84ce6837-548d-4a1f-919f-0bc0d9a3952f`; dataset `5c97eeeb-7e52-44b3-b010-b832b1f5424c`; version `d27fb144-f105-46c2-b36f-f51421f74e4e`. The versioned H5AD is available at `https://datasets.cellxgene.cziscience.com/d27fb144-f105-46c2-b36f-f51421f74e4e.h5ad`.                                                                                                                                                                                                                                  |
| Mouse mPFC snRNA-seq                           | Accession information is available from the corresponding author upon request during peer review.                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Human small-RNA and cohort VII mRNA sequencing | Accession information is available from the corresponding author upon request during peer review.                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Bulk DLPFC cohorts                             | Controlled-access data were obtained from the PsychENCODE Phase I Capstone Collection (`syn12080241`). Gene-level RSEM files were obtained from CMC Freeze 1 (`syn8241760`), CMC_HBCC Freeze 2 (`syn10623034`) and LIBD__szControl Freeze 2 (`syn10623021`). Diagnostic and file-mapping metadata are available from `syn12118229.1` and `syn8466658.1`; the identifier crosswalk is `syn21978133.1`. Stage 1 accepts matrices assembled from these files. The aggregate Figure 1b input is included under `source_data/`. |
| IISS and figure-level inputs                   | Cell-by-gene matrices, figure-level inputs, DAPI reference images, gene-labelled signal coordinates, cell centroids, cell-boundary polygons, auxiliary mask images and a sample manifest are included under `source_data/IISS/`. Integer gene values in the four analysis matrices are assigned rolling circle product counts. Raw sequencing-cycle images and the proprietary barcode library are not available to the authors.                                                                                           |

To download the versioned HBCC object:

```bash
mkdir -p inputs/hbcc
curl -L \
  https://datasets.cellxgene.cziscience.com/d27fb144-f105-46c2-b36f-f51421f74e4e.h5ad \
  -o inputs/hbcc/HBCC_Cohort.h5ad
```

The HBCC scripts start from five cell-class Seurat objects and a global UMAP
table prepared from this H5AD. The conversion retains raw counts, donor ID,
diagnosis, published cell-type annotations and `X_umap`, as described in the
Methods. Routine H5AD-to-Seurat conversion is not reimplemented here.
`00_build_lightweight_object.R` takes over from those five objects and merges
them into the single lightweight object that step 5 consumes.

## Environments

### Python

Four figure-generating modules are written in Python: the two peripheral-blood
modules, the mouse mPFC correlation and the IISS UMAP. The processed-IISS
validator uses only the Python standard library. The figure modules were run
with Python 3.12.10,
NumPy 2.4.3, pandas 2.3.3, SciPy 1.18.0, scikit-learn 1.8.0, umap-learn 0.5.12,
numba 0.64.0, llvmlite 0.46.0, statsmodels 0.14.6, seaborn 0.13.2,
openpyxl 3.1.5 and Matplotlib 3.10.8. `codes/requirements.txt` and
`codes/environment.yml` pin these and their transitive dependencies to the
same versions.

```bash
conda env create --file codes/environment.yml
conda activate il6-relay-analysis
```

### R

The large snRNA-seq workflows were run under R 4.5.0 or 4.5.1. Each workflow's
`environment/sessionInfo_step5.txt` is the sanitized `sessionInfo()` record of
the step 5 run that produced the reported figures and statistics. No separate
record is retained for steps 1 to 4. The three files describe three separate
runtimes, so their package versions differ.
`codes/R_PACKAGE_VERSIONS.tsv` provides the package versions captured from
those records, the locally validated small R modules and the documented
preprocessing workflows. No `renv.lock` is supplied because a lockfile was not
generated and validated in each large-object runtime.

The snRNA-seq scripts draw per-cell point layers with `geom_point` and do not
depend on `ggrastr` or `Cairo`. Point layers are therefore vector graphics,
which increases PDF size for panels with many cells.

The R modules use Seurat, SeuratObject, DoubletFinder, harmony, scMayoMap,
AUCell, SummarizedExperiment, BiocParallel, glmmTMB, Matrix, metafor, dplyr,
tidyr, tidyselect, purrr, tibble, stringr, readr, readxl, ggplot2, ggrepel,
scales, patchwork, gridExtra, openxlsx and edgeR, in addition to the base R packages.
Install the versions recorded for the relevant workflow before rerunning a
large snRNA-seq analysis.

## Small analysis modules

Create an output directory and run commands from the archive root.

### Bulk DLPFC meta-analysis

This analysis runs in two stages because the donor-level PsychENCODE inputs are
controlled-access and cannot be redistributed. Stage 1 derives the aggregate
cohort statistics and Supplementary Figure 1a-c from the restricted data.
Stage 2 fits the meta-analysis and draws Figure 1b from the aggregate table
released here. Stage 2 alone reproduces Figure 1b.

Stage 1 supports the original count-matrix workflow. Assemble one gene-by-donor
RSEM count matrix for each cohort and a metadata table with columns `name` and
`diagnosis`. The three matrices are converted to CPM with `edgeR::cpm`, after
which IL6 (`ENSG00000136244`) is transformed as `log10(CPM + 1)`. The same
metadata file can be used for all three cohorts, as in the original analysis:

```bash
Rscript codes/01_bulk_transcriptomics/01a_cohort_effect_sizes.R \
  --cmc-counts  <CMC_genecounts.matrix> \
  --hbcc-counts <CMC_HBCC_genecounts.matrix> \
  --libd-counts <LIBD_scControl_genecounts.matrix> \
  --metadata    <group.csv> \
  --output source_data/bulk_DLPFC_meta_analysis_input.csv \
  --donor-output-dir outputs/bulk_DLPFC/donor_tables \
  --plot-output-dir outputs/bulk_DLPFC
```

Alternatively, investigators who have already completed that preprocessing may
provide one delimited donor-level table per cohort. Each table must contain
`group` (`SCZ` or `CTR`) and `expr_norm` (`log10(CPM + 1)`-transformed IL6):

```bash
Rscript codes/01_bulk_transcriptomics/01a_cohort_effect_sizes.R \
  --cmc  <CMC_IL6_donor_level.csv> \
  --hbcc <HBCC_IL6_donor_level.csv> \
  --libd <LIBD_IL6_donor_level.csv> \
  --output source_data/bulk_DLPFC_meta_analysis_input.csv \
  --plot-output-dir outputs/bulk_DLPFC
```

For each cohort, the effect size is the log2 ratio of the group mean
`log10(CPM + 1)`-transformed IL6 expression. Following the original analysis,
the standard error is the standard deviation of transformed IL6 expression
across all retained donors, pooled over both groups, divided by the square root
of the total donor number. The P value is from a two-sided Wilcoxon rank-sum
test of SCZ against CTR. The Supplementary Figure 1a-c labels are generated from
these tests and display CMC and CMC_HBCC to four decimal places and LIBD as
`P < 2e-11` when the computed value is below that reporting threshold. Cohort
confidence limits use the 1.96 multiplier; pooled limits come from the fitted
model. Individual-level P values cannot be recomputed from the released
aggregate table.

Stage 2 runs from the released aggregate table:

```bash
Rscript codes/01_bulk_transcriptomics/01b_meta_analysis.R \
  --input source_data/bulk_DLPFC_meta_analysis_input.csv \
  --output-dir outputs/bulk_DLPFC
```

Required fields are `cohort`, `n_schizophrenia`, `n_control`,
`effect_size_log2_expression_ratio`, `standard_error`, `ci_95_lower`,
`ci_95_upper`, `weight_percent` and `p_value`. The optional `p_value_bound`
field records whether a P value is exact or an upper bound; the LIBD entry is
`upper` because that cohort is reported as P < 2e-11. The script recomputes the
restricted maximum-likelihood random-effects estimate, verifies it against the
supplied pooled row and against the supplied cohort confidence limits, and
writes the Figure 1b PDF and statistics table.

### Cohorts I-IV miRNA volcano plots

```bash
python codes/02_human_blood_expression/cohorts_I_IV_volcano_plots.py \
  --input source_data/cohorts_I_IV_miRNA_volcano_plot_data.xlsx \
  --output-dir outputs/cohorts_I_IV
```

The workbook contains sheets `Cohort I` through `Cohort IV`, each with `miRNA`,
`log2 fold change`, `-log10(P value)` and `Direction`.

Each sheet is already restricted to the differentially expressed miRNAs of that
contrast, defined as an absolute log2 fold change of at least 0.5 with P below
0.05 in both DESeq2 and edgeR. Selection is therefore not repeated in the
script: the script assigns colours, applies the display transformation and
draws the panels. Points below the plotted P = 0.05 guide line do not occur
because non-significant miRNAs are not part of the released table.

The horizontal axis is transformed for display only, as
`sign(x) * log2(|x| + 1) * 3`, implemented in `stretch_fc`. Tick labels show the
untransformed log2 fold change. The transformation does not alter
differential-expression values or thresholds. Outputs support Supplementary
Figure 4a.

### Cohort VII correlation

```bash
python codes/02_human_blood_expression/cohort_VII_IL6_let7b_correlation.py \
  --input source_data/cohort_VII_IL6_let7b_expression_data.xlsx \
  --output-dir outputs/cohort_VII
```

The `Expression data` sheet supplies group labels and matched let-7b-5p and
`IL6` `log2(CPM + 1)` values. The script applies the manuscript's
Shapiro-Wilk-guided Pearson-or-Spearman rule and writes Figure 2f and its
statistics table.

### Mouse mPFC correlation

```bash
python codes/03_mouse_expression/mouse_mPFC_Il6_let7b_correlation.py \
  --input source_data/mouse_mPFC_Il6_let7b_qPCR_data.tsv \
  --output-dir outputs/mouse_mPFC_correlation
```

Required fields are `sample_id`, `group`, `let7b5p_relative_expression` and
`Il6_relative_expression`. Outputs support Figure 2i.

### IISS cell-type summary and composition

```bash
IISS_SOURCE_DATA_DIR=source_data/IISS \
IISS_OUTPUT_DIR=outputs/IISS \
Rscript codes/05_IISS/iiss_celltype_summary_and_figure3.R
```

This workflow reads `con2.txt`, `con10.txt`, `ko3.txt` and `ko11.txt`. Each
matrix requires `cell.name`, `celltype` and gene-level assigned-count columns.
It writes the two Figure 3l composition PDFs, the Figure 3m marker dot plot and
their Source Data.

### IISS UMAP

```bash
python codes/05_IISS/iiss_umap.py \
  --input-dir source_data/IISS \
  --output-dir outputs/IISS
```

This workflow uses the same four matrices, a fixed seed and exactly five
cell-type classes. It writes Figure 3k and its cell-level Source Data.

### IISS Il6 and Il6st summaries

```bash
Rscript codes/05_IISS/iiss_IL6_IL6ST_heatmaps.R \
  --input source_data/IISS/IISS_IL6_IL6ST_cell_type_summary.xlsx \
  --output-dir outputs/IISS
```

The analysis is descriptive at sample level and writes Figure 3n,
Supplementary Figure 9a and their Source Data.

### IISS processed spatial data validation

```bash
python codes/05_IISS/validate_processed_spatial_data.py \
  --input-dir source_data/IISS \
  --output outputs/IISS/processed_spatial_validation.csv
```

This low-memory validator streams each gene-labelled signal table, checks cell
identifiers against the centroid and cell-boundary files, and confirms that the
counts for every retained cell and gene aggregate exactly to `con2.txt`,
`con10.txt`, `ko3.txt` and `ko11.txt`. File definitions, the four-sample mapping
and the commercial-data limitations are documented in
`source_data/IISS/README.md`.

## Environment variables

Every location and tunable used by the R workflows is read from an environment
variable. Variables marked required cause the script to stop with an explicit
error when unset.

| Variable                                                                                | Used by                                              | Required | Meaning                                                                                                                 |
| --------------------------------------------------------------------------------------- | ---------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------------------- |
| `SEURAT_RDS_PATH`                                                                       | Batiuk 02, 04, 05; HBCC 04, 05; mouse 04, 05; CAMK2A | yes      | Seurat object consumed by that step                                                                                     |
| `AUCELL_CHUNK_NCELLS`                                                                   | Batiuk 04; HBCC 04                                   | no       | Cells per AUCell chunk; alters memory use only, not within-cell ranks                                                   |
| `PLOT_FONT_FAMILY`                                                                      | both bulk scripts; both IISS R modules               | no       | Plot font family; defaults to Helvetica                                                                                 |
| `BATIUK_COUNTS_RDS`                                                                     | Batiuk 01                                            | yes      | Released Batiuk raw count matrices                                                                                      |
| `BATIUK_ANNOTATION_RDS`                                                                 | Batiuk 01                                            | no       | Released Batiuk annotation object                                                                                       |
| `BATIUK_SAMPLE_MANIFEST`                                                                | Batiuk 01, 03                                        | yes      | `source_data/snRNAseq/batiuk_sample_manifest.csv`                                                                       |
| `BATIUK_OUTPUT_DIR`                                                                     | Batiuk 01-04                                         | yes      | Workflow output root                                                                                                    |
| `BATIUK_RESULT_ROOT`                                                                    | Batiuk 05; CAMK2A                                    | yes      | Analysis root written by step 4                                                                                         |
| `BATIUK_EXTRA_TABLE_DIRS`                                                               | Batiuk 05                                            | no       | Colon- or semicolon-separated extra table directories                                                                   |
| `HBCC_RDS_AST`, `HBCC_RDS_EXCITNEU`, `HBCC_RDS_INHIBNEU`, `HBCC_RDS_MIC`, `HBCC_RDS_OL` | HBCC 00                                              | yes      | The five cell-class Seurat objects                                                                                      |
| `HBCC_GLOBAL_UMAP`                                                                      | HBCC 00                                              | yes      | Global UMAP table with `cell`, `UMAP_1`, `UMAP_2`                                                                       |
| `HBCC_LIGHTWEIGHT_OUTPUT`                                                               | HBCC 00                                              | yes      | Path for the merged lightweight object                                                                                  |
| `KEEP_GENES`                                                                            | HBCC 00                                              | no       | Comma-separated genes retained in the lightweight object                                                                |
| `HBCC_BASE_PATH`                                                                        | HBCC 04, 05                                          | no       | Workflow base directory; an alternative to setting `HBCC_RESULT_ROOT` directly                                          |
| `HBCC_RESULT_ROOT`                                                                      | HBCC 04, 05                                          | yes      | Analysis root shared by all five step 4 runs                                                                            |
| `HBCC_METADATA_PATH`                                                                    | HBCC 05                                              | yes      | `source_data/snRNAseq/hbcc_donor_filter_manifest.csv`                                                                   |
| `HBCC_MIN_AGE_YEARS`                                                                    | HBCC 05                                              | no       | Minimum control donor age; 21 in the reported analysis                                                                  |
| `HBCC_MIN_CELLS_PER_CT`                                                                 | HBCC 05                                              | no       | Minimum nuclei per donor per cell type; 50 in the reported analysis                                                     |
| `HBCC_EXCLUDE_COMORBID`                                                                 | HBCC 05                                              | no       | Whether to exclude donors flagged for exclusionary comorbidity                                                          |
| `HBCC_STEP5_LIGHT_MODE`                                                                 | HBCC 05                                              | no       | Skips the marker dot plot and the raw-count IL6 UMAP                                                                    |
| `HBCC_EXTRA_TABLE_DIRS`                                                                 | HBCC 05                                              | no       | Colon- or semicolon-separated extra table directories                                                                   |
| `MOUSE_INPUT_DIR`                                                                       | mouse 01                                             | yes      | Directory holding one subdirectory per sequencing library; each holds its 10x matrix under `filtered_cell_gene_matrix/` |
| `MOUSE_OUTPUT_DIR`                                                                      | mouse 01-04                                          | yes      | Workflow output root                                                                                                    |
| `MOUSE_RESULT_ROOT`                                                                     | mouse 05                                             | yes      | Analysis root written by step 4                                                                                         |
| `IISS_SOURCE_DATA_DIR`                                                                  | IISS composition module; IISS UMAP                   | no       | Defaults to the released `source_data/IISS`; the UMAP module also accepts `--input-dir`                                 |
| `IISS_OUTPUT_DIR`                                                                       | IISS composition module; IISS UMAP                   | no       | Output directory; the UMAP module also accepts `--output-dir`                                                           |

## snRNA-seq workflows

Large workflows write timestamped analysis roots. Record the root printed by
step 4 and pass it unchanged to step 5. Do not substitute newer package
versions for UMAP or object-generating steps without validating coordinates,
cell identities and figure outputs.

Every input and output location is taken from the environment variables listed
below. Each script stops with an explicit error when a required variable is
missing, no location is inferred from the local filesystem, and the scripts do
not change the working directory of the calling session.

### Batiuk workflow

The released cohort contains 24 sample matrices, including 14 control samples
from 14 donors and 10 schizophrenia samples from 9 donors. MB8 and MB8-2 are
retained as separate released analysis samples, matching the analyzed
structure.

```bash
export BATIUK_COUNTS_RDS=inputs/batiuk/snRNA-seq_raw_countmatrices.RDS
export BATIUK_ANNOTATION_RDS=inputs/batiuk/annotations_final.RDS
export BATIUK_SAMPLE_MANIFEST=source_data/snRNAseq/batiuk_sample_manifest.csv
export BATIUK_OUTPUT_DIR=outputs/snRNA/batiuk

Rscript codes/04_snRNAseq/human_batiuk/01_qc.R
Rscript codes/04_snRNAseq/human_batiuk/02_harmony.R
Rscript codes/04_snRNAseq/human_batiuk/03_annotation.R

export SEURAT_RDS_PATH=outputs/snRNA/batiuk/2-harmony_ALL/res_0.8/dim_35_res_0.8_harmony_Celltypes.annotated3.rds
Rscript codes/04_snRNAseq/human_batiuk/04_aucell_scoring.R

export BATIUK_RESULT_ROOT=outputs/snRNA/batiuk/results/human_IL6_<TIMESTAMP>
Rscript codes/04_snRNAseq/human_batiuk/05_postprocess_figures_statistics.R
```

`AUCELL_CHUNK_NCELLS` can be set to alter memory use without changing
within-cell ranks. `BATIUK_EXTRA_TABLE_DIRS` is optional and accepts
colon-separated directories only when validated AUCell tables are distributed
across analysis roots.

### CAMK2A AUCell validation

The positive-control analysis for Supplementary Figure 16 uses the final
Batiuk object and the corresponding step 4 tables.

```bash
export SEURAT_RDS_PATH=outputs/snRNA/batiuk/2-harmony_ALL/res_0.8/dim_35_res_0.8_harmony_Celltypes.annotated3.rds
export BATIUK_RESULT_ROOT=outputs/snRNA/batiuk/results/human_IL6_<TIMESTAMP>
Rscript codes/06_CAMK2A_AUCell_validation/06_CAMK2A_AUCell_validation.R
```

The five panel-level aggregate tables used for Supplementary Figure 16 are
provided under `source_data/snRNAseq/Batiuk/step6_CAMK2A_validation/`.

### HBCC workflow

Prepare five Seurat objects for astrocytes, excitatory neurons, inhibitory
neurons, microglia and oligodendrocytes plus a global UMAP table with columns
`cell`, `UMAP_1` and `UMAP_2`.

```bash
export HBCC_RDS_AST=inputs/hbcc/HBCC_Ast.rds
export HBCC_RDS_EXCITNEU=inputs/hbcc/HBCC_ExcitNeu.rds
export HBCC_RDS_INHIBNEU=inputs/hbcc/HBCC_InhibNeu.rds
export HBCC_RDS_MIC=inputs/hbcc/HBCC_Mic.rds
export HBCC_RDS_OL=inputs/hbcc/HBCC_OL.rds
export HBCC_GLOBAL_UMAP=inputs/hbcc/HBCC_X_umap.csv.gz
export HBCC_LIGHTWEIGHT_OUTPUT=inputs/hbcc/hbcc_lightweight.rds
export KEEP_GENES=AIF1,GFAP
export HBCC_RESULT_ROOT=outputs/snRNA/hbcc/results/hbcc_IL6_<RUN_ID>

Rscript codes/04_snRNAseq/human_hbcc/00_build_lightweight_object.R

export SEURAT_RDS_PATH=$HBCC_RDS_AST
Rscript codes/04_snRNAseq/human_hbcc/04_aucell_scoring.R
export SEURAT_RDS_PATH=$HBCC_RDS_EXCITNEU
Rscript codes/04_snRNAseq/human_hbcc/04_aucell_scoring.R
export SEURAT_RDS_PATH=$HBCC_RDS_INHIBNEU
Rscript codes/04_snRNAseq/human_hbcc/04_aucell_scoring.R
export SEURAT_RDS_PATH=$HBCC_RDS_MIC
Rscript codes/04_snRNAseq/human_hbcc/04_aucell_scoring.R
export SEURAT_RDS_PATH=$HBCC_RDS_OL
Rscript codes/04_snRNAseq/human_hbcc/04_aucell_scoring.R

export SEURAT_RDS_PATH=inputs/hbcc/hbcc_lightweight.rds
export HBCC_METADATA_PATH=source_data/snRNAseq/hbcc_donor_filter_manifest.csv
export HBCC_MIN_AGE_YEARS=21
export HBCC_MIN_CELLS_PER_CT=50
export HBCC_EXCLUDE_COMORBID=TRUE
Rscript codes/04_snRNAseq/human_hbcc/05_postprocess_figures_statistics.R
```

Step 4 is run separately on each full-feature cell-class object so that AUCell
rankings retain the original transcriptome-wide feature universe. The
lightweight merged object is constructed only for step 5 visualization and
metadata integration. All five step 4 runs must use the same
`HBCC_RESULT_ROOT`.

The release manifest stores the finalized donor inclusion decisions. Donor
filtering retains 84 controls and 50 schizophrenia donors before cell
availability is evaluated. One retained control donor has no nuclei in the
figure-level tables, which therefore contain 83 control and 50 schizophrenia
donors. The same script also accepts the original HBCC metadata schema for
compatibility.

### Mouse mPFC workflow

`MOUSE_INPUT_DIR` contains one subdirectory per sequencing library, and each
subdirectory holds its 10x matrix under `filtered_cell_gene_matrix/`. The script
lists these subdirectories with `list.dirs()` and strips a leading `sampleN_`
from each name; the resulting names must be `con1`, `con2`, `con3`, `ko1`, `ko2`
and `ko3`, which are the keys of the per-sample quality-control threshold table.

```bash
export MOUSE_INPUT_DIR=inputs/mouse_mPFC
export MOUSE_OUTPUT_DIR=outputs/snRNA/mouse_mPFC

Rscript codes/04_snRNAseq/mouse_mPFC/01_qc.R
Rscript codes/04_snRNAseq/mouse_mPFC/02_harmony.R
Rscript codes/04_snRNAseq/mouse_mPFC/03_annotation.R

export SEURAT_RDS_PATH=outputs/snRNA/mouse_mPFC/2-harmony_ALL/res_0.8/dim_35_res_0.8_harmony_CelltypesV3.rds
Rscript codes/04_snRNAseq/mouse_mPFC/04_aucell_scoring.R

export MOUSE_RESULT_ROOT=outputs/snRNA/mouse_mPFC/2-harmony_ALL/res_0.8/dim_35/results/mouse_IL6_<TIMESTAMP>
Rscript codes/04_snRNAseq/mouse_mPFC/05_postprocess_figures_statistics.R
```

The manifest `source_data/snRNAseq/mouse_mPFC_snRNAseq_sample_manifest.csv` records the
six sequencing libraries, that is three biological replicates per group, each prepared from
mPFC tissue pooled from two mice, together with group, sex, the number of mice pooled per
replicate and analyzed-nucleus counts.

## Figure mapping

| Manuscript output                                | Script                                                                                 | Main released input                                                     |
| ------------------------------------------------ | -------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Supplementary Figure 1a-c                        | `01_bulk_transcriptomics/01a_cohort_effect_sizes.R` (controlled-access inputs)         | PsychENCODE count matrices or preprocessed donor-level tables           |
| Figure 1b                                        | `01_bulk_transcriptomics/01b_meta_analysis.R`                                          | `bulk_DLPFC_meta_analysis_input.csv`                                    |
| Figure 1c-e; Supplementary Figure 1d-g           | `04_snRNAseq/human_hbcc/04_aucell_scoring.R` and `05_postprocess_figures_statistics.R` | HBCC H5AD-derived objects and `hbcc_donor_filter_manifest.csv`          |
| Figure 1f-i; Supplementary Figure 1h-k           | `04_snRNAseq/human_batiuk/01_qc.R` through `05_postprocess_figures_statistics.R`       | Zenodo 6921620 and `batiuk_sample_manifest.csv`                         |
| Supplementary Figure 16                          | `06_CAMK2A_AUCell_validation/06_CAMK2A_AUCell_validation.R`                            | Final Batiuk object and step 4 tables                                   |
| Supplementary Figure 4a                          | `02_human_blood_expression/cohorts_I_IV_volcano_plots.py`                              | `cohorts_I_IV_miRNA_volcano_plot_data.xlsx`                             |
| Figure 2f                                        | `02_human_blood_expression/cohort_VII_IL6_let7b_correlation.py`                        | `cohort_VII_IL6_let7b_expression_data.xlsx`                             |
| Figure 2i                                        | `03_mouse_expression/mouse_mPFC_Il6_let7b_correlation.py`                              | `mouse_mPFC_Il6_let7b_qPCR_data.tsv`                                    |
| Figure 3e-h; Figure 4a-b; Supplementary Figure 8 | `04_snRNAseq/mouse_mPFC/01_qc.R` through `05_postprocess_figures_statistics.R`         | Deposited mouse snRNA-seq and `mouse_mPFC_snRNAseq_sample_manifest.csv` |
| Figure 3k                                        | `05_IISS/iiss_umap.py`                                                                 | Four IISS assigned-count matrices                                       |
| Figure 3l-m                                      | `05_IISS/iiss_celltype_summary_and_figure3.R`                                          | Four IISS assigned-count matrices                                       |
| Figure 3n; Supplementary Figure 9a               | `05_IISS/iiss_IL6_IL6ST_heatmaps.R`                                                    | `IISS_IL6_IL6ST_cell_type_summary.xlsx`                                 |

The numbered snRNA-seq step 5 scripts also generate figure-to-Source-Data
manifests within their output directories.

The table above lists the panels that are generated by scripts in this archive.
The remaining panels fall into two groups. Behavioural, imaging, qPCR, ELISA,
immunoblot and photometry measurements were plotted in GraphPad Prism from the
values reported in the Source Data workbook and Supplementary Table 8, and
involve no custom analysis code. The IISS tissue and spatial map panels were
rendered from the processed spatial files now released under
`source_data/IISS/processed_spatial/`; the public validator checks their exact
count correspondence to the analysis matrices, but the commercial primary
image-decoding workflow is not reimplemented. Figure 2d and 2e and
Supplementary Figure 4b and 4c belong to the first group, and Figure 2d and 2e
use the same released workbook,
`cohort_VII_IL6_let7b_expression_data.xlsx`, as Figure 2f. The Source Data
workbook accompanying the manuscript provides the underlying values for every
plotted panel. Per-analysis-unit proportion tables are released for the cell
types that have a dedicated per-unit panel; for the remaining forest-plot rows
the point estimate, confidence interval, P value and q value are given in the
forest-plot source data.

## Validation

Stage 2 of the meta-analysis, both correlation analyses, all three IISS figure
modules and the processed-spatial validator were rerun from the released inputs
in this archive. Stage 1 of the meta-analysis was tested in both input modes
with synthetic count and donor-level fixtures and checked against the original
analysis code. A full cohort-level Stage 1 rerun requires controlled-access
PsychENCODE inputs and is therefore not claimed here.

Rerunning the three IISS modules regenerates `Fig3k_IISS_UMAP_source_data.csv`,
the two Figure 3l composition tables and `Fig3m_IISS_marker_dotplot_source_data.csv`
byte-for-byte, including the UMAP embedding, which is deterministic under the
fixed seed. The spatial validator confirms exact aggregation of the four
`CellGene.csv` tables to all retained cell-by-gene values in the released
matrices and verifies one-to-one centroid-to-boundary identifiers. Stage 2 of
the meta-analysis reproduces the pooled estimate, its
confidence limits, the heterogeneity statistics and the per-cohort significance
marks. The two correlation modules reproduce the reported coefficients, the
Shapiro-Wilk-guided choice of test and the reported P values. All Python files
pass syntax checks and all R files parse successfully.

The large snRNA-seq scripts have been checked for syntax, parameter interfaces,
fixed seeds, thresholds and factor order. The retained `sessionInfo()` files
document the runs that produced the reported results. Full
end-to-end reruns require the public HBCC and Batiuk inputs and the deposited
mouse dataset; this archive does not claim a clean-system large-object rerun
where those inputs were unavailable locally.

`RELEASE_MANIFEST.tsv` lists each released file, its byte size and SHA-256
digest. `SHA256SUMS.tsv` provides the same digests in checksum-first form. The
two manifest files exclude themselves to avoid recursive checksums.

## Statistical conventions

- Randomized operations use fixed seeds.
- Human binary snRNA-seq outcomes account for within-subject dependence using
  donor or released-sample identity as appropriate. The HBCC tables report
  donor counts as `n_donors_*`; the Batiuk tables, whose analysis units are
  released sample matrices, report `n_samples_CTR`, `n_samples_SCZ` and
  `n_samples_per_group`.
- Mouse snRNA-seq outcomes account for biological-replicate identity. Each replicate is
  mPFC tissue pooled from two mice and is the experimental unit in all models. The legacy
  column names `n_animals_per_group`, `n_animals_syn_control` and
  `n_animals_syn_let7b_ko` denote these replicates.
- IISS heatmaps and fold changes are descriptive sample-level summaries.
- Correlation analyses follow the prespecified Shapiro-Wilk-guided
  Pearson-or-Spearman rule.
- Figure-level sample sizes and tests are reported in Supplementary Table 8.

## Citation

Cite the associated manuscript, the original datasets and the software used by
each workflow.
