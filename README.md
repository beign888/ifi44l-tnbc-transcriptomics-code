# IFI44L in triple-negative breast cancer: reproducibility code

Version: `v1.0.0` manuscript code snapshot (2 November 2026)

This repository contains the analysis code used for a compartment-aware study of IFI44L in triple-negative breast cancer (TNBC). It covers receptor-defined TCGA-BRCA bulk transcriptomics, GSE58812 survival validation, GSE176078 single-cell RNA sequencing, GSE210616 10x Visium spatial transcriptomics and an explicitly exploratory GSE169246 treatment-cohort analysis.

The repository intentionally contains no author list, contribution statement, local patient material, immunohistochemistry, identifiable information, third-party raw data, cell-line experimental measurements or manuscript source. Public data must be retrieved from the original repositories. Cell-line source data and uncropped blots are supplied separately with the journal submission, not here.

## Scientific boundary

- The inferential unit is the tumour, patient, donor or independent biological experiment, as appropriate—not an individual cell, Visium spot, technical well or image field.
- IFI44L is excluded from every expression signature tested for correlation with IFI44L.
- Visium spots are treated as multicellular spatial niches, not single cells.
- GSE169246 is an underpowered exploratory analysis and is not evidence for or against an immunotherapy-predictive biomarker.
- Archived western blots are qualitative siRNA-screening documentation only; no cell-line-specific or quantitative immunoblot inference is made.

## Repository layout

```text
.
├── scripts/                 analysis, figure and validation scripts
├── environment/             recorded Python/R package versions
├── PUBLIC_DATA_ACCESS.md    accessions and input staging notes
├── SOURCE_DATA_MAP.tsv      figure-to-source-to-script map
├── FILE_MANIFEST_SHA256.tsv repository integrity manifest
└── LICENSE_STATUS.md        current reuse-rights status
```

Runtime-created directories (`source_metadata/`, `results/`, `figures/`, `logs/` and `bin/`) are excluded from version control.

## Software environment

Python 3.11 or later, R 4.3, a C compiler and zlib are required. The recorded package versions are in `environment/`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/python-requirements.txt
Rscript environment/install_R_packages.R
```

Exact historical versions may not remain installable on every operating system. If a compatible environment is used instead, record the resolved versions in the preflight report.

## Input staging

Obtain all source files from the repositories listed in `PUBLIC_DATA_ACCESS.md`. Stage them as follows:

```text
RAW_ROOT/
└── 数据/差异富集分析/
    ├── 1 乳腺癌分型/3_乳腺癌分型.xlsx
    └── 3 差异分析/
        ├── symbol.BRCA.txt
        └── 三阴.txt

REPOSITORY_ROOT/
└── source_metadata/
    ├── GSE58812/GSE58812_series_matrix.txt.gz
    ├── Wu2021_NatGenet_Supplementary_Tables.xlsx
    ├── GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad
    ├── spatial/GSE210616/
    │   ├── GSE210616_family.soft.gz
    │   ├── GSE210616_primary28_samples.tsv
    │   └── files/GSM*/...
    └── single_cell/GSE169246/
        ├── GSE169246_TNBC_RNA.barcode.tsv.gz
        ├── GSE169246_TNBC_RNA.feature.tsv.gz
        ├── GSE169246_TNBC_RNA.counts.mtx.gz
        └── Zhang2021_CancerCell_mmc2.xlsx
```

The three TCGA helper files are frozen inputs used by the submitted analysis. They are not included in this code-only repository. The sample-level receptor-definition audit is performed by `21_audit_tcga_tnbc_definition.py` before modelling.

Prepare the 28 primary GSE210616 sections from the GEO SOFT record, optionally downloading their files:

```bash
python scripts/08_prepare_gse210616_spatial.py \
  --soft source_metadata/spatial/GSE210616/GSE210616_family.soft.gz \
  --output-dir source_metadata/spatial/GSE210616 \
  --download
```

## Reproduction commands

Run the preflight first. It reports missing inputs and packages without fabricating replacements.

```bash
python scripts/17_reproducibility_preflight.py \
  --project "$PWD" \
  --raw-root "/absolute/path/to/RAW_ROOT" \
  --deep-hash
```

Run the public-data analyses and main figures:

```bash
python scripts/18_run_nonexperimental_pipeline.py \
  --project "$PWD" \
  --raw-root "/absolute/path/to/RAW_ROOT" \
  --stage all \
  --permutations 499
```

Add `--run-census-audit` only when network access to the frozen CELLxGENE Census version is available. Logs are written to `logs/reproduction/`. The pipeline stops at the first failed step and finishes by rerunning the preflight against the expected outputs.

The optional `36_make_supplementary_figures_v8.py` script is included for auditability but is not part of the public-data pipeline because its cell-line measurements and image assets are supplied separately with the manuscript.

## Validation status

The public snapshot was checked for Python syntax, R parseability, C compilation, forbidden local paths, author/contact information, secrets and internal-manifest consistency. These checks establish code-package integrity only. They do not substitute for a full numerical rerun from independently downloaded raw data or for scientific peer review. See `VALIDATION_REPORT.md`.

## Citation and authorship

Citation metadata and author attribution are intentionally deferred until the final authorship record has been approved. The GitHub account or commit identity hosting this snapshot is not, by itself, evidence of scholarly authorship or contribution.

## Licence

No reuse licence has yet been assigned. See `LICENSE_STATUS.md` before copying, modifying or redistributing the code.
