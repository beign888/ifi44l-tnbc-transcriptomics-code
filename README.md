# IFI44L in triple-negative breast cancer: reproducibility code

Version: `v1.1.0` public-multiomics manuscript code snapshot (8 November 2026)

This repository contains the analysis code used for a compartment-aware study of IFI44L in triple-negative breast cancer (TNBC). It covers receptor-defined TCGA-BRCA bulk transcriptomics, TCGA overall-survival (OS) and progression-free-interval (PFI) analyses, GSE58812 survival validation, GSE176078 single-cell RNA sequencing and GSE210616 10x Visium spatial transcriptomics.

The repository intentionally contains no author list, contribution statement, local patient material, immunohistochemistry, identifiable information, third-party raw data, cell-line measurements, immunoblots or manuscript source. Public inputs must be retrieved from the original repositories.

## Scientific boundary

- The inferential unit is the tumour, patient or donor—not an individual cell or Visium spot.
- IFI44L is excluded from every expression signature tested for correlation with IFI44L.
- Visium spots are treated as multicellular spatial niches, not single cells.
- Continuous standardized IFI44L is the primary survival predictor; median groups are used only for graphical sensitivity analyses.
- TCGA OS/PFI results are retained regardless of significance. The directionally favourable but non-significant TCGA estimates do not establish an independent prognostic biomarker.
- The conceptual synthesis is descriptive and does not assert a causal pathway or immunotherapy-predictive utility.

## Repository layout

```text
.
├── scripts/                 analysis, figure and validation scripts
├── environment/             recorded Python/R package versions
├── PUBLIC_DATA_ACCESS.md    accessions and input staging notes
├── SOURCE_DATA_MAP.tsv      figure-to-source-to-script map
├── FILE_MANIFEST_SHA256.tsv repository integrity manifest
├── LICENSE                  MIT License for the repository code
└── LICENSE_STATUS.md        licence scope and third-party exclusions
```

Runtime-created directories (`source_metadata/`, `results/`, `figures/` and `logs/`) are excluded from version control.

## Software environment

Python 3.11 or later and R 4.3 are required. The recorded package versions are in `environment/`.

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
    ├── TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx
    ├── GSE58812/GSE58812_series_matrix.txt.gz
    ├── Wu2021_NatGenet_Supplementary_Tables.xlsx
    ├── GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad
    └── spatial/GSE210616/
        ├── GSE210616_family.soft.gz
        ├── GSE210616_primary28_samples.tsv
        └── files/GSM*/...
```

The three TCGA expression/receptor helper files are frozen inputs used by the submitted analysis and are not included here. They reconstruct an archived ER-negative, PR-negative and HER2-negative intersection; they are not a PAM50 basal-like surrogate. `21_audit_tcga_tnbc_definition.py` checks the exact intersection and expression-matrix match before modelling.

Download the TCGA Clinical Data Resource from the official GDC file endpoint:

```bash
mkdir -p source_metadata/TCGA_CDR
curl --fail --location \
  --output source_metadata/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx \
  https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
```

Prepare the 28 primary GSE210616 sections from the GEO SOFT record, optionally downloading their files:

```bash
python scripts/08_prepare_gse210616_spatial.py \
  --soft source_metadata/spatial/GSE210616/GSE210616_family.soft.gz \
  --output-dir source_metadata/spatial/GSE210616 \
  --download
```

## Reproduction commands

Run the non-destructive preflight first:

```bash
python scripts/17_reproducibility_preflight.py \
  --project "$PWD" \
  --raw-root "/absolute/path/to/RAW_ROOT" \
  --deep-hash
```

Run all public-data analyses and Figures 1–5:

```bash
python scripts/18_run_nonexperimental_pipeline.py \
  --project "$PWD" \
  --raw-root "/absolute/path/to/RAW_ROOT" \
  --stage all \
  --permutations 499
```

Add `--run-census-audit` only when network access to the frozen CELLxGENE Census version is available. Logs are written to `logs/reproduction/`. The pipeline stops at the first failed step and finishes by rerunning the preflight against expected outputs.

## Validation status

The v1.1.0 snapshot was checked for Python syntax, R parseability, required-input declarations, forbidden local paths, author/contact information, secrets and internal-manifest consistency. The TCGA OS/PFI script and revised public-multiomics Figures 1, 2 and 5 were executed against the preserved public-data workspace. These checks establish a bounded local validation, not an independent numerical rerun from newly downloaded inputs or scientific peer review. See `VALIDATION_REPORT.md`.

## Citation and authorship

Citation metadata and author attribution are intentionally deferred until the final authorship record has been approved. The GitHub account or commit identity hosting this snapshot is not, by itself, evidence of scholarly authorship or contribution.

## Licence

The repository code is released under the MIT License; see `LICENSE`. The licence does not cover third-party datasets, publications or software dependencies, which remain subject to their own terms. See `LICENSE_STATUS.md` for the scope.
