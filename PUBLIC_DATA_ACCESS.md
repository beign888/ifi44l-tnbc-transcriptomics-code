# Public data access and redistribution boundary

| Resource | Identifier | Role |
|---|---|---|
| GDC / TCGA | TCGA-BRCA | Receptor-defined bulk TNBC discovery cohort |
| GDC publication supplement | `1b5f413e-a8d1-4d10-92eb-7c4ae739ed81` | Curated TCGA OS/PFI endpoints |
| GEO | GSE58812 | Independent bulk and survival validation |
| GEO | GSE176078 | Single-cell compartment localization |
| GEO | GSE210616 | Visium spatial analysis |

Repository terms apply. Retrieve raw data directly from the source repositories.

## Required frozen inputs

The analysis expects the exact filenames and directory layout documented in `README.md`. Raw public matrices, histology images and third-party supplementary workbooks are not redistributed in this repository. The TCGA expression matrix, archived receptor-status source table and receptor-defined TNBC sample list are likewise not included in this code-only snapshot.

The standardized TCGA clinical outcomes are from Liu et al., *Cell* 2018 (doi: `10.1016/j.cell.2018.02.052`) and the official NCI GDC publication page. Download `TCGA-CDR-SupplementalTableS1.xlsx` from `https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81`.

GSE210616 preparation begins from the GEO family SOFT record and is performed by `scripts/08_prepare_gse210616_spatial.py`. The CELLxGENE metadata audit is version-pinned to Census `2025-11-08` and dataset ID `9fddb063-056d-4202-8b8a-4b0ee531d3ce`.

The repository contains no locally collected patient cohort, immunohistochemistry, tissue images, cell-line measurements or identifiable participant information.
