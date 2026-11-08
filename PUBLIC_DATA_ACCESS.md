# Public data access and provenance

Snapshot: `v1.2.0`, prepared 9 November 2026. Original access and redistribution terms apply; the MIT code licence does not relicense datasets.

| Resource | Identifier and source | Role |
|---|---|---|
| TCGA-derived archived expression | Additional files 3 and 4; exact frozen 151-column subset | Main bulk cohort; strict 110-patient subset |
| cBioPortal DataHub clinical table | Commit `4c975725c54df7e31e0610fad22a742ad01cdfa3`, `public/brca_tcga/data_clinical_patient.txt` | ER/PR/HER2 reconciliation independent of outcomes |
| TCGA Clinical Data Resource | [Official GDC publication page](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018) | Curated OS/PFI endpoints |
| GEO | [GSE58812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58812) | Independent bulk/survival analysis |
| GEO / CELLxGENE | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Single-cell donor/compartment analysis |
| GEO | [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) | Visium analysis, 28 sections from 14 patients |

## Corrected TCGA eligibility

`scripts/21_audit_tcga_tnbc_definition.py` matches the original 1,072-sample header to the pinned clinical table. Main eligibility requires explicit negative ER and PR, negative HER2 IHC or FISH evidence, and exclusion of any positive HER2 IHC/FISH or IHC score 3. This yields 151 unique patients. Strict sensitivity additionally requires explicit HER2 IHC negativity and yields 110. Neither uses PAM50 for inclusion. The old 187-sample helper list is superseded: intersecting its labels did not validate the underlying receptor statuses.

Pinned clinical source:

```text
https://raw.githubusercontent.com/cBioPortal/datahub/4c975725c54df7e31e0610fad22a742ad01cdfa3/public/brca_tcga/data_clinical_patient.txt
SHA-256: 973b6e03436330e0563a15989d04fb3566a2c7cdc66db9aaa3badbb32c58badc
Staged path: source_metadata/TCGA_clinical/brca_tcga_clinical_patient.txt
```

The file is included in Additional file 3. The audit script rejects a different SHA rather than silently using a newer table. This public clinical source differs from the TCGA-CDR outcome table and does not recover expression-processing history.

## Frozen expression: Additional files 3 and 4

Extract both ZIPs into `source_metadata/`. The first contains `TCGA_frozen/TCGA_BRCA_TNBC_frozen.tsv.gz.part01`, the original header, provenance JSON and pinned clinical source; the second contains `.part02`. Each binary part is 13,586,530 bytes. Run:

```bash
python scripts/00_restore_frozen_tcga_input.py --project "$PWD"
```

The script concatenates both parts and checks:

```text
Restored file: source_metadata/TCGA_frozen/TCGA_BRCA_TNBC_frozen.tsv.gz
Size: 27173060 bytes
SHA-256: 2af76ae44807e00e7bfe269655ba43cc7162f0b63ed0c1e036fc1d502cfe6904
Selected sample columns: 151
Gene rows before analysis filters: 56536
Original full-matrix SHA-256: f53fc5c55002b22472aca402ecc6a52946a79b91ec91db658288666e9939ccdd
```

Original decimal strings are preserved without rounding or transformation. The archived records do **not** establish upstream expression release, normalization or measurement units. Do not label these values TPM/FPKM/RSEM or interpret coefficients as log2 fold changes. Reproducing calculations on this frozen subset does not recover upstream processing. The corrected median-group coefficient is an archived-unit expression difference, with no fold-change threshold asserted.

## Other required inputs

TCGA outcomes are from Liu et al., *Cell* 2018 (doi: `10.1016/j.cell.2018.02.052`). Stage the official file:

```bash
mkdir -p source_metadata/TCGA_CDR
curl --fail --location \
  --output source_metadata/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx \
  https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
```

The survival script checks MD5 `a4591b2dcee39591f59e5e25a6ce75fa` before joining outcomes. Retain both source-paper citations and dataset accession/URL references where applicable.

Stage the official GSE58812 series-matrix gzip at `source_metadata/GSE58812/GSE58812_series_matrix.txt.gz`. Stage Wu et al.'s workbook at `source_metadata/Wu2021_NatGenet_Supplementary_Tables.xlsx` and the matching processed h5ad at `source_metadata/GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad`. The replay h5ad SHA-256 is `04245d73cdcc0b473aa88ec655d1bcb051dbbd5f9a3d1aec4c2f37b1f0645df4`. The optional metadata audit is pinned to Census `2025-11-08`, dataset `9fddb063-056d-4202-8b8a-4b0ee531d3ce`; it does not regenerate the h5ad from raw sequencing.

GSE210616 preparation starts from official `source_metadata/spatial/GSE210616/GSE210616_family.soft.gz`. `scripts/08_prepare_gse210616_spatial.py` creates `GSE210616_primary28_samples.tsv` and downloads section matrices, positions, images and scalefactors. Preserve exact source versions/checksums. The code repository itself does not bundle these third-party raw matrices, images or workbooks.

The [repository homepage](https://github.com/beign888/ifi44l-tnbc-transcriptomics-code) is linked until a corrected release commit is verified. The final package manifest and manuscript must identify the same exact commit; v1.1.0 does not identify v1.2.0.
