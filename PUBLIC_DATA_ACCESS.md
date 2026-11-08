# Public data access and provenance

Snapshot: `v1.2.1`, prepared 9 November 2026. Original access and redistribution terms apply; the MIT code licence does not relicense datasets.

| Resource | Identifier and source | Role |
|---|---|---|
| TCGA-derived archived expression | Additional files 3 and 4; exact frozen 151-column subset | Main bulk cohort; strict 110-patient subset |
| Documented TCGA expression representation | [UCSC Xena HiSeqV2](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz), version `2017-10-13`; exact 151-column subset in Additional file 5 | Processing sensitivity on the same main-cohort patients; not independent patient validation |
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

## Documented processing sensitivity: Additional file 5

This separate sensitivity uses UCSC Xena dataset `TCGA.BRCA.sampleMap/HiSeqV2`, processed version `2017-10-13`. The [dataset metadata](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.json) identifies the University of North Carolina TCGA genome characterization center as producer, Level-3 RSEM normalized results as the input, and `log2(norm_count+1)` as the distributed units. No further expression transformation was applied. Quantification, annotation and the available gene universe differ jointly from the archived input; this check does not isolate a single normalization effect or establish the archived matrix's upstream provenance.

```text
Source URL: https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz
Complete source file: TCGA_BRCA_HiSeqV2_RSEM_log2normcount_2017-10-13.tsv.gz
Complete source size: 64286337 bytes
Complete source SHA-256: 263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746
Complete source dimensions: 20530 gene rows; 1218 sample columns
Supplied subset: documented_HiSeqV2_matched151_expression.tsv.gz
Subset size: 8251822 bytes
Subset SHA-256: 53d57d8013a86933143964d7ffd68d4247f1bf8737c0c2ce5238e97c2a9a85f8
Subset dimensions: 20530 gene rows; 151 primary-tumour sample columns
Staged directory: source_metadata/TCGA_documented_processing
```

Extract Additional file 5 and copy the contents of `TCGA_processing_sensitivity/` into the staged directory above. The subset preserves source decimal strings exactly, all source gene rows and the frozen clinical-list sample order. All 151 requested sample barcodes match exactly, with no missing, ambiguous or averaged sample. The four `baseline_*` files contain the frozen sample list, original 13 module definitions, archived IFI44L values and archived module scores; they permit standalone comparison without another project directory. Retain `HiSeqV2.json`, the design, acquisition audit, subset provenance and validation JSON alongside the subset.

To run only this sensitivity:

```bash
Rscript scripts/22_documented_tcga_processing_sensitivity.R \
  --work-dir source_metadata/TCGA_documented_processing \
  --output-dir results/processing_sensitivity
```

The default input is the supplied subset. For a complete-source replay, download the file under its recorded filename in the staged directory and add `--input source_metadata/TCGA_documented_processing/TCGA_BRCA_HiSeqV2_RSEM_log2normcount_2017-10-13.tsv.gz` to the command. To reconstruct the subset from this exact source and an already reproduced main 151-patient project, run:

```bash
python scripts/23_prepare_documented_tcga_subset.py \
  --project "$PWD" --work-dir source_metadata/TCGA_documented_processing
```

The helper checks the complete-source SHA before selecting columns and preserves source strings without transformation. It must be used with the main 151-patient baseline results, not the separate strict 110-patient project. The supplied design and source metadata remain required by the analysis. The complete 64 MB input is unnecessary when using the 8 MB supplied subset.

Scoring retains the original 13 module definitions and excludes IFI44L. After invariant-gene filtering, TRAC is unavailable in this legacy source, leaving 7/8 CD8 T-cell genes; all 19 IFN-alpha-response genes are retained. Coverage for every member is output explicitly. No GSEA, survival or adjusted checkpoint model is repeated. The complete-source and subset-only runs both exited zero and generated ten byte-identical result TSVs. `HiSeqV2_source_samples_not_selected.tsv` intentionally differs because the complete source includes unselected patients; it is outside that ten-table identity claim. Independent SciPy correlation/BH verification is recorded in `PROCESSING_SENSITIVITY_VALIDATION.json`.

The Xena metadata records public access and redistribution permission. Its original data terms remain applicable; inclusion as supporting data does not make this expression matrix MIT-licensed.

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

The [repository homepage](https://github.com/beign888/ifi44l-tnbc-transcriptomics-code) is linked until the v1.2.1 release commit is verified. The final package manifest and manuscript must identify the same exact commit; earlier release identifiers do not identify this snapshot.
