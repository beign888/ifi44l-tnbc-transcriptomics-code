# Public sources and reproducibility inputs

Version `v1.3.0`; The immutable commit is recorded in the manuscript and separate release receipt. Original access/redistribution terms apply. The MIT code licence does not relicense datasets.

| Source | Role | Acquisition |
|---|---|---|
| UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2, 2017-10-13 | Current main TCGA expression | Exact all-gene main subset plus full source header/receipts in Additional file 3; complete source remains public |
| cBioPortal DataHub BRCA clinical table | Outcome-independent ER/PR/HER2 eligibility | Pinned commit/file in Additional file 3 |
| TCGA Clinical Data Resource | OS/PFI | Official GDC publication file, separately downloaded |
| GSE58812 | 107-case bulk/survival cohort | Official GEO series matrix |
| GSE21653 | Additional bulk/DFS cohort | Exact all-probe 87-case subset plus complete 266-sample metadata/header in Additional file 4 |
| GSE176078 / Wu et al. | Single-cell donor/compartment analysis | Original workbook and pinned processed h5ad, separately obtained |
| GSE210616 | 28 Visium sections from 14 patients | Official GEO metadata/section files, separately obtained |

## Additional file 3: documented TCGA primary input

The [matrix](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz) and [producer/processing metadata](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.json) document UNC Level-3 RSEM normalized results, distributed by Xena as `log2(norm_count+1)`, version `2017-10-13`. Use the supplied values directly; no further log or count normalization is performed. Metadata records public access and redistribution permission, subject to original terms.

```text
Complete source: TCGA_BRCA_HiSeqV2_RSEM_log2normcount_2017-10-13.tsv.gz
Source size: 64286337 bytes; 20530 gene rows; 1218 sample columns
Source SHA-256: 263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746
Exact supplied subset: TCGA_BRCA_documented_primary_subset.tsv.gz
Subset size: 8632958 bytes; 20530 rows; 157 columns
Subset SHA-256: e268fa9d9da652c6f7418b408f663fbe60534e62ab98cd3791a51d593b402999
Full-header SHA-256: 2554effaed790acec265d15d2fd2b76647ff26923cc5b0d83945e5cc98493d2f
Staged directory: source_metadata/TCGA_documented_primary
```

The complete source header defines the selection universe even when the exact subset is read. Type-01 primary samples are reconciled to the pinned clinical source; the lexicographically first primary barcode is chosen when a patient has more than one. Main selection requires explicit ER/PR IHC negativity, resolved negative HER2 IHC/FISH evidence and no positive HER2 evidence/IHC3+ result (157 patients); strict selection additionally requires HER2-IHC negativity (113). No PAM50 or outcome-based inclusion is used.

Pinned receptor source:

```text
https://media.githubusercontent.com/media/cBioPortal/datahub/4c975725c54df7e31e0610fad22a742ad01cdfa3/public/brca_tcga/data_clinical_patient.txt
SHA-256: 973b6e03436330e0563a15989d04fb3566a2c7cdc66db9aaa3badbb32c58badc
Staged path: source_metadata/TCGA_clinical/brca_tcga_clinical_patient.txt
```

The corresponding GitHub raw URL serves a Git-LFS pointer, not the clinical table bytes. Use the media URL above or the supplied verified clinical file. Extract Additional file 3 into `source_metadata/`; no binary-part concatenation is needed. It includes the complete header, metadata, acquisition receipt and `TCGA_documented_input_provenance.json`. The scripts reject other input hashes and do not fall back to an undocumented matrix. The 64 MB complete source is unnecessary for the supplied-subset route.

To prepare a bundle again from the full source, retain the provided acquisition receipt/metadata and run:

```sh
python scripts/24_prepare_documented_tcga_input.py --project "$PWD" \
  --expression-matrix /path/to/TCGA_BRCA_HiSeqV2_RSEM_log2normcount_2017-10-13.tsv.gz \
  --source-receipt source_metadata/TCGA_documented_primary/REACQUISITION_RECEIPT.json \
  --metadata-file source_metadata/TCGA_documented_primary/HiSeqV2.json \
  --clinical-file source_metadata/TCGA_clinical/brca_tcga_clinical_patient.txt
```

The legacy annotation lacks TRAC: CD8 T-cell coverage is 7/8; IFN-alpha response retains 19/19. Per-module coverage is output explicitly; module scores are not cell fractions.

## Additional file 4: GSE21653

[Official GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21653) provides 266 GPL570 samples. [The complete processed matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21653/matrix/GSE21653_series_matrix.txt.gz) contains 54,613 probes × 266 columns and reports RMA with quantile normalization. No second log transformation is used. Historical clinical receptor calls `er ihc=0`, `pr ihc=0`, `erbb2=0` select 87 cases, with no molecular-subtype substitution.

```text
Complete matrix size: 22510462 bytes
Complete SHA-256: f13e3bf972bd92a6fd53a1fbe42ad7c84caf4651df3a5e836f70dd3b22b5e76f
Supplied subset: GSE21653_TNBC87_frozen_all_probes.tsv.gz
Subset size: 7380481 bytes; 54613 probe rows; 87 columns
Subset SHA-256: 66cecd1e1fc29bccf1eec1cfa4631468441d04c1d812e1aa846cf4bcdf6c592a
Staged directory: source_metadata/GSE21653
```

The supplied directory also contains the complete 266-case metadata/header, selected list, probe annotation, fixed module definitions, design and source/provenance audits. All probe rows and source decimal strings are preserved. `27_prepare_gse21653_inputs.py --project "$PWD" --download` fetches the exact full matrix and official sample metadata, verifies source hashes and reconstructs the subset. It requires the supplied module-definition file in the staged directory. The R analysis accepts the full or supplied subset and checks their identities.

The fixed IFI44L probe is `204439_at` (Entrez10964); hgu133plus2.db 3.13.0 annotation and limma::avereps are used consistently with GSE58812. DFS was defined by the original study from diagnosis to first recurrence, or death when recurrence information was unavailable. The primary positive-time analysis is 83 patients/25 events; two unchanged zero-time events are included in a separate 85-patient/27-event sensitivity. Age/pT/pN/grade complete cases number 81/24 events, limiting that five-coefficient model to exploratory interpretation. Quantitative historical receptor cutoffs are not reconstructed from the available records.

The original paper is Sabatier et al., *Breast Cancer Research and Treatment* 2011;126:407–420, [doi:10.1007/s10549-010-0897-9](https://doi.org/10.1007/s10549-010-0897-9), PMID20490655. Its [PubMed-linked author manuscript](https://hal.science/hal-00583558/file/PEER_stage2_10.1007%252Fs10549-010-0897-9.pdf) was inspected. It reports IPC/Marseille cases from 1992–2004, whereas GSE58812 reports two other French centres. No GSM accession overlaps, but a public individual-identity crosswalk is unavailable, so this is reported separate-cohort validation rather than identity-verified patient independence.

## Remaining official inputs

Stage the [TCGA-CDR publication file](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018) from Liu et al., *Cell* 2018, doi:10.1016/j.cell.2018.02.052:

```sh
mkdir -p source_metadata/TCGA_CDR
curl --fail --location --output source_metadata/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx \
  https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
```

Its required MD5 is `a4591b2dcee39591f59e5e25a6ce75fa`.

Stage the [GSE58812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58812) matrix at `source_metadata/GSE58812/GSE58812_series_matrix.txt.gz`. Its deposited MAS5 intensities receive `log2(x+1)`; this differs from the already processed GSE21653 source. The inspected source records lack reliably sample-linked stage/grade/node covariates for GSE58812; aggregate summaries are not imputed onto individuals. The annotation audit accompanies Additional file 5.

```sh
mkdir -p source_metadata/GSE58812
curl --fail --location --output source_metadata/GSE58812/GSE58812_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58812/matrix/GSE58812_series_matrix.txt.gz
```

For [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078), stage `Wu2021_NatGenet_Supplementary_Tables.xlsx` and `GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad`. The latter SHA-256 is `04245d73cdcc0b473aa88ec655d1bcb051dbbd5f9a3d1aec4c2f37b1f0645df4`. The optional Census audit uses version `2025-11-08`, dataset `9fddb063-056d-4202-8b8a-4b0ee531d3ce`; it does not regenerate raw sequencing or the processed h5ad.

The installed Census1.17 API provides the versioned source-h5ad retrieval route below. It refuses to overwrite an existing destination. This is an acquisition instruction, not a claim that the large file was reacquired during this documentation step; verify the checksum before analysis.

```sh
python - <<'PY'
import cellxgene_census
cellxgene_census.download_source_h5ad(
    "9fddb063-056d-4202-8b8a-4b0ee531d3ce",
    "source_metadata/GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad",
    census_version="2025-11-08",
)
PY
```

For [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616), stage official `source_metadata/spatial/GSE210616/GSE210616_family.soft.gz` and run `08_prepare_gse210616_spatial.py` to obtain the 28-section manifest, matrices, positions, images and scalefactors. These large public inputs are not all redistributed with the code/supporting files. Retain their original terms and checksums.

```sh
mkdir -p source_metadata/spatial/GSE210616
curl --fail --location --output source_metadata/spatial/GSE210616/GSE210616_family.soft.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210616/soft/GSE210616_family.soft.gz
python scripts/08_prepare_gse210616_spatial.py \
  --soft source_metadata/spatial/GSE210616/GSE210616_family.soft.gz \
  --output-dir source_metadata/spatial/GSE210616 --download
```
