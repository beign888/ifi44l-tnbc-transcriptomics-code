# Public sources and reproducibility inputs

Version `v1.4.0`. Source preparation, full/subset executions, offline restoration and the six-job extension replay have completed with the bounded checks described in VALIDATION_REPORT.md. The immutable commit and final supporting-archive identity are recorded separately by the release process. Original access/redistribution terms apply. The MIT code licence does not relicense datasets.

| Source | Role | Acquisition |
|---|---|---|
| UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2, 2017-10-13 | Current main TCGA expression | Exact all-gene main subset plus full source header/receipts in Additional file 3; complete source remains public |
| cBioPortal DataHub BRCA clinical table | Outcome-independent ER/PR/HER2 eligibility | Pinned commit/file in Additional file 3 |
| TCGA Clinical Data Resource | OS/PFI | Official GDC publication file, separately downloaded |
| GSE58812 | 107-case bulk/survival cohort | Official GEO series matrix |
| GSE21653 | Additional bulk/DFS cohort | Exact all-probe 87-case subset plus complete 266-sample metadata/header in Additional file 4 |
| GSE176078 / Wu et al. | Single-cell donor/compartment analysis | Original workbook and pinned processed h5ad, separately obtained |
| GSE210616 | 28 Visium sections from 14 patients | Official GEO metadata/section files, separately obtained |
| GSE96058 / SCAN-B | Additional clinical-receptor-defined RNA-seq/OS cohort | Official full expression file and both platform metadata files; frozen exact all-gene subset designated as Additional files 7–8 |

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

## v16 reuse of core inputs

The conditional three-cohort extension reads the saved full-expression gene matrices, source clinical tables and existing 19-gene IFN-alpha scores generated by the core workflow. Cohort-wide standardization precedes endpoint-specific exclusions. The source matrices are already `log2(RSEM normalized count+1)` (TCGA), `log2(MAS5+1)` (GSE58812) and deposited RMA (GSE21653); no further logarithm or cross-platform pooling is applied. The newly scored 16-gene comparison module excludes IFI44L, IFI44, MX1 and ISG15 simultaneously. Source membership/coverage tables and all ten fixed expression comparisons per cohort remain inspectable.

The within-lineage single-cell extension reads `results/single_cell/GSE176078_TNBC_targeted_IFI44L_immune.h5ad` and the corresponding core donor-major summary. Its source IFN-response score is already present in the H5AD; it is not the bulk 19-gene score and is not reconstructed from an incompletely available selected-feature list. The original plan, feature-availability amendment and subsequent exact-permutation amendment retain the distinction. This reuses the existing cohort, rather than supplying a new single-cell validation cohort.

## Additional files 7–8: GSE96058 / SCAN-B frozen inputs

The [official GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) supplies the following expression and metadata files:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL11154_series_matrix.txt.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL18573_series_matrix.txt.gz
```

The source describes transcript FPKM summed by gene symbol and distributed as `log2(FPKM+0.1)`. Use the deposited transformed values directly. Both platform metadata records are required to reconcile all 3,409 expression columns: 3,273 original primary records and 136 explicitly named technical replicates. The fixed policy retains original primary F records and excludes `F*repl` records before expression/outcome association analysis. Source clinical fields `er status=0`, `pgr status=0` and `her2 status=0` define the TNBC subset; RNA-predicted MGC/SGC receptor calls and PAM50 do not replace them. Historical binary diagnostic calls do not justify reconstructing a common modern numerical receptor cutoff.

The frozen preparation design contains 143 eligible clinical records, with 143 positive-time OS records/26 events and 137 complete pathological records/26 events. These source-clinical feasibility counts are not evidence of an IFI44L association. Actual source-bound execution retains platform-specific linkage and case/event counts: HiSeq 2000 has 126 cases/25 OS events and NextSeq 500 has 17 cases/one event. The source-platform amendment was frozen before association fitting and accompanies the original plan. It retains all original specifications and adds an instrument-stratified nested pair, yielding nine Cox models and three fixed LRTs with a single BH3 family. Source metadata identify 126 HiSeq 2000 and 17 NextSeq 500 cases.

After staging the locked plans, the helper can acquire the official full sources and preserve the exact selected decimal strings without an additional transform:

```sh
python scripts/34_prepare_GSE96058_inputs.py --project "$PWD" --download \
  --analysis-plan "$PWD/protocols/v16/GSE96058_ANALYSIS_PLAN.md" \
  --analysis-plan-amendment "$PWD/protocols/v16/GSE96058_PLATFORM_AMENDMENT.md"
```

Its staged directory is `source_metadata/GSE96058/`. The designated frozen input bundle includes both complete official metadata files; `GSE96058_all3409_metadata.tsv`; the selected metadata; replicate pairs and full selection/exclusion evidence; `GSE96058_complete_source_header.tsv`; sample-to-expression linkage; `GSE96058_TNBC143_frozen_all_genes.tsv.gz`; fixed module definitions; processing text; and `GSE96058_frozen_input_provenance.json`. The completed source preparation verifies 30,865 gene rows and source-column linkage, reading the full gzip to EOF. Full-source, exact-subset, metadata and linkage hashes are recorded in the acquisition/provenance receipt and must be checked against the staged bytes; uncompleted transfers must not be treated as frozen input. The complete large expression file need not be redistributed when the verified exact all-gene subset is staged, but its public URL, full header and source hash must remain available.

The full official gene-expression gzip was actually acquired: 591,676,211 bytes, SHA-256 `3c717baf7960e1f1477a72744a399a42fabba113c9e7ecb3d55997574d7e9732`. Streaming comparison verified all 4,413,695 selected expression decimal strings and 30,865 gene identifiers. Final full-source and subset R routes matched all 20 compared analysis TSVs and the normalized gene-expression gzip byte-for-byte. This checks source preservation and same-script reproducibility, not independent biological replication.

The exact gzip is 25,592,924 bytes (SHA-256 `2e49a251c2f9129658a2924e78372ebaf4fb936466e1553d7a7e1e6410e5b5ba`), exceeding the per-file supporting-material limit. Its original bytes are split without decompression, rounding, filtering or re-encoding: Additional file 7 contains metadata/provenance plus `.part01` (16,777,216 bytes), and Additional file 8 contains `.part02` (8,815,708 bytes). Both ZIPs are required and must be extracted into the same `source_metadata/` directory. `GSE96058_INPUT_PARTS.json` identifies their order, individual hashes and complete-file hash. Part01 SHA-256 is `48b6efe979214168ae82d3de62c1ee3ebb763a5086b55adca3ff04e0cff2c6a0`; part02 SHA-256 is `0f352779992468d7a5b359451f6310473cf714faf9bdc7e057617e1b1988418c`. These are binary transport volumes, not separate matrices or independently decompressible gzip files.

```sh
mkdir -p source_metadata
unzip -n /path/to/Additional_file_7.zip -d source_metadata
unzip -n /path/to/Additional_file_8.zip -d source_metadata
python scripts/34_prepare_GSE96058_inputs.py --project "$PWD" --restore-parts
```

`--restore-parts` is an offline concatenation-and-SHA-verification operation. It rejects a size/hash mismatch, retains both parts and never overwrites a different complete source. It must not be combined with `--download`. The final extension wrapper actually completed offline restoration before fitting; a separate isolated restoration test also reproduced the complete size and hash. The explicit command above is available for reuse. No source file is deleted. The large original 3,409-column source can remain outside the supporting archive while its full hash/header/metadata are preserved.

Additional file 6 is designated for supplementary methods, Figures S2–S4 and complete extension result/plotting-source tables. It is distinct from the two frozen-input transport archives. The final package therefore has eight Additional files. Final ZIP CRCs, per-file manifests and release synchronization remain separate from completed source preparation and must be recorded before those archives are described as sealed. Original data terms remain applicable.
