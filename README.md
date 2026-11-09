# IFI44L in triple-negative breast cancer: reproducibility code

Version: `v1.3.0`, prepared 9 November 2026. The immutable commit is recorded in the manuscript and separate verified release receipt.

Repository: <https://github.com/beign888/ifi44l-tnbc-transcriptomics-code>. A complete local execution of the final all-stage runner passed on 9 November 2026. All 124 compared analysis TSVs matched the independently staged reference, and all 22 supplied scripts matched the frozen execution snapshot; 17 distinct scripts were invoked. Use the actual remotely verified commit identified by the manuscript/release receipt and verify FILE_MANIFEST_SHA256.tsv.

The primary TCGA analysis uses the documented UCSC Xena `TCGA.BRCA.sampleMap/HiSeqV2` matrix, version `2017-10-13`, distributed as `log2(RSEM normalized count + 1)`. Eligibility is reconstructed from its complete **1,218-column** header and pinned clinical metadata: **157 main patients** and **113 strict HER2-IHC-negative patients**. GSE58812 (107 patients) and the newly executed GSE21653 analysis (87 receptor-defined TNBC cases) provide separate bulk cohorts. GSE176078 and GSE210616 supply single-cell and multicellular spatial context.

## Scientific scope

- TCGA eligibility requires type-01 primary tumour, explicit ER/PR IHC negativity, negative HER2 IHC or FISH evidence and no positive HER2 IHC/FISH or IHC 3+ result. Strict sensitivity adds explicit HER2-IHC negativity. One primary barcode per patient is chosen lexicographically before expression/outcome analysis. PAM50 does not determine inclusion.
- The complete documented source contains 20,530 gene rows. Additional file 3 supplies the exact all-gene 157-column subset, the complete source header, source metadata/provenance and pinned receptor table. The complete source or that SHA-verified subset can reproduce the analysis. No undocumented expression fallback is accepted and no additional log transform is applied.
- IFI44L is removed before response-gene limma fitting/moderation/BH, from the ranked enrichment list and from every correlated module. GO results retain the complete eligible-family denominator and explicitly identify unresolved tests.
- TCGA checkpoint adjustment ranks all variables, regresses out four target-excluded module scores and correlates rank residuals, with an approximate two-sided t reference. Every adjustment module excludes IFI44L and the checkpoint being tested. Nine-gene external checkpoint extensions in GSE58812/GSE21653 are **TCGA-result-driven exploratory analyses**, not independently prespecified discoveries.
- GSE21653 retains historical clinical receptor calls, with no RNA-subtype substitute. Its deposited RMA values are used directly. Primary DFS models require positive time (83 patients/25 events); two recorded zero-time events remain unchanged in a named sensitivity (85/27). The age/pT/pN/grade sensitivity has 81 complete cases/24 events and five coefficients, so it cannot establish independent prognostic performance.
- Expression modules are not measured cell fractions or purity. Dissociated-cell means do not quantify contributions to bulk RNA. Visium spots represent multicellular niches. Spatial leave-one-patient-out and bootstrap checks use the patient as the independent unit.
- Cohort-specific associations, null findings and uncertainty are retained. Numerical reproducibility does not establish causal regulation, treatment prediction or clinical utility.

No locally collected patient/IHC material, cell-line data, immunoblots or manuscript source is included. Dataset-specific terms remain separate from the code licence.

## Layout and installation

```text
scripts/                         input preparation, analysis, figures and preflight
environment/                     recorded software versions and input manifest
PUBLIC_DATA_ACCESS.md             exact sources, hashes and staging
SOURCE_DATA_MAP.tsv               figure/input-to-code map
VALIDATION_REPORT.md              executed checks and remaining release boundary
LICENSE / LICENSE_STATUS.md       MIT code terms and third-party exclusions
```

Runtime `source_metadata/`, `results/`, `sensitivity/`, `figures/` and `logs/` directories are not bundled as code. Follow `environment/README.md`; package lists record actual versions and are not a tested universal lockfile.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/python-requirements.txt
Rscript environment/install_R_packages.R
```

## Stage inputs

Extract Additional file 3 and Additional file 4 into `source_metadata/`. The former must create `TCGA_documented_primary/` and `TCGA_clinical/`; the latter must create `GSE21653/`. Archive filenames can differ, so use their actual local paths:

```sh
mkdir -p source_metadata
unzip -n /path/to/Additional_file_3.zip -d source_metadata
unzip -n /path/to/Additional_file_4.zip -d source_metadata
```

Additional file 5 supplies supplementary methods, Supplementary Figure S1, robustness outputs and annotation audits; it is not a replacement for the analysis inputs. Obtain the remaining official TCGA-CDR, GSE58812, GSE176078/Wu workbook/h5ad and GSE210616 inputs as described in `PUBLIC_DATA_ACCESS.md`.

```text
source_metadata/
├── TCGA_documented_primary/
│   ├── TCGA_BRCA_documented_primary_subset.tsv.gz
│   ├── TCGA_HiSeqV2_full_sample_header.tsv
│   ├── HiSeqV2.json
│   ├── REACQUISITION_RECEIPT.json
│   └── TCGA_documented_input_provenance.json
├── TCGA_clinical/brca_tcga_clinical_patient.txt
├── TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx
├── GSE21653/...
├── GSE58812/GSE58812_series_matrix.txt.gz
├── Wu2021_NatGenet_Supplementary_Tables.xlsx
├── GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad
└── spatial/GSE210616/...
```

## Execute

```sh
python scripts/17_reproducibility_preflight.py \
  --project "$PWD" --raw-root "$PWD" --deep-hash
python scripts/18_run_nonexperimental_pipeline.py \
  --project "$PWD" --raw-root "$PWD" --stage all --permutations 499
```

The complete runner rebuilds main TCGA eligibility/bulk/survival and exploratory categorical-stage/stage-stratified sensitivities, both external bulk cohorts, exploratory external checkpoint checks, single-cell/spatial results, patient-level spatial robustness, a separate strict-IHC TCGA analysis and the main/supplementary figures. It stops on failed steps and records logs. Strict outputs are isolated under `sensitivity/strict_IHC`; they do not replace the main cohort. `--run-census-audit` additionally queries the pinned Census metadata and does not regenerate the h5ad from raw sequencing.

To run only the new GSE21653 analysis after staging its supplied inputs:

```sh
Rscript scripts/25_gse21653_validation.R --project "$PWD"
```

`scripts/27_prepare_gse21653_inputs.py --project "$PWD" --download` can reacquire the exact full GSE21653 matrix/metadata and reconstruct the subset. `scripts/24_prepare_documented_tcga_input.py` reconstructs the documented TCGA input bundle from the exact complete source and its acquisition/clinical metadata; see `PUBLIC_DATA_ACCESS.md` for its required arguments.

## Validation and attribution

New isolated TCGA main/strict models and the complete-source/subset GSE21653 executions have bounded retained checks. GSE21653 module statistics and all four Cox-model coefficients were independently recalculated; 16 result TSVs and its complete normalized gene matrix matched between input routes. External checkpoint checks have a separate exploratory audit. The complete final runner was subsequently executed successfully, including the added stage sensitivity and all figures. The independent final replay matched all 124 analysis TSVs and all frozen analysis scripts. See VALIDATION_REPORT.md for its local-environment and shared-input boundaries.

Author attribution and contributions remain deferred. The hosting account is not evidence of manuscript authorship. Original code/documentation uses MIT; third-party inputs, publications, images and dependencies retain their own terms. The final manuscript and package must pin the same externally verified immutable commit. This tree intentionally does not embed its own content-derived commit hash.

The TCGA stage sensitivities use 152 complete cases. The categorical-stage OS term violated PH, motivating stage-stratified analyses for both endpoints; all specifications are retained and do not replace the primary continuous models. See scripts/29_tcga_stage_sensitivity.R and complete source tables.
