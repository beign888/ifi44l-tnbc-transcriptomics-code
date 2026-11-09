# IFI44L in triple-negative breast cancer: reproducibility code

Version: `v1.4.0`, prepared 9 November 2026. The recorded extension replay completed six jobs successfully, with 29 scripts and eight protocol files matching the frozen execution tree. The immutable release commit and remote verification are recorded outside this tree in the manuscript and a separate release receipt; the tree does not embed its own content-derived commit.

Repository: <https://github.com/beign888/ifi44l-tnbc-transcriptomics-code>. The 22 core analysis/preparation scripts are retained byte-for-byte from v1.3.0. Its completed local all-stage execution and 124-table byte comparison are **existing v15/v1.3.0 evidence**. The v16 extension replay is a separate completed record: six jobs exited zero and all 58 compared analysis/plot TSVs matched byte-for-byte. Neither record is a clean installation or cross-platform test. Seven extension scripts (30–36) bring the supplied total to 29 scripts.

The primary TCGA analysis uses the documented UCSC Xena `TCGA.BRCA.sampleMap/HiSeqV2` matrix, version `2017-10-13`, distributed as `log2(RSEM normalized count + 1)`. Eligibility is reconstructed from its complete **1,218-column** header and pinned clinical metadata: **157 main patients** and **113 strict HER2-IHC-negative patients**. GSE58812 (107 patients) and the newly executed GSE21653 analysis (87 receptor-defined TNBC cases) provide separate bulk cohorts. GSE176078 and GSE210616 supply single-cell and multicellular spatial context.

## Scientific scope

- TCGA eligibility requires type-01 primary tumour, explicit ER/PR IHC negativity, negative HER2 IHC or FISH evidence and no positive HER2 IHC/FISH or IHC 3+ result. Strict sensitivity adds explicit HER2-IHC negativity. One primary barcode per patient is chosen lexicographically before expression/outcome analysis. PAM50 does not determine inclusion.
- The complete documented source contains 20,530 gene rows. Additional file 3 supplies the exact all-gene 157-column subset, the complete source header, source metadata/provenance and pinned receptor table. The complete source or that SHA-verified subset can reproduce the analysis. No undocumented expression fallback is accepted and no additional log transform is applied.
- IFI44L is removed before response-gene limma fitting/moderation/BH, from the ranked enrichment list and from every correlated module. GO results retain the complete eligible-family denominator and explicitly identify unresolved tests.
- TCGA checkpoint adjustment ranks all variables, regresses out four target-excluded module scores and correlates rank residuals, with an approximate two-sided t reference. Every adjustment module excludes IFI44L and the checkpoint being tested. Nine-gene external checkpoint extensions in GSE58812/GSE21653 are **TCGA-result-driven exploratory analyses**, not independently prespecified discoveries.
- GSE21653 retains historical clinical receptor calls, with no RNA-subtype substitute. Its deposited RMA values are used directly. Primary DFS models require positive time (83 patients/25 events); two recorded zero-time events remain unchanged in a named sensitivity (85/27). The age/pT/pN/grade sensitivity has 81 complete cases/24 events and five coefficients, so it cannot establish independent prognostic performance.
- Expression modules are not measured cell fractions or purity. Dissociated-cell means do not quantify contributions to bulk RNA. Visium spots represent multicellular niches. Spatial leave-one-patient-out and bootstrap checks use the patient as the independent unit.
- Cohort-specific associations, null findings and uncertainty are retained. Numerical reproducibility does not establish causal regulation, treatment prediction or clinical utility.

## Version 16 extension scope

- `30_IFI44L_IFN_conditional_context.R` preserves five primary cohort/endpoint comparisons: TCGA OS/PFI, GSE58812 MFS/OS and GSE21653 positive-time DFS. Within each comparison, four Cox specifications use the same complete cases. IFI44L and the existing IFI44L-excluded 19-gene IFN-alpha score are standardized over the full expression cohort. The addition of IFI44L to age plus IFN-alpha is tested by a one-degree-of-freedom likelihood-ratio test (LRT), with BH correction across all five tests. Three secondary clinical-context comparisons have a separate BH family of three; none replaces a primary model. All coefficients, exclusions, event counts, scaling, rank/VIF/condition diagnostics, PH tests and warnings are retained. These correlated, event-limited models are conditional association checks, not evidence of incremental predictive performance.
- Fixed expression redundancy uses IFI44L, IFI44, MX1 and ISG15. All six gene pairs and all four gene-versus-module comparisons are retained. The latter use a shared 16-gene IFN-alpha module excluding all four transcripts. Two-sided Spearman tests use BH correction across the ten fixed tests within each cohort; no survival-based comparator selection is performed. `31_make_IFI44L_IFN_context_figure.py` reads saved results without re-fitting models and exports Supplementary Figure S2 plus exact plotting-source tables.
- `32_single_cell_within_lineage_IFN_context.py` examines donor-level IFI44L and the existing source IFN-response score within fixed author-defined major lineages in GSE176078. This single-cell score is distinct from the bulk 19-gene module. Minimum-cell eligibility, unavailable feature handling and the full nine-lineage multiplicity family are explicit. An exact donor-label permutation sensitivity was added after an extreme small-sample asymptotic result was observed; the original results and dated amendment are retained. Supplementary Figure S3 describes within-cohort context, not a new independent single-cell cohort, causal regulation or absolute cellular contributions to bulk RNA.
- The additional SCAN-B RNA-seq cohort, GSE96058, is selected using original clinical ER/PGR/HER2-negative calls, excluding RNA-predicted calls and explicit technical replicates. Clinical-source feasibility identifies 143 eligible primary records (26 OS events) and 137 complete clinical records (26 events). The deposited `log2(FPKM+0.1)` values require no second log transform. Scripts 34–36 prepare the exact all-gene subset, execute the fixed validation and create Supplementary Figure S4. Complete source acquisition and exact all-gene subset preparation are recorded. The complete source, exact subset, all nine Cox models, three LRTs and 13 module tests were actually executed and source-bound. Independent software checks passed 302/302; Figure S4 passed 181 plotting/source/font checks and actual visual review. All IFI44L survival tests and additional-cohort LRTs were nonsignificant; this does not establish absence of information. The original seven Cox specifications are retained; a source-platform amendment frozen before association fitting adds two instrument-stratified models for the 126 HiSeq 2000 and 17 NextSeq 500 cases. The final design therefore has nine fixed Cox specifications and three additional-cohort LRTs with BH correction across all three, separate from the original five/three testing families.
- Locked plans and amendments are staged under `protocols/v16/`. These are local, prospectively frozen or explicitly diagnostic-driven revision plans as identified in each document; they are not public preregistration. The extension wrapper records execution separately from the unchanged core runner.

No locally collected patient/IHC material, cell-line data, immunoblots or manuscript source is included. Dataset-specific terms remain separate from the code licence.

## Layout and installation

```text
scripts/                         22 unchanged core scripts plus seven v16 extension scripts
protocols/v16/                   locked revision plans and explicit amendments
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

Additional file 5 supplies the core supplementary methods, Supplementary Figure S1, robustness outputs and annotation audits. Additional file 6 is designated for the v16 supplementary methods, Figures S2–S4 and complete extension result/source tables; Additional files 7–8 jointly carry the new frozen GSE96058 input bundle. The exact 25,592,924-byte gene gzip is transported losslessly as two binary volumes; both files are required. These final supporting-file archives require their own checksum and content verification; results files do not substitute for frozen analysis inputs. Obtain the remaining official TCGA-CDR, GSE58812, GSE176078/Wu workbook/h5ad and GSE210616 inputs as described in `PUBLIC_DATA_ACCESS.md`.

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
├── spatial/GSE210616/...
└── GSE96058/...                   additional frozen input bundle (Additional files 7–8)
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

## Execute v16 extensions after the core workflow

After the core results, all seven extension scripts and locked plans have been staged, run the extension-only wrapper. It does not re-run or silently substitute for the core workflow above:

```sh
python scripts/33_run_v16_context_extensions.py --project "$PWD" \
  --plan-dir "$PWD/protocols/v16" --stage context
python scripts/33_run_v16_context_extensions.py --project "$PWD" \
  --plan-dir "$PWD/protocols/v16" --stage scanb
# --stage all executes both extension groups after the same prerequisites.
```

Extract both Additional files 7 and 8 into `source_metadata/` so their contents share `source_metadata/GSE96058/`. Additional file 7 contains metadata/provenance and `.part01`; Additional file 8 contains `.part02`. Restore the exact original gzip offline before analysis:

```sh
python scripts/34_prepare_GSE96058_inputs.py --project "$PWD" --restore-parts
```

The helper verifies the reconstructed size/SHA-256, uses no network and retains source parts. The final extension wrapper prepares this restoration before its SCAN-B step. Alternatively, reacquire the official source and reconstruct the subset as described in PUBLIC_DATA_ACCESS.md. Script 34 is an input-preparation helper and is not inferred to have run merely because it is supplied. Logs and the extension execution receipt are written under `logs/v16_extensions/`. A completed extension wrapper does not by itself demonstrate a full core-plus-extension run on a clean machine.

## Validation and attribution

The existing v1.3.0 core replay matched all 124 compared analysis TSVs in the recorded local environment with shared verified inputs. The isolated v16 conditional analysis has a separate 421-check numerical audit, including independent statsmodels PHReg fits of all 29 Cox models. IFNa16 ssGSEA and R Schoenfeld tests are not independently reimplemented by that audit. A separate 38-check single-cell audit reconstructs source aggregates and exact donor-permutation counts using an independent dynamic-programming implementation. These are bounded computational checks, not independent biological validation. See VALIDATION_REPORT.md for observed evidence and items still requiring final release records.

SCAN-B full-source and exact-subset R routes both exited zero with the final script; all 20 compared analysis TSVs and the normalized gene-expression gzip were byte-identical. All 4,413,695 selected expression decimal strings and 30,865 gene identifiers matched the complete official source exactly. Offline restoration reproduced the exact 25,592,924-byte gzip. The complete extension wrapper subsequently passed six jobs and all 58 compared analysis/plot TSVs, with three path/output-manifest tables explicitly excluded. This used the existing macOS runtime and shared immutable core prerequisites, and was not a new full-core replay, clean installation, cross-platform test or independent biological validation. Final repository/remote identity and supporting-archive checks are recorded by the external release process.

Author attribution and contributions remain deferred. The hosting account is not evidence of manuscript authorship. Original code/documentation uses MIT; third-party inputs, publications, images and dependencies retain their own terms. The final manuscript and package must pin the same externally verified immutable commit. This tree intentionally does not embed its own content-derived commit hash.

The unchanged core TCGA stage sensitivities use 152 complete cases. The categorical-stage OS term violated PH, motivating stage-stratified analyses for both endpoints; all specifications are retained. The new conditional clinical sensitivities retain their own complete-case/event limits and PH diagnostics. Neither replaces the primary core continuous models or establishes independent prognosis.
