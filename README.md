# IFI44L in triple-negative breast cancer: reproducibility code

Version: `v1.2.1` corrected clinical-cohort snapshot with a documented TCGA expression-processing sensitivity, prepared 9 November 2026.

Repository: <https://github.com/beign888/ifi44l-tnbc-transcriptomics-code>. A new public commit is not asserted here. Use the exact release commit recorded in the final release/package manifest and manuscript; verify the sealed file hashes in `FILE_MANIFEST_SHA256.tsv`. The previous v1.1.0 commit does not identify these corrected analyses.

This code covers TCGA-BRCA bulk transcriptomics and survival, independent GSE58812 bulk/survival analysis, GSE176078 single-cell RNA sequencing and GSE210616 Visium spatial transcriptomics. The corrected main TCGA cohort contains **151 unique patients**; a separate **110-patient strict HER2-IHC-negative sensitivity cohort** is supported. The new processing sensitivity uses the **same 151 patients** in a separately documented UCSC Xena expression representation. It is not an independent cohort. The historical 187-sample helper list is superseded and must not be described as a verified TNBC cohort.

## Scientific and provenance boundaries

- Main eligibility requires explicit ER and PR negativity, negative HER2 IHC or FISH evidence, and no positive HER2 IHC/FISH or IHC 3+ result. Strict eligibility additionally requires explicit HER2 IHC negativity. PAM50 is not used for inclusion. Selection is generated before outcome modelling from a commit- and SHA-pinned public clinical table.
- The exact 151-column archived TCGA-derived expression subset is supplied in two journal supporting-data archives, Additional files 3 and 4. Upstream expression release, units and original normalization cannot be recovered from the archived records. No TPM, FPKM, RSEM or log2-fold-change provenance is claimed. The preserved values are used without an additional transformation.
- Additional file 5 supplies a separate all-gene, 151-column subset of UCSC Xena `TCGA.BRCA.sampleMap/HiSeqV2`, version `2017-10-13`, distributed as `log2(norm_count+1)` from RSEM normalized results. This input is used only to check module correlations under documented expression processing; it does not identify the original archived matrix's source or replace the main analyses. No new GSEA, survival or adjusted checkpoint analysis is run in this sensitivity.
- IFI44L is excluded from the response-gene family **before** limma fitting, empirical-Bayes moderation and BH adjustment, and from the ranked enrichment list and all correlated signatures.
- In TCGA, primary checkpoint adjustment ranks all variables before OLS residualization, correlates rank residuals using Pearson correlation, and uses an approximate two-sided t reference with `n-k-2` degrees of freedom. Each adjustment signature excludes IFI44L and the checkpoint currently tested (leave-one-target-out). Shared-module rank adjustment and the historical raw-scale residual/Spearman calculation are separate TCGA sensitivity outputs. GSE58812 checkpoint analyses are unadjusted Spearman correlations with BH correction; no partial-correlation or leave-one-target-out checkpoint adjustment was performed in that cohort.
- Expression-module scores are not measured cell fractions or tumour purity. Dissociated single-cell expression/detection does not quantify each compartment's contribution to bulk expression. Visium spots represent multicellular niches.
- The inferential unit is the tumour, patient or donor. Continuous standardized IFI44L is the primary survival predictor; median groups are graphical sensitivities. TCGA results are retained regardless of significance, without claiming independent prognostic or treatment-predictive utility.

The repository contains no local patient/IHC material, cell-line measurements, immunoblots or manuscript source. Public TCGA-derived supporting inputs are distinct from a locally collected clinical cohort. Original dataset terms continue to apply.

## Layout and environment

```text
scripts/                    restoration, analysis, figure and validation code
environment/                recorded Python/R versions and installation instructions
PUBLIC_DATA_ACCESS.md       sources, hashes, staging and provenance limitations
SOURCE_DATA_MAP.tsv         figure/input-to-source-to-script map
VALIDATION_REPORT.md        completed checks and remaining validation boundaries
FILE_MANIFEST_SHA256.tsv    file integrity manifest, regenerated when sealed
LICENSE                     MIT terms for original code/documentation
LICENSE_STATUS.md           third-party exclusions
```

Runtime directories `source_metadata/`, `results/`, `figures/` and `logs/` are excluded from version control. Follow the recorded versions in `environment/`; the audit replay used Python 3.12.13 and preserved packages. A clean installation on every operating system has not been established.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/python-requirements.txt
Rscript environment/install_R_packages.R
```

## Restore and stage inputs

Obtain both `Additional_file_3_TCGA_Frozen_Input_Part1.zip` and `Additional_file_4_TCGA_Frozen_Input_Part2.zip` from the manuscript's supporting files. From the repository root:

```bash
mkdir -p source_metadata
unzip -n /path/to/Additional_file_3_TCGA_Frozen_Input_Part1.zip -d source_metadata
unzip -n /path/to/Additional_file_4_TCGA_Frozen_Input_Part2.zip -d source_metadata
python scripts/00_restore_frozen_tcga_input.py --project "$PWD"
```

The script concatenates `.part01` and `.part02` in that order and checks SHA-256. The restored gzip is 27,173,060 bytes. Both parts are required; neither is independently readable as a gzip. Part 1 also supplies the pinned clinical table, original 1,072-sample header and provenance receipt.

```text
source_metadata/
├── TCGA_frozen/
│   ├── TCGA_BRCA_TNBC_frozen.tsv.gz.part01
│   ├── TCGA_BRCA_TNBC_frozen.tsv.gz.part02
│   ├── TCGA_BRCA_TNBC_frozen.tsv.gz
│   ├── TCGA_BRCA_expression_header_samples.tsv
│   └── TCGA_frozen_input_provenance.json
├── TCGA_clinical/brca_tcga_clinical_patient.txt
├── TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx
├── GSE58812/GSE58812_series_matrix.txt.gz
├── Wu2021_NatGenet_Supplementary_Tables.xlsx
├── GSE176078_Wu2021_CELLxGENE_2025-11-08.h5ad
└── spatial/GSE210616/
    ├── GSE210616_family.soft.gz
    ├── GSE210616_primary28_samples.tsv
    └── files/GSM*/...
```

Stage other public files as described in `PUBLIC_DATA_ACCESS.md`. Historical receptor workbooks and `三阴.txt` are not required by the corrected default workflow. If the original `symbol.BRCA.txt` is supplied under a separate `--raw-root`, the scripts prefer it over the restored subset; verify its recorded hash rather than silently substituting another matrix.

For the processing sensitivity, extract Additional file 5 and copy the contents of its `TCGA_processing_sensitivity/` directory into `source_metadata/TCGA_documented_processing/`. That directory must contain `documented_HiSeqV2_matched151_expression.tsv.gz`, the four `baseline_*` inputs and `PROCESSING_SENSITIVITY_DESIGN.json`. The supplied subset is 8,251,822 bytes (SHA-256 `53d57d8013a86933143964d7ffd68d4247f1bf8737c0c2ce5238e97c2a9a85f8`) and preserves all 20,530 source gene rows and the frozen 151-sample order. The complete 64,286,337-byte Xena source is optional for this replay. See `PUBLIC_DATA_ACCESS.md` for its identity and the subset-preparation helper.

```bash
python scripts/08_prepare_gse210616_spatial.py \
  --soft source_metadata/spatial/GSE210616/GSE210616_family.soft.gz \
  --output-dir source_metadata/spatial/GSE210616 --download
python scripts/17_reproducibility_preflight.py \
  --project "$PWD" --raw-root "$PWD" --deep-hash
```

## Reproduction commands

Run the main 151-patient workflow after staging every required input:

```bash
python scripts/18_run_nonexperimental_pipeline.py \
  --project "$PWD" --raw-root "$PWD" --stage all --permutations 499 \
  --run-processing-sensitivity
```

This rebuilds clinical selection, bulk/survival, single-cell and spatial results and Figures 1–5, then runs the separately staged processing sensitivity. It stops at a failed step and records logs in `logs/reproduction/`. Omit `--run-processing-sensitivity` to run the original core workflow without Additional file 5. The flag is available with `--stage all` or `--stage public-analyses`. Add `--run-census-audit` only when access to the pinned Census version is available; this optional metadata query does not regenerate raw sequencing or the processed h5ad.

The processing sensitivity can also be run on its own, without rerunning the core analyses:

```bash
Rscript scripts/22_documented_tcga_processing_sensitivity.R \
  --work-dir source_metadata/TCGA_documented_processing \
  --output-dir results/processing_sensitivity
```

Its default input is the supplied small subset. Use `--input` to select the complete documented source matrix instead. The supplied baseline files permit a standalone comparison; `--baseline-project "$PWD"` optionally reads them from `results/bulk`. The script uses the original 13 module definitions with IFI44L excluded. TRAC is absent in this legacy annotation, so the CD8 T-cell module has 7/8 members; the IFN-alpha-response module retains 19/19. Scores use the original ssGSEA settings and correlations use two-sided Spearman tests with BH correction across all 13 modules.

For the 110-patient strict-IHC sensitivity, use a **separate copy of the code project with staged inputs**, then run:

```bash
python scripts/21_audit_tcga_tnbc_definition.py \
  --project "$PWD" --raw-root "$PWD" --strict-ihc
Rscript scripts/01_bulk_tnbc_reanalysis.R --project "$PWD" --raw-root "$PWD"
Rscript scripts/06_tcga_tnbc_survival.R "$PWD"
```

Do not run `--stage all` after strict selection in the same results directory: its default selection restores the main cohort. Keep main and sensitivity results separate.

## Validation and citation

The core numerical evidence is inherited from the bounded v1.2.0 checks: corrected main/strict bulk snapshots, GSE58812, 19 single-cell/spatial tables, and comparison of the corrected-cohort Python Kaplan–Meier helper against R survival, log-log confidence intervals and risk sets. A separate frozen-input route restored the exact eligible sample set without the original matrix or helper list; raw expression, core limma results and the survival script outputs matched the original-input route. This frozen-route check did not rerun GSEA. Subsequent label/diagnostic-only bulk edits were checked separately; they are not a second full end-to-end replay. Some enrichment tests remain explicitly unresolved.

The new documented-processing sensitivity was separately executed against the complete source and the supplied subset: both exited zero, ten result TSVs matched byte for byte, and an independent SciPy recalculation verified all 13 correlations and BH values. Source-string comparison checked 3,100,030 values. This is a bounded processing check on overlapping patients, not a fresh end-to-end execution of v1.2.1 or independent biological validation. Core analysis algorithms and results from scripts 01–21 are unchanged; the runner adds an optional call to script 22, while script 23 reconstructs its documented subset when needed. See `VALIDATION_REPORT.md` for the exact execution and release-check boundaries.

Citation metadata, author attribution and contributions remain deferred. The GitHub hosting account is not evidence of manuscript authorship. Original code/documentation uses MIT; third-party data, publications, images and dependencies retain their own terms (`LICENSE_STATUS.md`).
