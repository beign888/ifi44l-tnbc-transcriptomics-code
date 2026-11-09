# Recorded analysis environment and v1.4.0 extension boundary

`python-requirements.txt` and `R-packages.tsv` record package versions used by the preserved core analysis environment. Version1.4.0 retains those files and the22 unchanged core scripts fromv1.3.0. The software/input manifests now record the actual SCAN-B source, split-input restoration and completed six-job extension replay, separately from historical core execution. None of these files is a universal lockfile.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/python-requirements.txt
Rscript environment/install_R_packages.R
```

The R helper installs missing packages through the configured CRAN/Bioconductor release and asks the user to compare resolved versions with the recorded table. It does not guarantee historical binaries. The core records use R 4.3.2, including the listed `survival`, `GSVA`, `limma` and annotation versions. Existing GSE21653 independent validation used Python 3.12.14, NumPy 2.5.2, pandas 2.3.3, SciPy 1.18.1 and statsmodels 0.15.0. Preserve the actual per-execution session records when comparing results.

## Core versus v16 extensions

The complete v1.3.0 core runner passed in the existing local environment: 124 compared numerical TSVs matched, and all 22 supplied scripts matched the execution snapshot (17 were invoked). It used verified shared source inputs and was not a fresh installation or cross-platform test. Keeping those script bytes does not transfer completion status to the new extension wrapper.

The conditional extension uses R data.table, survival, GSVA, jsonlite and digest; its saved-result figure uses Python NumPy, pandas and matplotlib. The single-cell extension uses anndata, h5py, NumPy, pandas, SciPy and matplotlib and reads the targeted H5AD produced by the core workflow. The GSE96058 preparation helper uses the Python standard library; the final R analysis and plotting code must be checked against the recorded package lists and their actual runtime reports before sealing. No desktop office application is needed to execute analyses or make figures. Optional CELLxGENE Census access is only required for the separate core Census audit.

After staging all core inputs/results, the locked `protocols/v16/` files and both frozen GSE96058 transport archives (Additional files 7–8), use:

```sh
python scripts/33_run_v16_context_extensions.py --project "$PWD" \
  --plan-dir "$PWD/protocols/v16" --stage all
```

The exact GSE96058 input gzip is distributed losslessly across two binary parts. Stage both archives in the same source_metadata tree; script 34 supports `--restore-parts` for SHA-verified offline restoration, and the final wrapper prepares this before SCAN-B analysis. Source parts are retained and no second log transformation is applied. The source-platform amendment retains the original seven models and adds two instrument-stratified models, with all three LRTs corrected together by BH3.

The wrapper records separate extension logs and an execution receipt under `logs/v16_extensions/`. A successful context-only or SCAN-B-only run does not imply the other stage completed. A successful extension-only run does not imply a new execution of the unchanged core or a clean-machine installation.

The final extension-only stage-all replay actually completed six jobs with exit0. Its29scripts and8protocols matched the frozen tree, and58 compared analysis/plot TSVs matched byte-for-byte; three path/output-manifest tables were excluded. This used the existing macOS runtime and shared immutable core prerequisites. SCAN-B also has complete-source/exact-subset execution and302 independent source/numerical checks. The historical124-table core replay was not repeated or relabelled as a new clean core-plus-extension installation. Runtime/session files, final manifests and the external release receipt define the observed scope. No cross-operating-system guarantee or new biological validation is implied.
