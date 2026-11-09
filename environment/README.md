# Recorded analysis environment

`python-requirements.txt` and `R-packages.tsv` record package versions used by the preserved analysis environment. `software-manifest.json` distinguishes actual bounded execution from an untested fresh installation; these files are not a universal lockfile. `input-manifest.json` identifies current source inputs and supporting-file roles.

The GSE21653 R runs used R4.3.2 and the recorded R packages. Its independent Python validation used Python3.12.14, NumPy2.5.2, pandas2.3.3, SciPy1.18.1 and statsmodels0.15.0. Other executed components retain their own session records. Match those versions for numerical comparison; automatically installing a newer package does not establish equivalence.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/python-requirements.txt
Rscript environment/install_R_packages.R
```

The R helper installs missing packages using the configured CRAN/Bioconductor release and asks the user to compare resolved versions with the recorded table. It does not silently promise an exact recreation of historical CRAN binaries. No desktop office application is required for analysis or figure generation. Optional CELLxGENE Census access is required only for `--run-census-audit`.

The complete final v1.3.0 all-stage runner passed in the recorded local environment; all 124 compared numerical TSVs and all 22 supplied script hashes matched the retained reference/execution snapshot. This is not a fresh installation or cross-platform test. No temporary workstation Python paths are required by the public code.
