# Public code snapshot validation report

Validation date: 2 November 2026

Snapshot: `v1.0.0`

## Checks completed

| Check | Result | Scope |
|---|---|---|
| Python syntax compilation | PASS | All 16 Python scripts |
| Python command-line entry points | PASS | `--help` executed for all 16 scripts in the recorded project virtual environment |
| R syntax parsing | PASS | Two analysis scripts plus the R package installer |
| C compilation | PASS | `12_extract_targeted_mtx.c` compiled with zlib |
| Input and environment preflight against the preserved local project | PASS | All declared source inputs present; 28/28 spatial sections complete; declared Python and R packages available |
| Existing expected-output inventory | PARTIAL | 11 of 12 public-pipeline output paths already present; the remaining Figure 5 path reflects a renamed public-snapshot output |
| Figure 5 public-snapshot build | PASS | PDF, SVG, PNG, TIFF and contract JSON generated successfully in an isolated temporary directory |
| Author/contact information scan | PASS | No deferred author names, contribution statement or correspondence email detected |
| Local-path scan | PASS | No local absolute filesystem path or T7 volume path detected |
| Common secret-pattern scan | PASS | No password, API-key, token or private-key pattern detected |
| Redistribution boundary | PASS | No third-party raw data, local patient/IHC material, cell-line measurements or manuscript source included |

## Interpretation boundary

These checks establish static executability, input/environment availability on the preserved development machine and code-package privacy/integrity. They are not a complete independent numerical rerun from newly downloaded public data. The full pipeline is designed to stop on missing inputs or failed steps and to rerun its output preflight at completion.

The existing-project preflight reported one absent path because the public snapshot writes the conceptual model as `Figure5_compartment_dependent_model_no_local_human_cohort.pdf`, whereas the preserved working project used an earlier filename. The public Figure 5 builder was therefore executed separately in an isolated directory and completed successfully; no scientific output in the preserved working project was overwritten.
