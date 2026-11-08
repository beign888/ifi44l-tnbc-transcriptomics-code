# Public code snapshot validation report

Validation date: 8 November 2026

Snapshot: `v1.1.0`

## Checks completed

| Check | Result | Scope |
|---|---|---|
| Python syntax compilation | PASS | All Python scripts |
| R syntax parsing | PASS | Analysis scripts and the R package installer |
| TCGA-CDR file integrity | PASS | MD5 `a4591b2dcee39591f59e5e25a6ce75fa` in preserved workspace |
| TCGA OS/PFI join | PASS | 187 expression samples mapped to 187 unique patients |
| TCGA survival models | PASS | OS 24 events; PFI 30 events; unadjusted and age-adjusted continuous models |
| Revised Figure 2 build | PASS | GSE58812 and TCGA survival results displayed together, including non-significant TCGA estimates |
| Revised Figure 5 build | PASS | Public-multiomics descriptive synthesis; no local or cell-line branch |
| Input and environment preflight against preserved workspace | PASS | Declared public inputs and required packages available |
| Author/contact information scan | PASS | No deferred author names, contribution statement or correspondence email detected |
| Local-path scan | PASS | No local absolute filesystem path or removable-volume path detected |
| Common secret-pattern scan | PASS | No password, API-key, token or private-key pattern detected |
| Redistribution boundary | PASS | No third-party raw data, local patient/IHC material, cell-line measurements or manuscript source included |

## Numerical result boundary

TCGA continuous IFI44L estimates were directionally favourable but non-significant: OS unadjusted HR 0.829 (95% CI 0.546–1.259; P=0.378), OS age-adjusted HR 0.835 (0.554–1.259; P=0.389), PFI unadjusted HR 0.839 (0.582–1.209; P=0.346), and PFI age-adjusted HR 0.837 (0.579–1.211; P=0.346). Median-split curves are graphical sensitivity analyses only.

Pathologic-stage adjustment was not promoted because stage-I OS events were absent and stage-I PFI events were sparse. The analysis therefore does not claim independent prognostic value.

## Interpretation boundary

These checks establish static executability and a bounded local rerun against preserved public inputs. They do not constitute a fully independent numerical reproduction from newly downloaded source files, scientific peer review, clinical validation or submission readiness.
