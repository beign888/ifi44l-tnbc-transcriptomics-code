# v16 targeted analysis plan

Prepared after reviewing v15 results and before fitting the new models. This is a local prospective analysis decision record for this revision, not public preregistration or an independent discovery design.

## Question

Does the existing association of IFI44L with clinical outcomes contain detectable information conditional on the broader interferon-associated expression state? Are the observed expression patterns shared with other fixed interferon-stimulated genes?

## Primary conditional association analyses

Use the frozen documented TCGA157, GSE58812 107 and GSE21653 87 expression cohorts. Retain prior endpoint/time rules and receptor definitions without outcome-based reselection. Five endpoint/cohort pairs are fixed: TCGA OS/PFI, GSE58812 MFS/OS, GSE21653 positive-time DFS. Standardize IFI44L and the existing IFI44L-excluded 19-gene IFN_ALPHA_RESPONSE ssGSEA score over each full eligible expression cohort, before outcome complete-case exclusion; do not pool platforms. Age remains per ten years.

For each endpoint, fit on exactly the same complete cases: age; age+IFI44L; age+IFN-alpha; age+IFN-alpha+IFI44L. Use unpenalized Efron-tie Cox regression and preserve all coefficients, confidence intervals, term/global Schoenfeld diagnostics, warnings, events and model definitions. The one-degree-of-freedom likelihood-ratio comparison is the last model against age+IFN-alpha. BH correction is applied across the five primary conditional tests. Retain estimates regardless of significance or sign. No new threshold, score optimization, mediation, causal or validated predictive claim.

Three named secondary comparisons are fixed: TCGA OS/PFI with age and categorical-stage strata (the previous stage-PH issue is already known); GSE21653 DFS with age, pT, pN and grade. Add IFN-alpha then IFI44L to the same clinical base; use common complete cases and retain these event-limited results. BH correction is across these three secondary comparisons separately. They do not replace the primary models. No further interaction or subgroup tests are selected by results.

Report predictor Pearson/Spearman correlations, standardized-design condition numbers, variance-inflation factors and rank. State that a nonsignificant conditional test, particularly under near-collinearity and sparse events, is not evidence of no information. Likelihood-ratio improvement is a within-cohort association check, not external predictive validation. Do not report unvalidated AUC/C-index as clinical performance.

## Fixed expression comparators

Comparators are IFI44, MX1 and ISG15, chosen now for biological relationship/representative interferon response, without screening their survival significance. Compare all six gene pairs in each cohort. Also rescore a common IFN module removing IFI44L and all three comparators from the original 19-gene list (16 remaining genes), then compare each of the four transcripts with that module. Keep all ten correlations per cohort, with BH across ten within each cohort; preserve gene coverage and scores. These are exploratory redundancy checks, not a search for a superior prognostic marker.

## Additional-data feasibility

Assess METABRIC receptor provenance before accepting it as TNBC validation; do not silently equate HER2 copy-number surrogate with negative clinical HER2 IHC/FISH. Assess author-provided GSE161529 cell annotations and donor/sample definitions before treating it as independent cellular replication. Any executed new-cohort analysis needs its own recorded design, source hashes and inference-unit rules. Acquisition/feasibility alone is not a completed validation.

## Evidence and deliverables

Preserve v15 unchanged. New scripts, actual outputs, diagnostics and independent audits go in the versioned v16 directories; all results, including nulls, must be reflected in the new manuscript/supplement and mirrored to T7. Author information/contributions remain deferred. No local IHC/cell-line evidence is reintroduced and no manuscript is submitted.
