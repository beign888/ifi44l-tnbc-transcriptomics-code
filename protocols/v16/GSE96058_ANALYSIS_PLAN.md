# Additional independent RNA-seq cohort GSE96058

Decision recorded after clinical-metadata feasibility and before acquiring/inspecting IFI44L expression or fitting expression-outcome associations in this cohort. This is an exploratory revision extension, not public preregistration. It does not alter the already fixed five primary/three secondary conditional tests in the existing three cohorts.

## Inclusion and processing

Use both official platform metadata files and the official gene-expression file GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz. Link expression by original F sample title, retain the primary F record and exclude F*repl technical replicates without using expression/outcomes. Use only source clinical er status=0, pgr status=0 and her2 status=0, excluding missing/positive clinical calls. Do not substitute RNA-predicted MGC/SGC receptor status or PAM50. Clinical feasibility identifies 143 eligible primary cases, with 143 positive-time OS records/26 events and 137 complete pathological cases/26 events. Retain an all3409 decision ledger, sample linkage, duplicates/exclusions and missingness.

The source documents log2(FPKM+0.1) values, summed by symbol from transcripts. Use them directly with no second log transform or cross-platform pooling. Preserve an exact all-gene subset as source decimal strings; keep its hash, full-source hash/header, both complete metadata sources, fixed membership and download provenance. Historical source diagnostic receptor calls are used; do not infer a uniform modern <1% hormone-receptor cutoff from the binary GEO fields.

## Fixed analyses

Score all 13 existing IFI44L-excluded modules from fixed gene definitions using the same ssGSEA method, with cohort-specific gene coverage. Report all 13 Spearman correlations with BH13. No pathway/cutpoint selection from this cohort. Standardize IFI44L and IFN-alpha scores over all143 eligible expression samples.

Fit seven Efron-tie Cox specifications for OS (positive days and event0/1): (1) unadjusted IFI44L, (2) age/10+IFI44L, (3) clinical base+IFI44L, (4) age/10+IFN-alpha, (5) age/10+IFN-alpha+IFI44L, (6) clinical base+IFN-alpha and (7) clinical base+IFN-alpha+IFI44L. The clinical base is age/10, tumour size(mm)/10, node positive versus negative and NHG G3 versus G1/2. No treatment-response or treatment-benefit analysis. Models4/5 use common complete cases; models3/6/7 use common complete clinical cases. Retain all coefficients, PH, rank/VIF, warnings, exclusions, scaling and sample inputs. Clinical models have only26 events for5–6 coefficients and are explicitly event-limited.

Two one-df nested LRTs compare addition of IFI44L to the age+IFN-alpha base and to the clinical+IFN-alpha base. BH adjustment is across these two fixed, additional-cohort exploratory tests. They are separate from the original5/3 families; all nominal and adjusted P values remain available. The cohort was chosen for its source-linked clinical annotation and independent RNA-seq platform, not observed IFI44L significance. No inference from nonsignificance to absence, no validated incremental prediction, no optimized risk score or new high/low threshold.

## Deliverables and boundaries

Provide the actual input preparation/download/validation helper, R model code, full/subset route verification where feasible, full source hashes/metadata, complete result tables and a source-bound Supplementary Figure S4. Existing v15, the new three-cohort analyses and local experimental records remain unchanged. The cohort may be added only after successful input identity/processing/clinical provenance checks; otherwise preserve failure evidence and do not claim completed validation.
