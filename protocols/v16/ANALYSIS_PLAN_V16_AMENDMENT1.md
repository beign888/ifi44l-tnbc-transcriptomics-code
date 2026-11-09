# Clarification before new model fitting

The original plan remains unchanged at SHA256 `8d24b3a4971f981ad4651d7f332e229da8ec2e25b9695af67e3b8961d947c97e`.

The ten fixed expression-redundancy correlations per cohort use two-sided Spearman tests with the same asymptotic/no-exact-ties convention as v15, followed by BH across ten tests per cohort. Pearson and Spearman values reported for the Cox predictor diagnostics do not form another significance-screening family. This clarification was made before script30 implementation/model fitting; it does not alter targets, cohorts, endpoints or covariance specifications.
