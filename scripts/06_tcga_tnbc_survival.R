#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(survival)
  library(digest)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(flag, default = NULL) {
  position <- match(flag, args)
  if (is.na(position)) return(default)
  if (position == length(args)) stop("Missing value after ", flag)
  args[[position + 1L]]
}
positional_project <- if (length(args) && !startsWith(args[[1]], "--")) args[[1]] else "."
project_dir <- normalizePath(argument_value("--project", positional_project), mustWork = TRUE)

expression_file <- file.path(project_dir, "results", "bulk", "tnbc_sample_metadata.tsv")
cdr_file <- file.path(project_dir, "source_metadata", "TCGA_CDR", "TCGA-CDR-SupplementalTableS1.xlsx")
result_dir <- file.path(project_dir, "results", "bulk")
log_dir <- file.path(project_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expression_file)) stop("Missing TCGA TNBC expression metadata: ", expression_file)
if (!file.exists(cdr_file)) stop("Missing TCGA-CDR source file: ", cdr_file)
expected_cdr_md5 <- "a4591b2dcee39591f59e5e25a6ce75fa"
observed_cdr_md5 <- unname(tools::md5sum(cdr_file))
if (!identical(observed_cdr_md5, expected_cdr_md5)) {
  stop("TCGA-CDR MD5 mismatch: expected ", expected_cdr_md5, ", observed ", observed_cdr_md5)
}

expr <- fread(expression_file)
required_expr <- c("sample_id", "IFI44L_expression", "IFI44L_z", "median_group")
if (!all(required_expr %in% names(expr))) stop("Expression metadata schema mismatch")
documented_dir <- file.path(project_dir, "source_metadata", "TCGA_documented_primary")
source_receipt_file <- argument_value("--source-receipt", file.path(documented_dir, "TCGA_documented_input_provenance.json"))
if (!file.exists(source_receipt_file)) stop("Missing documented expression-source receipt")
provenance <- fromJSON(source_receipt_file)
if (!identical(provenance$source_matrix_sha256, "263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746") ||
    !identical(provenance$units, "log2(RSEM normalized count + 1)")) stop("Documented expression-source or units mismatch")
matrix_file <- argument_value("--expression-matrix", NULL)
if (is.null(matrix_file)) {
  subset_file <- file.path(documented_dir, provenance$subset_file)
  full_file <- file.path(documented_dir, provenance$source_matrix_file)
  matrix_file <- if (file.exists(subset_file)) subset_file else full_file
}
if (!file.exists(matrix_file)) stop("Missing documented full expression matrix or exact subset")
matrix_sha256 <- digest(file = matrix_file, algo = "sha256")
if (!matrix_sha256 %in% c(provenance$source_matrix_sha256, provenance$subset_sha256)) stop("Documented expression matrix SHA-256 mismatch")
manifest <- fread(file.path(result_dir, "TCGA_TNBC_sample_definition_manifest.tsv"))
included <- tolower(as.character(manifest$included_in_analysis)) %in% c("true", "1")
if (!setequal(expr$sample_id, manifest$sample_id[included])) stop("Survival rows differ from the documented clinical-definition manifest")
matrix_header <- names(fread(matrix_file, nrows = 0, check.names = FALSE))
source_values <- fread(matrix_file, select = c(matrix_header[[1]], expr$sample_id), check.names = FALSE)
source_ifi <- source_values[as.character(source_values[[1]]) == "IFI44L"]
if (nrow(source_ifi) != 1L) stop("Documented matrix must contain exactly one IFI44L row")
source_sample_order <- expr$sample_id
expected_ifi <- as.numeric(source_ifi[, ..source_sample_order])
if (!isTRUE(all.equal(expr$IFI44L_expression, expected_ifi, tolerance = 1e-10)) ||
    !isTRUE(all.equal(expr$IFI44L_z, as.numeric(scale(expected_ifi)), tolerance = 1e-10))) stop("Survival IFI44L values do not match the documented expression matrix")
expected_group <- ifelse(expected_ifi > median(expected_ifi), "High", "Low")
if (!identical(as.character(expr$median_group), expected_group)) stop("Median groups do not match documented IFI44L values")
rm(source_values, source_ifi)
expr[, patient_id := substr(sample_id, 1, 12)]
if (anyDuplicated(expr$patient_id)) stop("Multiple TNBC aliquots map to the same TCGA patient")

cdr <- as.data.table(suppressWarnings(read_excel(cdr_file, sheet = "TCGA-CDR")))
required_cdr <- c(
  "bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis",
  "ajcc_pathologic_tumor_stage", "OS", "OS.time", "PFI", "PFI.time"
)
if (!all(required_cdr %in% names(cdr))) stop("TCGA-CDR schema mismatch")
cdr <- cdr[type == "BRCA"]

merged <- merge(
  expr,
  cdr[, ..required_cdr],
  by.x = "patient_id",
  by.y = "bcr_patient_barcode",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(merged) != nrow(expr)) stop("Patient merge changed the TNBC row count")
if (anyNA(merged$type)) stop("At least one receptor-defined TNBC patient was absent from TCGA-CDR")

merged[, age_years := as.numeric(age_at_initial_pathologic_diagnosis)]
merged[, stage_group := fifelse(
  grepl("^Stage I($|A|B)", ajcc_pathologic_tumor_stage), "I",
  fifelse(
    grepl("^Stage II($|A|B)", ajcc_pathologic_tumor_stage), "II",
    fifelse(grepl("^Stage III($|A|B|C)", ajcc_pathologic_tumor_stage), "III", NA_character_)
  )
)]
merged[, stage_ordinal := fifelse(stage_group == "I", 1,
                                  fifelse(stage_group == "II", 2,
                                          fifelse(stage_group == "III", 3, NA_real_)))]
merged[, age_z := as.numeric(scale(age_years))]
merged[, `:=`(
  OS_event = as.integer(as.numeric(OS)),
  OS_days = as.numeric(`OS.time`),
  PFI_event = as.integer(as.numeric(PFI)),
  PFI_days = as.numeric(`PFI.time`)
)]

if (!all(merged$OS_event %in% c(0L, 1L)) || !all(merged$PFI_event %in% c(0L, 1L))) {
  stop("Unexpected event coding in TCGA-CDR")
}
if (any(merged$OS_days <= 0) || any(merged$PFI_days <= 0)) stop("Non-positive survival time detected")

cox_row <- function(data, endpoint, model_label, formula) {
  fit <- coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE)
  term <- "IFI44L_z"
  coef_table <- summary(fit)$coefficients
  conf_table <- summary(fit)$conf.int
  if (!term %in% rownames(coef_table)) stop("IFI44L term missing from Cox model")
  zph <- cox.zph(fit, transform = "km")
  zph_table <- zph$table
  data.table(
    cohort = "TCGA receptor-defined TNBC",
    endpoint = endpoint,
    model = model_label,
    n = fit$n,
    events = fit$nevent,
    hazard_ratio_per_IFI44L_SD = unname(conf_table[term, "exp(coef)"]),
    ci95_low = unname(conf_table[term, "lower .95"]),
    ci95_high = unname(conf_table[term, "upper .95"]),
    p_value = unname(coef_table[term, "Pr(>|z|)"]),
    IFI44L_PH_test_p = unname(zph_table[term, "p"]),
    global_PH_test_p = unname(zph_table["GLOBAL", "p"]),
    ties = "Efron",
    estimand = "relative hazard for the specified endpoint per 1-SD higher IFI44L expression"
  )
}

cox_rows <- list()
for (endpoint in c("OS", "PFI")) {
  time_col <- paste0(endpoint, "_days")
  event_col <- paste0(endpoint, "_event")
  unadjusted <- merged[is.finite(get(time_col)) & !is.na(get(event_col))]
  adjusted <- unadjusted[is.finite(age_z)]
  unadjusted_formula <- as.formula(sprintf("Surv(%s, %s) ~ IFI44L_z", time_col, event_col))
  adjusted_formula <- as.formula(sprintf("Surv(%s, %s) ~ IFI44L_z + age_z", time_col, event_col))
  cox_rows[[length(cox_rows) + 1L]] <- cox_row(
    unadjusted, endpoint, "Unadjusted", unadjusted_formula
  )
  cox_rows[[length(cox_rows) + 1L]] <- cox_row(
    adjusted, endpoint, "Age-adjusted", adjusted_formula
  )
}
cox_results <- rbindlist(cox_rows)

logrank_rows <- rbindlist(lapply(c("OS", "PFI"), function(endpoint) {
  time_col <- paste0(endpoint, "_days")
  event_col <- paste0(endpoint, "_event")
  dat <- merged[is.finite(get(time_col)) & !is.na(get(event_col))]
  formula <- as.formula(sprintf("Surv(%s, %s) ~ median_group", time_col, event_col))
  test <- survdiff(formula, data = dat, rho = 0)
  data.table(
    endpoint = endpoint,
    n = nrow(dat),
    events = sum(dat[[event_col]]),
    IFI44L_high_n = sum(dat$median_group == "High"),
    IFI44L_low_n = sum(dat$median_group == "Low"),
    logrank_chisq = unname(test$chisq),
    logrank_df = length(test$n) - 1L,
    logrank_p = pchisq(test$chisq, df = length(test$n) - 1L, lower.tail = FALSE),
    role = "graphical median-split sensitivity analysis; continuous Cox model is primary"
  )
}))

audit <- rbindlist(list(
  data.table(metric = "TCGA_TNBC_expression_rows", value = nrow(expr)),
  data.table(metric = "TCGA_CDR_MD5", value = observed_cdr_md5),
  data.table(metric = "documented_expression_source_sha256", value = provenance$source_matrix_sha256),
  data.table(metric = "analysed_expression_matrix_sha256", value = matrix_sha256),
  data.table(metric = "expression_units", value = provenance$units),
  data.table(metric = "unique_TCGA_patients", value = uniqueN(expr$patient_id)),
  data.table(metric = "matched_to_TCGA_CDR", value = sum(!is.na(merged$type))),
  data.table(metric = "OS_evaluable_n", value = sum(is.finite(merged$OS_days) & !is.na(merged$OS_event))),
  data.table(metric = "OS_events", value = sum(merged$OS_event, na.rm = TRUE)),
  data.table(metric = "PFI_evaluable_n", value = sum(is.finite(merged$PFI_days) & !is.na(merged$PFI_event))),
  data.table(metric = "PFI_events", value = sum(merged$PFI_event, na.rm = TRUE)),
  data.table(metric = "stage_I_n", value = sum(merged$stage_group == "I", na.rm = TRUE)),
  data.table(metric = "stage_II_n", value = sum(merged$stage_group == "II", na.rm = TRUE)),
  data.table(metric = "stage_III_n", value = sum(merged$stage_group == "III", na.rm = TRUE)),
  data.table(metric = "stage_unavailable_or_nonstandard_n", value = sum(is.na(merged$stage_group))),
  data.table(metric = "OS_events_stage_I", value = sum(merged[stage_group == "I"]$OS_event, na.rm = TRUE)),
  data.table(metric = "PFI_events_stage_I", value = sum(merged[stage_group == "I"]$PFI_event, na.rm = TRUE)),
  data.table(metric = "continuous_primary_predictor", value = "standardized IFI44L expression"),
  data.table(metric = "adjusted_model_status", value = "age-adjusted sensitivity only"),
  data.table(metric = "stage_adjustment_status", value = "not fitted; exploratory age-adjusted sensitivity only because overall and stage-specific event counts were limited"),
  data.table(metric = "patient_level_unit_of_inference", value = "true")
), fill = TRUE)

public_columns <- c(
  "sample_id", "patient_id", "IFI44L_expression", "IFI44L_z", "median_group",
  "age_years", "ajcc_pathologic_tumor_stage", "stage_group", "stage_ordinal",
  "OS_event", "OS_days", "PFI_event", "PFI_days"
)
setorder(merged, patient_id)
fwrite(merged[, ..public_columns], file.path(result_dir, "TCGA_TNBC_survival_merged.tsv"), sep = "\t")
fwrite(cox_results, file.path(result_dir, "TCGA_TNBC_continuous_IFI44L_Cox.tsv"), sep = "\t")
fwrite(logrank_rows, file.path(result_dir, "TCGA_TNBC_median_split_logrank_sensitivity.tsv"), sep = "\t")
fwrite(audit, file.path(result_dir, "TCGA_TNBC_survival_audit.tsv"), sep = "\t")
model_definitions <- copy(cox_results[, .(cohort, endpoint, model, n, events, ties, estimand)])
model_definitions[, `:=`(
  predictor = "documented log2(RSEM normalized count + 1) IFI44L values standardized using sample SD within the clinically verified cohort",
  adjustment = ifelse(model == "Age-adjusted", "age standardized using sample SD within cohort", "none"),
  endpoint_source = "TCGA Pan-Cancer Clinical Data Resource",
  outcome_merge_key = "first 12 characters of TCGA expression barcode = bcr_patient_barcode",
  median_split_role = "graphical sensitivity only; High strictly above median, Low at or below median"
)]
fwrite(model_definitions, file.path(result_dir, "TCGA_TNBC_Cox_model_definitions.tsv"), sep = "\t")
stage_events <- merged[, .(n = .N, OS_events = sum(OS_event), PFI_events = sum(PFI_event)), by = stage_group]
fwrite(stage_events, file.path(result_dir, "TCGA_TNBC_stage_event_counts.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "TCGA_TNBC_survival_R_sessionInfo.txt"))
print(cox_results)
print(logrank_rows)

