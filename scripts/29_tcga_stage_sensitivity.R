#!/usr/bin/env Rscript
# Added after the documented-source cohort update. Exploratory sensitivity,
# retaining the original continuous/age-adjusted core models without alteration.
suppressPackageStartupMessages({library(data.table);library(survival);library(digest);library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i+1L]]
}
project <- normalizePath(argument_value("--project", "."), mustWork = TRUE)
directory <- file.path(project, "results", "bulk")
input <- file.path(directory, "TCGA_TNBC_survival_merged.tsv")
provenance_file <- file.path(project, "source_metadata", "TCGA_documented_primary", "TCGA_documented_input_provenance.json")
provenance <- fromJSON(provenance_file)
if (!identical(provenance$source_matrix_sha256, "263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746")) stop("Stage sensitivity requires the documented TCGA source")
data <- fread(input)
required <- c("patient_id", "sample_id", "IFI44L_expression", "IFI44L_z", "age_years", "stage_group", "OS_event", "OS_days", "PFI_event", "PFI_days")
if (!all(required %in% names(data)) || anyDuplicated(data$patient_id)) stop("Merged outcome input schema or patient identity mismatch")
manifest <- fread(file.path(directory, "TCGA_TNBC_sample_definition_manifest.tsv"))
included <- tolower(as.character(manifest$included_in_analysis)) %in% c("true", "1")
if (!setequal(manifest$sample_id[included], data$sample_id)) stop("Stage sensitivity input differs from the selected primary cohort")
if (nrow(data) != provenance$primary_cohort_n) stop("Stage sensitivity is specified for primary TNBC only, not the strict sparse subset")
if (!isTRUE(all.equal(data$IFI44L_z, as.numeric(scale(data$IFI44L_expression)), tolerance = 1e-10))) stop("IFI44L z differs from full-primary-cohort standardization")
data[, age_per10_years := age_years / 10]
data[, stage_category := factor(stage_group, levels = c("I", "II", "III"))]
role <- "post-data-update exploratory age-and-stage sensitivity; sparse events; no variable selection; core Cox models retained"
results <- list();all_coefficients <- list();PH <- list();warnings_table <- list();definitions <- list();exclusions <- list();stage_counts <- list()
capture_warning <- function(endpoint, step) {
  function(w) warnings_table[[length(warnings_table)+1L]] <<- data.table(endpoint=endpoint,step=step,warning=conditionMessage(w))
}
for (endpoint in c("OS", "PFI")) {
  time <- paste0(endpoint, "_days");event <- paste0(endpoint, "_event")
  eligible <- is.finite(data$IFI44L_z) & is.finite(data$age_per10_years) & !is.na(data$stage_category) &
    is.finite(data[[time]]) & data[[time]] > 0 & data[[event]] %in% c(0,1)
  cases <- data[eligible]
  dropped <- copy(data[!eligible, .(patient_id,sample_id,age_years,stage_group)])
  dropped[, `:=`(endpoint=endpoint,reason="missing/nonstandard stage, non-finite covariate or unevaluable endpoint; complete-case rule")]
  exclusions[[endpoint]] <- dropped
  if (!all(c("I","II","III") %in% as.character(cases$stage_category))) stop("Stage categories incomplete for exploratory sensitivity")
  for (model_kind in c("categorical_stage", "stage_stratified")) {
  stratified <- model_kind == "stage_stratified"
  model_label <- if (stratified) "Age-adjusted stage-stratified exploratory sensitivity" else "Age-and-categorical-stage-adjusted exploratory sensitivity"
  model_role <- if (stratified) "diagnostic-driven exploratory addition after categorical-stage OS PH violation; applied to both endpoints; no significance-based model selection" else role
  key <- paste(endpoint,model_kind,sep="_")
  stage_formula <- if (stratified) "strata(stage_category)" else "stage_category"
  formula <- as.formula(sprintf("Surv(%s, %s) ~ IFI44L_z + age_per10_years + %s", time,event,stage_formula))
  fit <- withCallingHandlers(coxph(formula,data=cases,ties="efron",x=TRUE,y=TRUE),warning=capture_warning(endpoint,paste(model_kind,"coxph")))
  zph <- withCallingHandlers(cox.zph(fit,transform="km"),warning=capture_warning(endpoint,paste(model_kind,"cox.zph")))
  coefficients <- summary(fit)$coefficients;CI <- summary(fit)$conf.int
  expected_coefficients <- if (stratified) 2L else 4L
  if (nrow(coefficients)!=expected_coefficients) stop("Unexpected stage sensitivity coefficient count")
  all_coefficients[[key]] <- data.table(endpoint=endpoint,model=model_kind,term=rownames(coefficients),
    beta=coefficients[,"coef"],SE=coefficients[,"se(coef)"],wald_z=coefficients[,"z"],
    hazard_ratio=CI[,"exp(coef)"],ci95_low=CI[,"lower .95"],ci95_high=CI[,"upper .95"],p_value=coefficients[,"Pr(>|z|)"],n=fit$n,events=fit$nevent)
  results[[key]] <- data.table(cohort="documented-source TCGA primary TNBC",endpoint=endpoint,
    model=model_label,model_kind=model_kind,source_cohort_n=nrow(data),n=fit$n,excluded_n=nrow(data)-fit$n,events=fit$nevent,
    hazard_ratio_per_IFI44L_SD=CI["IFI44L_z","exp(coef)"],ci95_low=CI["IFI44L_z","lower .95"],ci95_high=CI["IFI44L_z","upper .95"],p_value=coefficients["IFI44L_z","Pr(>|z|)"],
    IFI44L_PH_test_p=zph$table["IFI44L_z","p"],global_PH_test_p=zph$table["GLOBAL","p"],fitted_coefficients=nrow(coefficients),ties="Efron",analysis_role=model_role)
  PH[[key]] <- data.table(endpoint=endpoint,model=model_kind,term=rownames(zph$table),chisq=zph$table[,"chisq"],df=zph$table[,"df"],p_value=zph$table[,"p"])
  definitions[[key]] <- data.table(endpoint=endpoint,model=model_kind,formula=paste(deparse(formula),collapse=" "),n=fit$n,events=fit$nevent,
    IFI44L_scale="one sample SD in the full 157-patient documented primary cohort, retained from core models",
    age_scale="per 10 years; no data-driven cutpoint",stage=if (stratified) "separate baseline hazards for I/II/III; no common stage hazard-ratio coefficient" else "categorical I/II/III; stage I reference; unknown/nonstandard stage excluded",
    complete_case_rule="finite IFI44L/age; I/II/III stage; positive endpoint time and event in {0,1}",
    inference="nominal two-sided Wald P and 95% CI; sparse exploratory sensitivity, not independent prognostic validation",
    proportional_hazards="Schoenfeld residual tests, km transform; term/global tests retained; global non-significance does not establish every term's PH assumption",analysis_role=model_role)
  }
  counts <- cases[, .(n=.N,events=sum(get(event))), by=.(stage=as.character(stage_category))]
  counts[, endpoint := endpoint]
  stage_counts[[endpoint]] <- counts
}
if (!length(warnings_table)) warnings_table <- list(data.table(endpoint=character(),step=character(),warning=character()))
files <- list(TCGA_TNBC_stage_adjusted_IFI44L_Cox=rbindlist(results),TCGA_TNBC_stage_sensitivity_all_coefficients=rbindlist(all_coefficients),
  TCGA_TNBC_stage_sensitivity_PH_diagnostics=rbindlist(PH),TCGA_TNBC_stage_sensitivity_warnings=rbindlist(warnings_table),
  TCGA_TNBC_stage_sensitivity_model_definitions=rbindlist(definitions),TCGA_TNBC_stage_sensitivity_exclusions=rbindlist(exclusions),TCGA_TNBC_stage_sensitivity_stage_events=rbindlist(stage_counts))
for (name in names(files)) fwrite(files[[name]],file.path(directory,paste0(name,".tsv")),sep="\t")
receipt <- list(analysis_role=role,input_file=basename(input),input_sha256=digest(file=input,algo="sha256"),
  documented_source_sha256=provenance$source_matrix_sha256,source_cohort_n=nrow(data),models=length(results),outcomes=c("OS","PFI"),
  stage_reference="I for categorical model; independent baseline hazards for stratified model",fitted_coefficient_n=list(categorical_stage=4L,stage_stratified=2L),
  stage_stratified_analysis_role="diagnostic-driven addition after categorical-stage OS PH violation, uniformly applied to OS/PFI; all model results retained",
  core_models_modified=FALSE,strict_cohort_model_fitted=FALSE,warnings_n=nrow(rbindlist(warnings_table)))
write_json(receipt,file.path(directory,"TCGA_TNBC_stage_sensitivity_provenance.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(directory,"TCGA_TNBC_stage_sensitivity_R_sessionInfo.txt"))
print(rbindlist(results))
