#!/usr/bin/env Rscript
# TCGA-result-driven exploratory external consistency analysis. This analysis
# is not pre-registered and does not estimate an independent causal effect.
suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(digest)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}
required_argument <- function(flag) {
  value <- argument_value(flag)
  if (is.null(value)) stop("Required argument: ", flag)
  value
}
cohort <- required_argument("--cohort")
expression_file <- normalizePath(required_argument("--expression-matrix"), mustWork = TRUE)
metadata_file <- normalizePath(required_argument("--sample-metadata"), mustWork = TRUE)
definition_file <- normalizePath(required_argument("--signature-definitions"), mustWork = TRUE)
source_scale <- required_argument("--source-scale")
output_dir <- required_argument("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260831)
checkpoints <- c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TNFRSF4", "TNFRSF9")
covariates <- c("IFN_ALPHA_RESPONSE", "IFN_GAMMA_RESPONSE", "T_CELL_INFLAMED", "EPITHELIAL_CONTENT")
classification <- "TCGA-result-driven exploratory external consistency analysis; not pre-registered"

metadata <- fread(metadata_file)
if (!"sample_id" %in% names(metadata) || anyDuplicated(metadata$sample_id)) stop("Metadata requires unique sample_id")
values <- fread(expression_file, check.names = FALSE)
gene <- as.character(values[[1L]])
if (anyNA(gene) || any(!nzchar(gene)) || anyDuplicated(gene)) stop("Expression matrix requires unique nonempty gene symbols")
if (!setequal(names(values)[-1L], metadata$sample_id)) stop("Expression columns and metadata sample universe differ")
sample_order <- as.character(metadata$sample_id)
expression <- as.matrix(values[, ..sample_order]);storage.mode(expression) <- "double";rownames(expression) <- gene
input_gene_rows <- nrow(expression)
if (any(!is.finite(expression))) stop("External expression matrix contains non-finite values")
expression <- expression[apply(expression, 1L, sd) > 0, , drop = FALSE]
if (!"IFI44L" %in% rownames(expression)) stop("IFI44L is absent or constant")
ifi <- as.numeric(expression["IFI44L", ])
metadata_ifi_column <- argument_value("--ifi-column")
if (is.null(metadata_ifi_column)) {
  matches <- intersect(c("IFI44L_log2", "IFI44L_expression"), names(metadata))
  if (length(matches) != 1L) stop("Specify --ifi-column when metadata has no unique IFI44L expression field")
  metadata_ifi_column <- matches[[1L]]
}
if (!metadata_ifi_column %in% names(metadata)) stop("IFI44L metadata field absent")
if (!isTRUE(all.equal(as.numeric(metadata[[metadata_ifi_column]]), ifi, tolerance = 1e-10))) stop("External IFI44L values differ from the clinical-analysis metadata")

definitions <- fread(definition_file)
if (!all(c("signature", "gene") %in% names(definitions))) stop("Signature definition schema mismatch")
definitions <- unique(definitions[, .(signature, gene)])
if (!all(covariates %in% definitions$signature)) stop("Missing one of the four predefined modules")
# Intentionally ignore the TCGA-specific present column. Each external matrix
# supplies its own platform coverage and its own full transcriptomic background.
signatures <- lapply(covariates, function(s) unique(definitions[signature == s, gene]))
names(signatures) <- covariates
partial_rank <- function(x, y, covs) {
  matrix_values <- cbind(x = x, y = y, as.matrix(covs))
  if (any(!is.finite(matrix_values))) stop("Non-finite partial-rank input")
  ranks <- apply(matrix_values, 2L, rank, ties.method = "average")
  design <- cbind(Intercept = 1, ranks[, -(1:2), drop = FALSE])
  decomposition <- qr(design)
  xr <- qr.resid(decomposition, ranks[, "x"])
  yr <- qr.resid(decomposition, ranks[, "y"])
  rho <- unname(cor(xr, yr))
  df <- nrow(matrix_values) - decomposition$rank - 1L
  if (!is.finite(rho) || df <= 0L) stop("Non-estimable partial-rank model")
  statistic <- rho * sqrt(df / max(1 - rho^2, .Machine$double.eps))
  list(rho = rho, p = 2 * pt(abs(statistic), df = df, lower.tail = FALSE),
       n = nrow(matrix_values), k = decomposition$rank - 1L, df = df)
}
rows <- list();coverage <- list();models <- list();inputs <- list();warnings_reported <- character()
for (target in checkpoints) {
  model_sets <- lapply(signatures, function(members) intersect(setdiff(members, c("IFI44L", target)), rownames(expression)))
  coverage[[target]] <- rbindlist(lapply(covariates, function(s) {
    members <- signatures[[s]]
    data.table(cohort = cohort, checkpoint = target, signature = s, gene = members,
               present_in_external_matrix = members %in% gene,
               retained_variable_gene = members %in% rownames(expression),
               excluded_target_or_IFI44L = members %in% c("IFI44L", target),
               used_for_adjustment = members %in% model_sets[[s]])
  }))
  if (any(lengths(model_sets) < 2L)) stop("Insufficient target-excluded signature coverage for ", target)
  if (!target %in% rownames(expression)) {
    rows[[target]] <- data.table(cohort = cohort, gene = target, status = "unavailable or constant checkpoint",
      raw_spearman_rho = NA_real_, raw_p_value = NA_real_, adjusted_spearman_rho = NA_real_,
      adjusted_p_value = NA_real_, n = ncol(expression), covariate_rank = NA_integer_, degrees_of_freedom = NA_integer_)
    next
  }
  scores <- withCallingHandlers(gsva(ssgseaParam(expression, model_sets, minSize = 2, maxSize = 500, normalize = TRUE), verbose = FALSE),
    warning = function(w) warnings_reported <<- c(warnings_reported, paste(target, conditionMessage(w))))
  covs <- as.data.frame(t(scores[covariates, , drop = FALSE]))
  y <- as.numeric(expression[target, ])
  result <- partial_rank(ifi, y, covs)
  raw <- suppressWarnings(cor.test(ifi, y, method = "spearman", exact = FALSE))
  rows[[target]] <- data.table(cohort = cohort, gene = target, status = "estimable",
    raw_spearman_rho = unname(raw$estimate), raw_p_value = raw$p.value,
    adjusted_spearman_rho = result$rho, adjusted_p_value = result$p,
    n = result$n, covariate_rank = result$k, degrees_of_freedom = result$df)
  inputs[[target]] <- cbind(data.table(cohort = cohort, checkpoint = target, sample_id = sample_order,
    IFI44L_expression = ifi, checkpoint_expression = y), as.data.table(covs))
  models[[target]] <- data.table(cohort = cohort, checkpoint = target, n = result$n,
    covariates = paste(covariates, collapse = ";"),
    model = "rank all variables; OLS residualization; Pearson correlation of rank residuals; approximate two-sided t reference",
    degrees_of_freedom = result$df, covariate_rank = result$k,
    adjustment_signature_policy = "full predefined lists; platform-specific coverage; exclude IFI44L and current checkpoint from all four modules",
    signature_scoring = "ssGSEA on full external normalized gene matrix, normalize=TRUE; no additional expression transformation",
    multiplicity = "BH across all nine checkpoints separately within this external cohort; unavailable p-values stay NA and remain in family",
    analysis_role = classification)
}
results <- rbindlist(rows, fill = TRUE)
results[, `:=`(raw_fdr_bh = NA_real_, adjusted_fdr_bh = NA_real_)]
valid <- is.finite(results$raw_p_value)
results[valid, raw_fdr_bh := p.adjust(raw_p_value, method = "BH", n = length(checkpoints))]
valid <- is.finite(results$adjusted_p_value)
results[valid, adjusted_fdr_bh := p.adjust(adjusted_p_value, method = "BH", n = length(checkpoints))]
results[, `:=`(analysis_role = classification, BH_family_size = length(checkpoints))]
coverage_all <- rbindlist(coverage)
coverage_counts <- coverage_all[, .(predefined_genes = .N, present_in_external_matrix = sum(present_in_external_matrix),
  variable_genes = sum(retained_variable_gene), used_for_adjustment = sum(used_for_adjustment)), by = .(cohort, checkpoint, signature)]
fwrite(results, file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_partial_correlations.tsv")), sep = "\t")
fwrite(rbindlist(models), file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_model_definitions.tsv")), sep = "\t")
fwrite(coverage_all, file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_gene_coverage.tsv")), sep = "\t")
fwrite(coverage_counts, file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_coverage_counts.tsv")), sep = "\t")
fwrite(rbindlist(inputs), file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_model_inputs.tsv")), sep = "\t")
writeLines(if (length(warnings_reported)) warnings_reported else "No warnings reported", file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_warnings.txt")))
receipt <- list(cohort = cohort, analysis_role = classification, source_scale = source_scale, additional_transformation = "none",
  expression_file = basename(expression_file), expression_sha256 = digest(file = expression_file, algo = "sha256"),
  sample_metadata_file = basename(metadata_file), sample_metadata_sha256 = digest(file = metadata_file, algo = "sha256"),
  signature_definitions_file = basename(definition_file), signature_definitions_sha256 = digest(file = definition_file, algo = "sha256"),
  sample_n = ncol(expression), input_gene_rows = input_gene_rows, retained_variable_gene_rows = nrow(expression),
  IFI44L_metadata_column = metadata_ifi_column, checkpoint_family_size = length(checkpoints),
  estimable_checkpoint_tests = sum(results$status == "estimable"), covariate_module_n = length(covariates),
  outcomes_used = FALSE, warnings_n = length(warnings_reported), sample_order = "clinical-analysis metadata sample_id order")
write_json(receipt, file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_provenance.json")), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, paste0(cohort, "_checkpoint_exploratory_R_sessionInfo.txt")))
print(results)
