#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(limma)
  library(GSVA)
  library(survival)
  library(survminer)
  library(data.table)
  library(ggplot2)
})

command_args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(flag, default = NULL) {
  position <- match(flag, command_args)
  if (is.na(position)) return(default)
  if (position == length(command_args)) stop("Missing value after ", flag)
  command_args[[position + 1]]
}

project_dir <- normalizePath(argument_value("--project", "."), mustWork = TRUE)
source_file <- file.path(project_dir, "source_metadata", "GSE58812", "GSE58812_series_matrix.txt.gz")
result_dir <- file.path(project_dir, "results", "validation_GSE58812")
figure_dir <- file.path(project_dir, "figures", "validation_GSE58812")
log_dir <- file.path(project_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) stop("Missing downloaded GSE58812 series matrix: ", source_file)
gse <- getGEO(filename = source_file, getGPL = FALSE)
probe_expr <- exprs(gse)
clinical <- as.data.table(pData(gse), keep.rownames = "sample_id")
fwrite(data.table(
  original_field = names(clinical),
  nonmissing_n = vapply(clinical, function(x) sum(!is.na(x) & nzchar(as.character(x))), integer(1)),
  source = "GSE58812 GEO series matrix pData"
), file.path(result_dir, "GSE58812_clinical_field_availability.tsv"), sep = "\t")

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
clinical[, `:=`(
  age_years = num(`age at diag:ch1`),
  metastasis_event = num(`meta:ch1`),
  mfs_days = num(`mfs (days):ch1`),
  death_event = num(`death:ch1`),
  os_days = num(`os (days):ch1`)
)]
if (!all(clinical$metastasis_event %in% c(0, 1)) ||
    !all(clinical$death_event %in% c(0, 1))) stop("Unexpected GSE58812 survival-event coding")
if (any(!is.finite(clinical$mfs_days) | clinical$mfs_days <= 0) ||
    any(!is.finite(clinical$os_days) | clinical$os_days <= 0)) stop("Invalid GSE58812 survival time")

# MAS5 values are on a linear intensity scale in this series matrix.
expr_log2 <- log2(pmax(probe_expr, 0) + 1)
probe_annotation <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr_log2),
  keytype = "PROBEID",
  columns = c("PROBEID", "SYMBOL", "ENTREZID")
)
probe_annotation <- unique(probe_annotation[!is.na(probe_annotation$SYMBOL) & nzchar(probe_annotation$SYMBOL), ])
mapped <- expr_log2[match(probe_annotation$PROBEID, rownames(expr_log2)), , drop = FALSE]
rownames(mapped) <- probe_annotation$SYMBOL
gene_expr <- avereps(mapped, ID = rownames(mapped))
if (!"IFI44L" %in% rownames(gene_expr)) stop("IFI44L mapping failed")
fwrite(as.data.table(gene_expr, keep.rownames = "gene"),
       file.path(result_dir, "GSE58812_normalized_gene_expression.tsv.gz"), sep = "\t", compress = "gzip")

ifi <- as.numeric(gene_expr["IFI44L", clinical$sample_id])
clinical[, `:=`(
  IFI44L_log2 = ifi,
  IFI44L_z = as.numeric(scale(ifi)),
  age_z = as.numeric(scale(age_years)),
  IFI44L_group = factor(ifelse(ifi >= median(ifi), "High", "Low"), levels = c("Low", "High"))
)]
fwrite(clinical[, .(sample_id, age_years, metastasis_event, mfs_days, death_event, os_days,
                    IFI44L_log2, IFI44L_z, IFI44L_group)],
       file.path(result_dir, "GSE58812_clinical_IFI44L.tsv"), sep = "\t")

fit_cox <- function(time, event, adjusted = FALSE) {
  vars <- c(time, event, "IFI44L_z", if (adjusted) "age_z")
  dat <- na.omit(as.data.frame(clinical[, ..vars]))
  formula <- if (adjusted) {
    as.formula(sprintf("Surv(%s, %s) ~ IFI44L_z + age_z", time, event))
  } else {
    as.formula(sprintf("Surv(%s, %s) ~ IFI44L_z", time, event))
  }
  model <- coxph(formula, data = dat, ties = "efron", x = TRUE)
  s <- summary(model)
  row <- which(rownames(s$coefficients) == "IFI44L_z")
  zph <- cox.zph(model)
  data.table(
    endpoint = ifelse(time == "mfs_days", "Metastasis-free survival", "Overall survival"),
    model = ifelse(adjusted, "Age-adjusted", "Unadjusted"),
    n = model$n,
    events = model$nevent,
    hazard_ratio_per_IFI44L_SD = s$conf.int[row, "exp(coef)"],
    ci95_low = s$conf.int[row, "lower .95"],
    ci95_high = s$conf.int[row, "upper .95"],
    p_value = s$coefficients[row, "Pr(>|z|)"],
    proportional_hazards_p = zph$table["IFI44L_z", "p"]
  )
}

cox_results <- rbindlist(list(
  fit_cox("mfs_days", "metastasis_event", FALSE),
  fit_cox("mfs_days", "metastasis_event", TRUE),
  fit_cox("os_days", "death_event", FALSE),
  fit_cox("os_days", "death_event", TRUE)
))
cox_results[, fdr_bh := p.adjust(p_value, method = "BH")]
fwrite(cox_results, file.path(result_dir, "GSE58812_continuous_IFI44L_Cox.tsv"), sep = "\t")
model_definitions <- copy(cox_results[, .(endpoint, model, n, events)])
model_definitions[, `:=`(
  predictor = "IFI44L gene-level log2(MAS5+1) expression standardized using sample SD across 107 patients",
  adjustment = ifelse(model == "Age-adjusted", "patient age standardized using sample SD", "none"),
  tie_method = "Efron", PH_diagnostic = "Schoenfeld residual test (survival::cox.zph)",
  group_rule = "median split for graphical sensitivity only; High >= median, Low < median",
  unavailable_covariates = "patient-level tumour size; nodal status; TNM stage; histological grade absent from deposited series matrix"
)]
fwrite(model_definitions, file.path(result_dir, "GSE58812_Cox_model_definitions.tsv"), sep = "\t")

# Reuse the predeclared signatures, excluding IFI44L from every score to prevent
# circular self-correlation.
sig_def <- fread(file.path(project_dir, "results", "bulk", "immune_signature_definitions.tsv"))
signature_list <- split(sig_def[gene != "IFI44L", gene], sig_def[gene != "IFI44L", signature])
signature_list <- lapply(signature_list, intersect, y = rownames(gene_expr))
signature_list <- signature_list[lengths(signature_list) >= 2]
gse_signature_members <- copy(sig_def)
gse_signature_members[, `:=`(
  present_in_GSE58812 = gene %in% rownames(gene_expr),
  used_in_GSE58812 = gene != "IFI44L" & gene %in% rownames(gene_expr)
)]
fwrite(gse_signature_members,
       file.path(result_dir, "GSE58812_immune_signature_definitions.tsv"), sep = "\t")
gsva_param <- ssgseaParam(exprData = gene_expr, geneSets = signature_list, normalize = TRUE)
scores <- gsva(gsva_param, verbose = FALSE)
fwrite(as.data.table(t(scores), keep.rownames = "sample_id"),
       file.path(result_dir, "GSE58812_immune_signature_ssGSEA_scores.tsv"), sep = "\t")

cor_results <- rbindlist(lapply(rownames(scores), function(sig) {
  test <- suppressWarnings(cor.test(ifi, as.numeric(scores[sig, clinical$sample_id]), method = "spearman", exact = FALSE))
  data.table(signature = sig, n = length(ifi), spearman_rho = unname(test$estimate), p_value = test$p.value)
}))
cor_results[, fdr_bh := p.adjust(p_value, method = "BH")]
setorder(cor_results, -spearman_rho)
fwrite(cor_results, file.path(result_dir, "GSE58812_IFI44L_immune_signature_correlations.tsv"), sep = "\t")

checkpoint_genes <- intersect(c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TNFRSF4", "TNFRSF9"), rownames(gene_expr))
checkpoint_results <- rbindlist(lapply(checkpoint_genes, function(gene) {
  test <- suppressWarnings(cor.test(ifi, gene_expr[gene, clinical$sample_id], method = "spearman", exact = FALSE))
  data.table(gene = gene, n = length(ifi), spearman_rho = unname(test$estimate), p_value = test$p.value)
}))
checkpoint_results[, fdr_bh := p.adjust(p_value, method = "BH")]
setorder(checkpoint_results, -spearman_rho)
fwrite(checkpoint_results, file.path(result_dir, "GSE58812_IFI44L_checkpoint_correlations.tsv"), sep = "\t")

plot_km <- function(time, event, endpoint_label, filename) {
  dat <- as.data.frame(clinical[, .(time = get(time), event = get(event), IFI44L_group)])
  dat <- na.omit(dat)
  fit <- survfit(Surv(time, event) ~ IFI44L_group, data = dat)
  p <- ggsurvplot(
    fit, data = dat, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
    palette = c("#2166AC", "#B2182B"), xlab = "Time (days)", ylab = endpoint_label,
    legend.title = "IFI44L", legend.labs = c("Low", "High"), ggtheme = theme_classic(base_size = 11)
  )
  pdf(file.path(figure_dir, paste0(filename, ".pdf")), width = 6.5, height = 6.4)
  print(p)
  dev.off()
  png(file.path(figure_dir, paste0(filename, ".png")), width = 2600, height = 2500, res = 400)
  print(p)
  dev.off()
}
plot_km("mfs_days", "metastasis_event", "Metastasis-free survival probability", "GSE58812_IFI44L_MFS_KM_median_sensitivity")
plot_km("os_days", "death_event", "Overall survival probability", "GSE58812_IFI44L_OS_KM_median_sensitivity")

p_cor <- ggplot(cor_results, aes(reorder(signature, spearman_rho), spearman_rho, fill = fdr_bh < .05)) +
  geom_col(width = .72) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70"), guide = "none") +
  labs(x = NULL, y = "Spearman correlation with IFI44L", title = "Independent TNBC validation: GSE58812 (n=107)") +
  theme_classic(base_size = 11)
ggsave(file.path(figure_dir, "GSE58812_IFI44L_immune_correlations.pdf"), p_cor, width = 7.2, height = 5.8)
ggsave(file.path(figure_dir, "GSE58812_IFI44L_immune_correlations.png"), p_cor, width = 7.2, height = 5.8, dpi = 400)

audit <- data.table(
  dataset = "GSE58812", platform = annotation(gse), n_samples = ncol(gene_expr),
  n_genes = nrow(gene_expr), IFI44L_probe = "204439_at",
  preprocessing = "GEO series matrix MAS5 intensities; log2(x+1); probes collapsed by limma::avereps",
  primary_inference = "continuous IFI44L per SD; Cox model",
  median_split = "visual sensitivity analysis only"
)
fwrite(audit, file.path(result_dir, "GSE58812_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "GSE58812_R_sessionInfo.txt"))
message("GSE58812 validation completed successfully.")
