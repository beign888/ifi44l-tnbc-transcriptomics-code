#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(GSVA)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(digest)
  library(jsonlite)
})

set.seed(20260831)

command_args <- commandArgs(trailingOnly = TRUE)
argument_value <- function(flag, default = NULL) {
  position <- match(flag, command_args)
  if (is.na(position)) return(default)
  if (position == length(command_args)) stop("Missing value after ", flag)
  command_args[[position + 1]]
}

project_dir <- normalizePath(argument_value("--project", "."), mustWork = TRUE)
documented_dir <- file.path(project_dir, "source_metadata", "TCGA_documented_primary")
result_dir <- file.path(project_dir, "results", "bulk")
figure_dir <- file.path(project_dir, "figures", "bulk")
log_dir <- file.path(project_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

source_receipt_file <- argument_value("--source-receipt", file.path(documented_dir, "TCGA_documented_input_provenance.json"))
if (!file.exists(source_receipt_file)) stop("Missing documented-source receipt; run 24_prepare_documented_tcga_input.py or restore the supplied documented input bundle")
provenance <- fromJSON(source_receipt_file)
expected_source_sha256 <- "263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746"
if (!identical(provenance$source_matrix_sha256, expected_source_sha256) ||
    !identical(provenance$units, "log2(RSEM normalized count + 1)")) stop("Documented expression source or units mismatch")
header_file <- file.path(documented_dir, provenance$full_header_file)
if (!file.exists(header_file) || !identical(digest(file = header_file, algo = "sha256"), provenance$full_header_sha256)) stop("Documented full-header SHA-256 mismatch")
clinical_file <- file.path(project_dir, "source_metadata", "TCGA_clinical", "brca_tcga_clinical_patient.txt")
if (!file.exists(clinical_file) || !identical(digest(file = clinical_file, algo = "sha256"), provenance$clinical_source_sha256)) stop("Clinical source SHA-256 mismatch")
matrix_file <- argument_value("--expression-matrix", NULL)
if (is.null(matrix_file)) {
  subset_file <- file.path(documented_dir, provenance$subset_file)
  full_file <- file.path(documented_dir, provenance$source_matrix_file)
  matrix_file <- if (file.exists(subset_file)) subset_file else full_file
}
if (!file.exists(matrix_file)) stop("Missing documented expression matrix; archived matrices are not a fallback")
observed_matrix_sha256 <- digest(file = matrix_file, algo = "sha256")
if (!observed_matrix_sha256 %in% c(provenance$source_matrix_sha256, provenance$subset_sha256)) stop("Supplied matrix does not match the documented full matrix or exact subset")
manifest_file <- file.path(result_dir, "TCGA_TNBC_sample_definition_manifest.tsv")
if (!file.exists(manifest_file)) stop("Run 21_audit_tcga_tnbc_definition.py first to generate the clinically verified TNBC manifest")
manifest <- fread(manifest_file, colClasses = c(sample_type_code = "character"))
if (!all(c("sample_id", "included_in_analysis") %in% names(manifest))) stop("TNBC definition-manifest schema mismatch")
is_included <- tolower(as.character(manifest$included_in_analysis)) %in% c("true", "1")
tnbc_requested <- unique(as.character(manifest$sample_id[is_included]))
if (!length(tnbc_requested) || anyDuplicated(tnbc_requested)) stop("Empty or duplicated TNBC inclusion list")
if (!all(manifest$sample_type_code[is_included] == "01") ||
    anyDuplicated(manifest$patient_id[is_included])) stop("Analysis manifest must contain one type01 sample per patient")

message("[1/8] Reading header and matching the TNBC sample list")
header <- names(fread(matrix_file, nrows = 0, check.names = FALSE))
all_samples <- header[-1]
tnbc_samples <- tnbc_requested[tnbc_requested %in% all_samples]
tnbc_missing <- setdiff(tnbc_requested, all_samples)
if (length(tnbc_missing)) stop("The clinically eligible TNBC list contains samples absent from the supplied expression matrix")

audit <- data.table(
  metric = c(
    "matrix_samples", "tnbc_requested", "tnbc_matched", "tnbc_missing",
    "analysis_seed"
  ),
  value = c(length(all_samples), length(tnbc_requested), length(tnbc_samples),
            length(tnbc_missing), 20260831)
)
fwrite(audit, file.path(result_dir, "cohort_audit.tsv"), sep = "\t")
fwrite(data.table(sample_id = tnbc_missing),
       file.path(result_dir, "tnbc_samples_missing_from_matrix.tsv"), sep = "\t")

message("[2/8] Loading only the ", length(tnbc_samples), " clinically verified TNBC columns")
dt <- fread(matrix_file, select = c(header[1], tnbc_samples), check.names = FALSE)
gene_id <- as.character(dt[[1]])
expr <- as.matrix(dt[, -1, with = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- gene_id
rm(dt)
gc()

# Collapse duplicated gene symbols by their mean, matching limma's documented behavior.
expr <- avereps(expr, ID = rownames(expr))
expr <- expr[rowSums(is.finite(expr)) == ncol(expr), , drop = FALSE]
expr <- expr[apply(expr, 1, sd) > 0, , drop = FALSE]
stopifnot("IFI44L" %in% rownames(expr), ncol(expr) == length(tnbc_samples))
response_expr <- expr[rownames(expr) != "IFI44L", , drop = FALSE]
input_audit <- data.table(
  metric = c("source_matrix_basename", "analysed_samples", "retained_variable_genes",
             "response_genes_excluding_IFI44L", "documented_expression_min", "documented_expression_q25",
             "documented_expression_median", "documented_expression_q75", "documented_expression_max",
             "expression_units", "additional_transformation_applied", "upstream_release_and_normalization",
             "source_matrix_sha256", "analysed_matrix_sha256", "source_receipt_sha256", "full_source_columns"),
  value = c(basename(matrix_file), ncol(expr), nrow(expr), nrow(response_expr),
            min(expr), quantile(expr, .25), median(expr), quantile(expr, .75), max(expr),
            provenance$units, "none; supplied log2 RSEM values used directly; no additional log or count normalization",
            "UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2; 2017-10-13; UNC Level3 RSEM normalized_results",
            provenance$source_matrix_sha256, observed_matrix_sha256,
            digest(file = source_receipt_file, algo = "sha256"), provenance$full_header_samples)
)
fwrite(input_audit, file.path(result_dir, "TCGA_TNBC_expression_input_audit.tsv"), sep = "\t")

ifi <- as.numeric(expr["IFI44L", ])
names(ifi) <- colnames(expr)
ifi_z <- as.numeric(scale(ifi))
median_cut <- median(ifi)
group <- factor(ifelse(ifi > median_cut, "High", "Low"), levels = c("Low", "High"))

sample_metadata <- data.table(
  sample_id = colnames(expr),
  IFI44L_expression = ifi,
  IFI44L_z = ifi_z,
  median_group = as.character(group)
)
fwrite(sample_metadata, file.path(result_dir, "tnbc_sample_metadata.tsv"), sep = "\t")

message("[3/8] Continuous IFI44L association model with BH-FDR")
design_cont <- model.matrix(~ ifi_z)
fit_cont <- eBayes(lmFit(response_expr, design_cont), trend = TRUE, robust = TRUE)
cont <- as.data.table(topTable(fit_cont, coef = "ifi_z", number = Inf,
                              adjust.method = "BH", sort.by = "none"), keep.rownames = "gene")
setnames(cont, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
         c("beta_per_IFI44L_SD", "mean_expression", "moderated_t", "p_value", "fdr_bh", "lods"))
cont[, abs_moderated_t := abs(moderated_t)]
setorder(cont, fdr_bh, -abs_moderated_t)
cont[, abs_moderated_t := NULL]
fwrite(cont, file.path(result_dir, "continuous_IFI44L_gene_associations.tsv"), sep = "\t")
fwrite(cont[fdr_bh < 0.05],
       file.path(result_dir, "continuous_IFI44L_gene_associations_FDR05.tsv"), sep = "\t")

message("[4/8] Median high-versus-low sensitivity analysis with BH-FDR")
design_group <- model.matrix(~ 0 + group)
colnames(design_group) <- levels(group)
fit_group <- contrasts.fit(lmFit(response_expr, design_group), makeContrasts(High - Low, levels = design_group))
fit_group <- eBayes(fit_group, trend = TRUE, robust = TRUE)
de <- as.data.table(topTable(fit_group, number = Inf, adjust.method = "BH", sort.by = "none"),
                    keep.rownames = "gene")
setnames(de, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
         c("expression_difference_high_minus_low", "mean_expression", "moderated_t", "p_value", "fdr_bh", "lods"))
de[, abs_moderated_t := abs(moderated_t)]
setorder(de, fdr_bh, -abs_moderated_t)
de[, abs_moderated_t := NULL]
fwrite(de, file.path(result_dir, "median_group_differential_expression.tsv"), sep = "\t")

message("[5/8] GO Biological Process GSEA using the full continuous-model ranking")
rank_dt <- cont[gene != "IFI44L" & is.finite(moderated_t), .(gene, moderated_t)]
mapping <- suppressMessages(bitr(rank_dt$gene, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db))
rank_mapped <- merge(rank_dt, as.data.table(mapping), by.x = "gene", by.y = "SYMBOL")
rank_mapped <- rank_mapped[order(-abs(moderated_t))][!duplicated(ENTREZID)]
mapping_audit <- data.table(
  metric = c("ranked_symbols_excluding_IFI44L", "symbols_mapped_to_ENTREZID",
             "symbols_unmapped", "fraction_symbols_unmapped", "unique_ENTREZID_in_ranking",
             "multiple_mapping_resolution", "GO_gene_set_size_min", "GO_gene_set_size_max",
             "fgsea_nPermSimple", "random_seed"),
  value = c(nrow(rank_dt), uniqueN(mapping$SYMBOL),
            nrow(rank_dt) - uniqueN(mapping$SYMBOL),
            (nrow(rank_dt) - uniqueN(mapping$SYMBOL)) / nrow(rank_dt), nrow(rank_mapped),
            "retain largest absolute moderated t per ENTREZID", 15, 500, 10000, 20260831)
)
fwrite(mapping_audit, file.path(result_dir, "GSEA_GO_BP_mapping_audit.tsv"), sep = "\t")
gene_list <- rank_mapped$moderated_t
names(gene_list) <- rank_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_warnings <- character()
gsea_go <- withCallingHandlers(suppressMessages(gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  by = "fgsea",
  nPermSimple = 10000,
  verbose = FALSE,
  seed = TRUE
)), warning = function(w) {
  gsea_warnings <<- c(gsea_warnings, conditionMessage(w))
})
# DOSE/clusterProfiler filters unresolved p-values before returning gseaResult.
# Reconstruct the eligible test population from its stored gene sets and the
# exact ranking; retain omitted eligible IDs as explicitly unresolved rows.
go_tab <- as.data.table(gsea_go@result)
eligible_set_sizes <- vapply(gsea_go@geneSets, function(members) {
  length(intersect(members, names(gene_list)))
}, integer(1))
eligible_ids <- names(eligible_set_sizes)[eligible_set_sizes >= 15 & eligible_set_sizes <= 500]
unresolved_ids <- setdiff(eligible_ids, go_tab$ID)
go_tab[, result_status := "resolved"]
if (length(unresolved_ids)) {
  unresolved <- data.table(
    ID = unresolved_ids, Description = paste("Unresolved pathway", unresolved_ids),
    setSize = unname(eligible_set_sizes[unresolved_ids]),
    enrichmentScore = NA_real_, NES = NA_real_, pvalue = NA_real_, p.adjust = NA_real_,
    qvalue = NA_real_, rank = NA_integer_, leading_edge = NA_character_,
    core_enrichment = NA_character_,
    result_status = "eligible pathway omitted by fgsea/DOSE because p-value or NES was unresolved"
  )
  go_tab <- rbindlist(list(go_tab, unresolved), use.names = TRUE, fill = TRUE)
}
# Count every eligible pathway in the primary multiplicity family. An unresolved
# p-value remains NA, but cannot reduce the BH denominator for estimable tests.
go_tab[, backend_p_adjust := p.adjust]
go_tab[, p.adjust := NA_real_]
bh_family_size <- nrow(go_tab)
estimable_pathways <- is.finite(go_tab$pvalue)
go_tab[estimable_pathways, p.adjust := stats::p.adjust(pvalue, method = "BH", n = bh_family_size)]
fwrite(go_tab, file.path(result_dir, "GSEA_GO_BP_all.tsv"), sep = "\t")
fwrite(go_tab[!is.finite(NES) | !is.finite(pvalue) | !is.finite(p.adjust)],
       file.path(result_dir, "GSEA_GO_BP_unresolved.tsv"), sep = "\t")
writeLines(if (length(gsea_warnings)) gsea_warnings else "No fgsea warnings reported",
           file.path(result_dir, "GSEA_GO_BP_warnings.txt"))
fwrite(data.table(
  metric = c("all_pathway_rows", "rows_with_finite_NES_and_adjusted_p", "unresolved_pathway_rows",
             "BH_FDR05_pathways", "warnings_count", "BH_family_size", "BH_policy", "qvalue_role"),
  value = c(nrow(go_tab), sum(is.finite(go_tab$NES) & is.finite(go_tab$p.adjust)),
            sum(!is.finite(go_tab$NES) | !is.finite(go_tab$pvalue) | !is.finite(go_tab$p.adjust)),
            sum(go_tab$p.adjust < 0.05, na.rm = TRUE), length(gsea_warnings), bh_family_size,
            "BH denominator includes every eligible pathway; unresolved p-values remain NA",
            "backend qvalue retained for transparency; primary multiplicity result is full-family BH p.adjust")
), file.path(result_dir, "GSEA_GO_BP_test_audit.tsv"), sep = "\t")
if (nrow(go_tab)) {
  fwrite(go_tab[p.adjust < 0.05], file.path(result_dir, "GSEA_GO_BP_FDR05.tsv"), sep = "\t")
}

message("[6/8] Predeclared immune and interferon ssGSEA signatures")
signatures <- list(
  # IFI44L is intentionally excluded to avoid self-correlation when this score is
  # tested against IFI44L expression.
  IFN_ALPHA_RESPONSE = c("IFI6", "IFI27", "IFI35", "IFI44", "IFIT1", "IFIT2", "IFIT3",
                         "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "STAT1", "IRF7", "IRF9", "USP18"),
  IFN_GAMMA_RESPONSE = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP2",
                         "TAP1", "TAP2", "PSMB8", "PSMB9", "HLA-DRA", "HLA-DRB1", "B2M"),
  T_CELL_INFLAMED = c("CCL5", "CD27", "CD274", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
                      "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"),
  CYTOLYTIC_ACTIVITY = c("GZMA", "PRF1"),
  CD8_T_CELL = c("CD3D", "CD3E", "CD8A", "CD8B", "TRAC", "LCK", "NKG7", "GZMK"),
  TREG = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "CCR8", "TIGIT", "TNFRSF18"),
  NK_CELL = c("NKG7", "GNLY", "KLRD1", "NCR1", "PRF1", "XCL1", "XCL2", "FCGR3A"),
  T_CELL_EXHAUSTION = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CD244"),
  M2_MACROPHAGE = c("CD163", "MRC1", "MSR1", "C1QA", "C1QB", "C1QC", "IL10", "TGFB1"),
  DENDRITIC_CELL = c("FCER1A", "CD1C", "CLEC10A", "HLA-DRA", "ITGAX", "CLEC9A"),
  ANTIGEN_PRESENTATION = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "HLA-DRA"),
  CHECKPOINT_AXIS = c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"),
  EPITHELIAL_CONTENT = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1")
)
signature_members <- rbindlist(lapply(names(signatures), function(nm) {
  data.table(signature = nm, gene = signatures[[nm]], present = signatures[[nm]] %in% rownames(expr))
}))
fwrite(signature_members, file.path(result_dir, "immune_signature_definitions.tsv"), sep = "\t")

signatures_present <- lapply(signatures, function(x) intersect(x, rownames(expr)))
signatures_present <- signatures_present[lengths(signatures_present) >= 2]
ssgsea_scores <- gsva(ssgseaParam(expr, signatures_present, minSize = 2, maxSize = 500,
                                  normalize = TRUE), verbose = FALSE)
fwrite(data.table(signature = rownames(ssgsea_scores), as.data.table(ssgsea_scores)),
       file.path(result_dir, "immune_signature_ssGSEA_scores_wide.tsv"), sep = "\t")

sig_cor <- rbindlist(lapply(rownames(ssgsea_scores), function(sig) {
  ct <- suppressWarnings(cor.test(ifi, ssgsea_scores[sig, ], method = "spearman", exact = FALSE))
  data.table(signature = sig, spearman_rho = unname(ct$estimate), p_value = ct$p.value,
             n = length(ifi))
}))
sig_cor[, fdr_bh := p.adjust(p_value, method = "BH")]
sig_cor[, abs_spearman_rho := abs(spearman_rho)]
setorder(sig_cor, fdr_bh, -abs_spearman_rho)
sig_cor[, abs_spearman_rho := NULL]
fwrite(sig_cor, file.path(result_dir, "IFI44L_immune_signature_correlations.tsv"), sep = "\t")

message("[7/8] Checkpoint correlations with rank-based, target-excluded adjustment")
checkpoints <- intersect(c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TNFRSF4", "TNFRSF9"),
                         rownames(expr))
cov_names <- intersect(c("IFN_ALPHA_RESPONSE", "IFN_GAMMA_RESPONSE", "T_CELL_INFLAMED", "EPITHELIAL_CONTENT"),
                       rownames(ssgsea_scores))
cov_df <- as.data.frame(t(ssgsea_scores[cov_names, , drop = FALSE]))

# Preserve the historical raw-expression residual/Spearman calculation as a
# sensitivity analysis only. It is not a conventional partial Spearman test.
ifi_resid <- resid(lm(ifi ~ ., data = cov_df))
legacy_checkpoint_cor <- rbindlist(lapply(checkpoints, function(gene) {
  raw <- suppressWarnings(cor.test(ifi, expr[gene, ], method = "spearman", exact = FALSE))
  gene_resid <- resid(lm(expr[gene, ] ~ ., data = cov_df))
  adj <- suppressWarnings(cor.test(ifi_resid, gene_resid, method = "spearman", exact = FALSE))
  data.table(gene = gene,
             raw_spearman_rho = unname(raw$estimate), raw_p_value = raw$p.value,
             adjusted_spearman_rho = unname(adj$estimate), adjusted_p_value = adj$p.value,
             adjusted_for = paste(cov_names, collapse = ";"))
}))
legacy_checkpoint_cor[, raw_fdr_bh := p.adjust(raw_p_value, "BH")]
legacy_checkpoint_cor[, adjusted_fdr_bh := p.adjust(adjusted_p_value, "BH")]
legacy_checkpoint_cor[, method := "legacy Spearman correlation of raw-scale OLS residuals; uncorrected cor.test degrees of freedom"]
fwrite(legacy_checkpoint_cor,
       file.path(result_dir, "IFI44L_checkpoint_legacy_residual_spearman_sensitivity.tsv"), sep = "\t")

# Partial Spearman is Pearson correlation of rank residuals. The t reference
# uses n-k-2 degrees of freedom, where k is the rank of the covariate design
# excluding its intercept. This is an approximate, two-sided partial-rank test.
partial_spearman <- function(x, y, covariates) {
  values <- cbind(x = x, y = y, as.matrix(covariates))
  values <- values[apply(values, 1, function(z) all(is.finite(z))), , drop = FALSE]
  ranked <- apply(values, 2, rank, ties.method = "average")
  design <- cbind(Intercept = 1, ranked[, -(1:2), drop = FALSE])
  qr_design <- qr(design)
  x_residual <- qr.resid(qr_design, ranked[, "x"])
  y_residual <- qr.resid(qr_design, ranked[, "y"])
  rho <- unname(cor(x_residual, y_residual, method = "pearson"))
  df <- nrow(values) - qr_design$rank - 1L
  if (!is.finite(rho) || df <= 0) stop("Invalid partial Spearman model")
  statistic <- rho * sqrt(df / max(1 - rho^2, .Machine$double.eps))
  list(rho = rho, p_value = 2 * pt(abs(statistic), df = df, lower.tail = FALSE),
       n = nrow(values), covariate_rank = qr_design$rank - 1L, df = df)
}

adjustment_definitions <- list()
model_definitions <- list()
shared_rank_sensitivity <- list()
checkpoint_rows <- list()
for (gene in checkpoints) {
  target_signatures <- lapply(signatures[cov_names], function(members) {
    intersect(setdiff(members, c("IFI44L", gene)), rownames(expr))
  })
  if (any(lengths(target_signatures) < 2)) stop("Insufficient genes in target-excluded adjustment signature")
  target_scores <- gsva(ssgseaParam(expr, target_signatures, minSize = 2,
                                   maxSize = 500, normalize = TRUE), verbose = FALSE)
  target_covariates <- as.data.frame(t(target_scores[cov_names, , drop = FALSE]))
  adjusted <- partial_spearman(ifi, expr[gene, ], target_covariates)
  shared_adjusted <- partial_spearman(ifi, expr[gene, ], cov_df)
  raw <- suppressWarnings(cor.test(ifi, expr[gene, ], method = "spearman", exact = FALSE))
  checkpoint_rows[[gene]] <- data.table(
    gene = gene, raw_spearman_rho = unname(raw$estimate), raw_p_value = raw$p.value,
    adjusted_spearman_rho = adjusted$rho, adjusted_p_value = adjusted$p_value,
    adjusted_for = paste(cov_names, collapse = ";"), n = adjusted$n,
    covariate_rank = adjusted$covariate_rank, degrees_of_freedom = adjusted$df,
    method = "rank-all-variables; OLS residualization; Pearson correlation of rank residuals; two-sided approximate t test",
    adjustment_signature_policy = "exclude IFI44L and current checkpoint gene from every adjustment signature"
  )
  shared_rank_sensitivity[[gene]] <- data.table(
    gene = gene, adjusted_spearman_rho = shared_adjusted$rho,
    adjusted_p_value = shared_adjusted$p_value, n = shared_adjusted$n,
    covariate_rank = shared_adjusted$covariate_rank, degrees_of_freedom = shared_adjusted$df,
    method = "rank-based partial Spearman with shared modules before checkpoint-specific exclusion"
  )
  adjustment_definitions[[gene]] <- rbindlist(lapply(cov_names, function(signature) {
    members <- signatures[[signature]]
    data.table(checkpoint = gene, signature = signature, gene = members,
               present = members %in% rownames(expr),
               used_for_adjustment = members %in% target_signatures[[signature]])
  }))
  model_definitions[[gene]] <- data.table(
    checkpoint = gene, predictor = "IFI44L", n = adjusted$n,
    covariates = paste(cov_names, collapse = ";"), covariate_rank = adjusted$covariate_rank,
    degrees_of_freedom = adjusted$df,
    excluded_signature_members = paste(intersect(gene, unique(unlist(signatures[cov_names]))), collapse = ";"),
    signature_scoring = "ssGSEA on full expression matrix; signature members exclude IFI44L and current checkpoint",
    test = "Pearson correlation of rank residuals; t reference df=n-k-2",
    multiplicity = "BH across the nine checkpoint-specific primary tests"
  )
}
checkpoint_cor <- rbindlist(checkpoint_rows)
checkpoint_cor[, raw_fdr_bh := p.adjust(raw_p_value, "BH")]
checkpoint_cor[, adjusted_fdr_bh := p.adjust(adjusted_p_value, "BH")]
fwrite(checkpoint_cor, file.path(result_dir, "IFI44L_checkpoint_partial_correlations.tsv"), sep = "\t")
shared_rank_cor <- rbindlist(shared_rank_sensitivity)
shared_rank_cor[, adjusted_fdr_bh := p.adjust(adjusted_p_value, "BH")]
fwrite(shared_rank_cor, file.path(result_dir, "IFI44L_checkpoint_shared_module_rank_partial_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(adjustment_definitions),
       file.path(result_dir, "IFI44L_checkpoint_adjustment_signature_definitions.tsv"), sep = "\t")
fwrite(rbindlist(model_definitions),
       file.path(result_dir, "IFI44L_checkpoint_model_definitions.tsv"), sep = "\t")

message("[8/8] Publication-oriented diagnostic figures and session information")
volcano <- copy(de)
volcano[, status := fifelse(fdr_bh < 0.05 & expression_difference_high_minus_low > 0, "Higher",
                            fifelse(fdr_bh < 0.05 & expression_difference_high_minus_low < 0, "Lower", "Not significant"))]
volcano[, neglog10_fdr := -log10(pmax(fdr_bh, .Machine$double.xmin))]
p1 <- ggplot(volcano, aes(expression_difference_high_minus_low, neglog10_fdr, color = status)) +
  geom_point(size = 0.7, alpha = 0.55) +
  scale_color_manual(values = c(Higher = "#B2182B", Lower = "#2166AC", `Not significant` = "grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35) +
  labs(x = "Expression difference: high minus low (log2 RSEM units)", y = expression(-log[10](BH-FDR)), color = NULL) +
  theme_classic(base_size = 11)
ggsave(file.path(figure_dir, "bulk_median_group_volcano.pdf"), p1, width = 6.2, height = 5.2)
ggsave(file.path(figure_dir, "bulk_median_group_volcano.png"), p1, width = 6.2, height = 5.2, dpi = 400)

p2 <- ggplot(sig_cor, aes(reorder(signature, spearman_rho), spearman_rho, fill = fdr_bh < 0.05)) +
  geom_col(width = 0.75) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70"), guide = "none") +
  labs(x = NULL, y = "Spearman correlation with IFI44L") +
  theme_classic(base_size = 11)
ggsave(file.path(figure_dir, "IFI44L_immune_signature_correlations.pdf"), p2, width = 6.8, height = 5.8)
ggsave(file.path(figure_dir, "IFI44L_immune_signature_correlations.png"), p2, width = 6.8, height = 5.8, dpi = 400)

if (nrow(go_tab)) {
  go_sig <- go_tab[p.adjust < 0.05 & !is.na(NES)]
  go_plot <- rbind(
    go_sig[NES > 0][order(p.adjust)][1:min(.N, 10)],
    go_sig[NES < 0][order(p.adjust)][1:min(.N, 10)]
  )
  if (nrow(go_plot)) {
    go_plot[, direction := ifelse(NES >= 0, "Positive", "Negative")]
    go_plot[, Description_wrapped := vapply(
      Description,
      function(x) paste(strwrap(x, width = 45), collapse = "\n"),
      character(1)
    )]
    p3 <- ggplot(go_plot, aes(NES, reorder(Description_wrapped, NES), size = -log10(p.adjust), color = direction)) +
      geom_point(alpha = 0.85) +
      scale_color_manual(values = c(Positive = "#B2182B", Negative = "#2166AC")) +
      labs(x = "Normalized enrichment score", y = NULL, size = expression(-log[10](BH-FDR)), color = NULL) +
      theme_classic(base_size = 10)
    ggsave(file.path(figure_dir, "GSEA_GO_BP_top20.pdf"), p3, width = 10.5, height = 8.2)
    ggsave(file.path(figure_dir, "GSEA_GO_BP_top20.png"), p3, width = 10.5, height = 8.2, dpi = 400)
  }
}

writeLines(capture.output(sessionInfo()), file.path(log_dir, "bulk_R_sessionInfo.txt"))
message("Bulk TNBC reanalysis completed successfully.")
