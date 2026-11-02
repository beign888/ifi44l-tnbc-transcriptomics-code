#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(GSVA)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
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
raw_root <- normalizePath(
  argument_value("--raw-root", dirname(project_dir)),
  mustWork = TRUE
)
input_dir <- normalizePath(
  file.path(raw_root, "数据", "差异富集分析", "3 差异分析"),
  mustWork = TRUE
)
result_dir <- file.path(project_dir, "results", "bulk")
figure_dir <- file.path(project_dir, "figures", "bulk")
log_dir <- file.path(project_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(input_dir, "symbol.BRCA.txt")
tnbc_file <- file.path(input_dir, "三阴.txt")

message("[1/8] Reading header and matching the TNBC sample list")
header <- names(fread(matrix_file, nrows = 0, check.names = FALSE))
all_samples <- header[-1]
tnbc_requested <- unique(trimws(readLines(tnbc_file, warn = FALSE)))
tnbc_requested <- tnbc_requested[nzchar(tnbc_requested)]
tnbc_samples <- tnbc_requested[tnbc_requested %in% all_samples]
tnbc_missing <- setdiff(tnbc_requested, all_samples)

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

message("[2/8] Loading only the 187 matched TNBC columns")
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
fit_cont <- eBayes(lmFit(expr, design_cont), trend = TRUE, robust = TRUE)
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
fit_group <- contrasts.fit(lmFit(expr, design_group), makeContrasts(High - Low, levels = design_group))
fit_group <- eBayes(fit_group, trend = TRUE, robust = TRUE)
de <- as.data.table(topTable(fit_group, number = Inf, adjust.method = "BH", sort.by = "none"),
                    keep.rownames = "gene")
setnames(de, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
         c("log2_difference_high_minus_low", "mean_expression", "moderated_t", "p_value", "fdr_bh", "lods"))
de[, abs_moderated_t := abs(moderated_t)]
setorder(de, fdr_bh, -abs_moderated_t)
de[, abs_moderated_t := NULL]
fwrite(de, file.path(result_dir, "median_group_differential_expression.tsv"), sep = "\t")
fwrite(de[fdr_bh < 0.05 & abs(log2_difference_high_minus_low) >= 0.585],
       file.path(result_dir, "median_group_DE_FDR05_absLFC0585.tsv"), sep = "\t")

message("[5/8] GO Biological Process GSEA using the full continuous-model ranking")
rank_dt <- cont[gene != "IFI44L" & is.finite(moderated_t), .(gene, moderated_t)]
mapping <- suppressMessages(bitr(rank_dt$gene, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db))
rank_mapped <- merge(rank_dt, as.data.table(mapping), by.x = "gene", by.y = "SYMBOL")
rank_mapped <- rank_mapped[order(-abs(moderated_t))][!duplicated(ENTREZID)]
gene_list <- rank_mapped$moderated_t
names(gene_list) <- rank_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go <- suppressMessages(gseGO(
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
))
go_tab <- as.data.table(as.data.frame(gsea_go))
fwrite(go_tab, file.path(result_dir, "GSEA_GO_BP_all.tsv"), sep = "\t")
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

message("[7/8] Checkpoint correlations before and after interferon/immune-state adjustment")
checkpoints <- intersect(c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TNFRSF4", "TNFRSF9"),
                         rownames(expr))
cov_names <- intersect(c("IFN_ALPHA_RESPONSE", "IFN_GAMMA_RESPONSE", "T_CELL_INFLAMED", "EPITHELIAL_CONTENT"),
                       rownames(ssgsea_scores))
cov_df <- as.data.frame(t(ssgsea_scores[cov_names, , drop = FALSE]))
ifi_resid <- resid(lm(ifi ~ ., data = cov_df))
checkpoint_cor <- rbindlist(lapply(checkpoints, function(gene) {
  raw <- suppressWarnings(cor.test(ifi, expr[gene, ], method = "spearman", exact = FALSE))
  gene_resid <- resid(lm(expr[gene, ] ~ ., data = cov_df))
  adj <- suppressWarnings(cor.test(ifi_resid, gene_resid, method = "spearman", exact = FALSE))
  data.table(gene = gene,
             raw_spearman_rho = unname(raw$estimate), raw_p_value = raw$p.value,
             adjusted_spearman_rho = unname(adj$estimate), adjusted_p_value = adj$p.value,
             adjusted_for = paste(cov_names, collapse = ";"))
}))
checkpoint_cor[, raw_fdr_bh := p.adjust(raw_p_value, "BH")]
checkpoint_cor[, adjusted_fdr_bh := p.adjust(adjusted_p_value, "BH")]
fwrite(checkpoint_cor, file.path(result_dir, "IFI44L_checkpoint_partial_correlations.tsv"), sep = "\t")

message("[8/8] Publication-oriented diagnostic figures and session information")
volcano <- copy(de)
volcano[, status := fifelse(fdr_bh < 0.05 & log2_difference_high_minus_low >= 0.585, "Higher",
                            fifelse(fdr_bh < 0.05 & log2_difference_high_minus_low <= -0.585, "Lower", "Not selected"))]
volcano[, neglog10_fdr := -log10(pmax(fdr_bh, .Machine$double.xmin))]
p1 <- ggplot(volcano, aes(log2_difference_high_minus_low, neglog10_fdr, color = status)) +
  geom_point(size = 0.7, alpha = 0.55) +
  scale_color_manual(values = c(Higher = "#B2182B", Lower = "#2166AC", `Not selected` = "grey75")) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = 2, linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35) +
  labs(x = "Expression difference: IFI44L-high minus low", y = expression(-log[10](BH-FDR)), color = NULL) +
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
