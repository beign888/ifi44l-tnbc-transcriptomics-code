#!/usr/bin/env Rscript

cran_packages <- c("data.table", "ggplot2", "survival", "survminer", "readxl")
bioconductor_packages <- c(
  "limma", "GSVA", "clusterProfiler", "org.Hs.eg.db", "GEOquery",
  "Biobase", "AnnotationDbi", "hgu133plus2.db"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

message("R package installation check complete. Compare resolved versions with environment/R-packages.tsv.")
