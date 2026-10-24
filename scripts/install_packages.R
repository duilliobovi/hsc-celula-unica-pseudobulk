# Instala os pacotes usados no projeto (CRAN e Bioconductor).
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("data.table", "Matrix", "glmnet", "pROC", "ggplot2", "ggrepel", "scales",
          "cluster", "nnls", "jsonlite", "testthat")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)
bioc <- c("limma", "SummarizedExperiment")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
message("Pacotes prontos.")
