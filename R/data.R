# Download, leitura e agregação dos dados de célula única.

#' Baixa um arquivo para data/, se ainda não estiver no disco.
download_if_missing <- function(url, dest_dir = CONFIG$data_dir) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(dest_dir, basename(url))
  if (!file.exists(dest)) {
    old <- options(timeout = max(1800, getOption("timeout")))
    on.exit(options(old))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

#' Lê o objeto dos autores e devolve contagens esparsas e tabela de células.
#'
#' O arquivo é um RangedSummarizedExperiment; só a matriz de contagens e a
#' anotação das células são usadas. O doador sai do nome da amostra
#' ("BMMC_D1T1" -> "D1"); T1 e T2 são duas corridas técnicas do mesmo doador.
load_granja <- function(path) {
  se <- readRDS(path)
  cells <- as.data.frame(SummarizedExperiment::colData(se))
  cells$donor <- sub("T[0-9]+$", "", sub("^[A-Z0-9]+_", "", cells$Group))
  cells$tissue <- sub("_.*$", "", cells$Group)
  counts <- SummarizedExperiment::assay(se, "counts")
  list(counts = counts, cells = cells[, c("Group", "donor", "tissue", "nUMI",
                                          "nGene", "BioClassification", "UMAP1",
                                          "UMAP2")])
}

#' Restringe ao compartimento CD34+ (HSC + progenitores) dos doadores com HSC.
cd34_compartment <- function(data) {
  keep <- data$cells$BioClassification %in% c(CONFIG$hsc_cluster,
                                              CONFIG$progenitor_clusters) &
    data$cells$donor %in% CONFIG$donors
  cells <- data$cells[keep, ]
  cells$is_hsc <- cells$BioClassification == CONFIG$hsc_cluster
  list(counts = data$counts[, keep], cells = cells)
}

#' Soma as contagens das células de cada grupo (pseudo-bulk).
#'
#' Somar contagens de células do mesmo doador e tipo produz uma "amostra"
#' com a mesma natureza de um RNA-seq bulk. É a unidade correta para testes:
#' a replicação biológica está nos doadores, não nas células.
pseudobulk <- function(counts, groups) {
  groups <- factor(groups)
  design <- Matrix::sparse.model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  as.matrix(counts %*% design)
}

#' Expressão normalizada por célula: log1p(contagens / nUMI x 10.000).
log_normalize <- function(counts) {
  size <- Matrix::colSums(counts)
  scaled <- counts %*% Matrix::Diagonal(x = 1e4 / size)
  scaled@x <- log1p(scaled@x)
  colnames(scaled) <- colnames(counts)
  scaled
}

#' log2(CPM + 1) de uma matriz densa de contagens (genes x amostras).
log_cpm <- function(counts) {
  log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
}
