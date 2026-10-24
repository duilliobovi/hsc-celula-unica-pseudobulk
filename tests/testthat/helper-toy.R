# Compartimento de brinquedo: 3 doadores, HSC e dois grupos de progenitores.
toy_compartment <- function(n_per_group = 20L, n_genes = 30L, seed = 1L) {
  set.seed(seed)
  groups <- c(CONFIG$hsc_cluster, CONFIG$progenitor_clusters)
  cells <- expand.grid(i = seq_len(n_per_group), cluster = groups, donor = CONFIG$donors,
                       stringsAsFactors = FALSE)
  cells$BioClassification <- cells$cluster
  cells$is_hsc <- cells$cluster == CONFIG$hsc_cluster
  rate <- ifelse(cells$is_hsc, 6, 2)
  counts <- sapply(rate, function(r) stats::rpois(n_genes, r))
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  list(counts = Matrix::Matrix(counts, sparse = TRUE), cells = cells)
}
