# Misturas pseudo-bulk montadas com células reais e fração de HSC conhecida.
#
# Diferente de simular perfis, aqui cada mistura é a soma das contagens de
# células de verdade de um doador. A variação entre células, a esparsidade
# do 10x e as diferenças de profundidade entram de graça na mistura.

#' Sorteia índices sem recorrer a sample(x) com x de tamanho 1 (armadilha do R).
pick <- function(idx, size) {
  if (size == 0) return(integer(0))
  idx[sample.int(length(idx), size, replace = size > length(idx))]
}

#' Gera n misturas de `cells_per` células a partir dos doadores indicados.
#'
#' Fração de HSC uniforme na faixa pedida; os progenitores restantes são
#' divididos entre os grupos por pesos Dirichlet, para variar a composição.
#' Devolve as contagens (genes x misturas), a fração real e o doador.
simulate_pseudobulks <- function(compartment, donors, n, seed,
                                 cells_per = CONFIG$cells_per_mixture,
                                 fraction_range = CONFIG$hsc_fraction_range) {
  set.seed(seed)
  cells <- compartment$cells
  clusters <- CONFIG$progenitor_clusters
  rows <- vector("list", n)
  fraction <- numeric(n)
  donor_of <- character(n)
  for (i in seq_len(n)) {
    d <- donors[sample.int(length(donors), 1)]
    n_hsc <- round(stats::runif(1, fraction_range[1], fraction_range[2]) * cells_per)
    weights <- stats::rgamma(length(clusters), shape = 1)
    n_k <- as.vector(stats::rmultinom(1, cells_per - n_hsc, weights / sum(weights)))
    chosen <- pick(which(cells$donor == d & cells$is_hsc), n_hsc)
    for (k in seq_along(clusters)) {
      pool <- which(cells$donor == d & cells$BioClassification == clusters[k])
      chosen <- c(chosen, pick(pool, n_k[k]))
    }
    rows[[i]] <- chosen
    fraction[i] <- n_hsc / cells_per
    donor_of[i] <- d
  }
  # Matriz células x misturas com o número de vezes que cada célula entrou.
  indicator <- Matrix::sparseMatrix(i = unlist(rows), j = rep(seq_len(n), lengths(rows)),
                                    x = 1, dims = c(ncol(compartment$counts), n))
  counts <- as.matrix(compartment$counts %*% indicator)
  list(counts = counts, fraction = fraction, donor = donor_of)
}
