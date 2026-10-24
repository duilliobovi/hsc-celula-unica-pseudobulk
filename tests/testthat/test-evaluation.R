test_that("métricas em casos conhecidos", {
  y <- c(0.1, 0.2, 0.3, 0.4)
  expect_equal(ccc(y, y), 1)
  expect_equal(c_index(y, y), 1)
  expect_equal(c_index(y, -y), 0)
  expect_equal(rmse(y, y + 0.1), 0.1)
  expect_lt(ccc(y, 2 * y), 0.9)
})

test_that("NNLS com perfis por célula recupera a fração de células", {
  set.seed(3)
  groups <- c(CONFIG$hsc_cluster, CONFIG$progenitor_clusters)
  S <- matrix(stats::rgamma(40 * length(groups), 2, 0.1), nrow = 40,
              dimnames = list(paste0("g", 1:40), groups))
  n_cells <- c(30, rep(10, length(groups) - 1))
  model <- fit_nnls(NULL, list(profiles = S, genes = rownames(S)))
  expect_equal(model$predict(S %*% n_cells), 30 / sum(n_cells), tolerance = 1e-6)
})

test_that("normalização log por célula preserva zeros", {
  m <- Matrix::Matrix(matrix(c(0, 2, 3, 0), 2), sparse = TRUE)
  out <- log_normalize(m)
  expect_equal(out[1, 1], 0)
  expect_equal(out[2, 1], log1p(1e4))
})
