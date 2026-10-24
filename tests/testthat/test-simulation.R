test_that("pick funciona com um só índice (armadilha do sample)", {
  set.seed(1)
  expect_equal(pick(7L, 3), c(7L, 7L, 7L))
  expect_length(pick(1:5, 0), 0)
})

test_that("cada pseudo-bulk soma exatamente as células sorteadas", {
  cp <- toy_compartment()
  sim <- simulate_pseudobulks(cp, "D1", n = 5, seed = 2, cells_per = 50,
                              fraction_range = c(0.2, 0.2))
  expect_equal(sim$fraction, rep(0.2, 5))
  expect_true(all(sim$donor == "D1"))
  # a soma das contagens fica entre 50 células "frias" e 50 "quentes"
  totals <- colSums(sim$counts)
  expect_true(all(totals > 50 * 30 * 1 & totals < 50 * 30 * 8))
})

test_that("pseudo-bulk agrega por grupo sem perder contagens", {
  cp <- toy_compartment()
  pb <- pseudobulk(cp$counts, cp$cells$donor)
  expect_equal(sum(pb), sum(cp$counts))
  expect_equal(colnames(pb), CONFIG$donors)
})
