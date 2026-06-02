test_that("entropy() builds an EntropyRegularizer", {
  expect_true(S7::S7_inherits(entropy(), EntropyRegularizer))
})

test_that("kl() builds a KLRegularizer carrying the prior", {
  p <- kl(prior = c(0.2, 0.3, 0.5))
  expect_true(S7::S7_inherits(p, KLRegularizer))
  expect_equal(p@prior, c(0.2, 0.3, 0.5))
})

test_that("sum_squares() builds a SumSquaresRegularizer", {
  expect_true(S7::S7_inherits(sum_squares(), SumSquaresRegularizer))
})
