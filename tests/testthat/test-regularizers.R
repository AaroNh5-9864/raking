test_that("ZeroRegularizer prox returns weights unchanged", {
  reg <- ZeroRegularizer()
  expect_identical(prox(reg, c(-0.5, 0, 0.5), lam = 0.02), c(-0.5, 0, 0.5))
})

test_that("EntropyRegularizer prox matches Python rsw (no limit)", {
  reg <- EntropyRegularizer()
  expect_equal(
    prox(reg, c(-0.5, 0, 0.5), lam = 0.02),
    c(0.000000000005, 0.042952983411, 0.494100333561)
  )
})

test_that("EntropyRegularizer limit clips toward uniform", {
  reg <- EntropyRegularizer(limit = 3)
  expect_equal(
    prox(reg, c(-0.5, 0, 0.5), lam = 0.02),
    c(0.111111111111, 0.111111111111, 0.494100333561)
  )
})

test_that("EntropyRegularizer rejects limit <= 1", {
  expect_error(EntropyRegularizer(limit = 0.5), "greater than 1")
})

test_that("KLRegularizer prox matches Python rsw", {
  reg <- KLRegularizer(prior = c(0.2, 0.3, 0.5))
  expect_equal(
    prox(reg, c(-0.5, 0, 0.5), lam = 0.02),
    c(0.000000000001, 0.027669230107, 0.480783807865)
  )
})

test_that("SumSquaresRegularizer shrinks toward zero", {
  reg <- SumSquaresRegularizer()
  expect_equal(prox(reg, c(0.2, 0.4), lam = 0.5), c(0.1, 0.2))
})
