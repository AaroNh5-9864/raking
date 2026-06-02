test_that("EqualityLoss stores a numeric target", {
  loss <- EqualityLoss(target = c(0.4, 0.6))
  expect_identical(loss@target, c(0.4, 0.6))
})

test_that("EqualityLoss prox returns the target unchanged", {
  loss <- EqualityLoss(target = c(0.4, 0.6))
  expect_identical(prox(loss, v = c(0.1, 0.9), lam = 0.02), c(0.4, 0.6))
})

test_that("LeastSquaresLoss defaults diag_weight to 1", {
  loss <- LeastSquaresLoss(target = c(1, 2))
  expect_equal(loss@diag_weight, 1)
})

test_that("LeastSquaresLoss prox matches the closed-form solution", {
  loss <- LeastSquaresLoss(target = c(2, 4))
  expect_equal(prox(loss, v = c(0, 0), lam = 1), c(1, 2))
})

test_that("LeastSquaresLoss evaluate sums weighted squared deviations", {
  loss <- LeastSquaresLoss(target = c(2, 4))
  expect_equal(evaluate(loss, v = c(3, 4)), 1)
})

test_that("InequalityLoss rejects lower > upper", {
  expect_error(
    InequalityLoss(target = c(0.5), lower = c(0.2), upper = c(0.1)),
    "lower"
  )
})

test_that("InequalityLoss prox clips into the band (matches Python rsw)", {
  loss <- InequalityLoss(
    target = c(0.5, 0.5),
    lower  = c(-0.1, -0.1),
    upper  = c(0.1, 0.1)
  )
  expect_equal(prox(loss, v = c(0.8, 0.45), lam = 1), c(0.6, 0.45))
})

test_that("KLLoss prox matches Python rsw (Lambert-W)", {
  loss <- KLLoss(target = c(0.4, 0.6), scale = 0.5)
  expect_equal(
    prox(loss, v = c(0.5, 0.5), lam = 1),
    c(0.245033929401, 0.317782008182)
  )
})

test_that("KLLoss evaluate is zero at the target and positive elsewhere", {
  loss <- KLLoss(target = c(0.4, 0.6), scale = 0.5)
  expect_equal(evaluate(loss, c(0.4, 0.6)), 0)
  expect_equal(evaluate(loss, c(0.5, 0.5)), 0.010205498630063758)
})

test_that("L1Loss prox is soft-thresholding around the target", {
  loss <- L1Loss(target = 0.5)
  expect_equal(prox(loss, v = 0.6, lam = 0.02), 0.58)  # 0.5 + soft(0.1, 0.02)
  expect_equal(prox(loss, v = 0.51, lam = 0.02), 0.5)  # within threshold -> target
})

test_that("L1Loss evaluate is the sum of absolute deviations", {
  loss <- L1Loss(target = c(0.4, 0.6))
  expect_equal(evaluate(loss, c(0.5, 0.5)), 0.2)
})
