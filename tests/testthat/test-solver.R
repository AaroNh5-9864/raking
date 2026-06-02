test_that("projection_simplex matches Python rsw", {
  expect_equal(
    projection_simplex(c(0.3, -0.2, 1.1, 0.4)),
    c(0.033333333333, 0, 0.833333333333, 0.133333333333)
  )
  expect_equal(projection_simplex(c(5, 1, 2)), c(1, 0, 0))
})

test_that("projection_simplex returns a valid distribution", {
  p <- projection_simplex(c(-1, -2, -3))
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0))
})

test_that("projection_bounded_simplex matches scipy projection", {
  expect_equal(projection_bounded_simplex(c(0.3, -0.2, 1.1, 0.4), 0.1, 0.5),
               c(0.15, 0.1, 0.5, 0.25), tolerance = 1e-7)
  expect_equal(projection_bounded_simplex(c(5, 1, 2), 0.2, 0.5),
               c(0.5, 0.2, 0.3), tolerance = 1e-7)
})

test_that("projection_bounded_simplex respects bounds and sums to z", {
  p <- projection_bounded_simplex(c(0.3, -0.2, 1.1, 0.4), 0.1, 0.5)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0.1 - 1e-9 & p <= 0.5 + 1e-9))
})

test_that("admm solves equality + entropy (matches Python rsw)", {
  F <- matrix(c(1, 0, 1, 0, 1, 0,
                1, 1, 0, 0, 1, 1), nrow = 2, byrow = TRUE)
  sol <- admm(F, list(EqualityLoss(target = c(0.5, 0.5))),
              EntropyRegularizer(), lam = 1, rho = 50,
              eps_abs = 1e-8, eps_rel = 1e-8)
  expect_equal(sol$w_best,
               c(0.125, 0.125, 0.25, 0.25, 0.125, 0.125),
               tolerance = 1e-5)
  expect_equal(sum(sol$w_best), 1, tolerance = 1e-6)
})

test_that("admm: least-squares + inequality + entropy (matches Python rsw)", {
  F <- matrix(c(1, 0, 1, 0, 1, 0,
                0, 1, 0, 1, 0, 1,
                1, 1, 1, 0, 0, 0,
                0, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
  losses <- list(
    LeastSquaresLoss(target = c(0.5, 0.5)),
    InequalityLoss(target = c(0.5, 0.5),
                   lower = c(-0.1, -0.1), upper = c(0.1, 0.1))
  )
  sol <- admm(F, losses, EntropyRegularizer(), lam = 1, rho = 50,
              eps_abs = 1e-8, eps_rel = 1e-8)
  expect_equal(sol$w_best, rep(1 / 6, 6), tolerance = 1e-5)
})

test_that("evaluate() returns a loss's objective value", {
  loss <- LeastSquaresLoss(target = c(0.5, 0.5))
  # sum of squared deviations: (0.6-0.5)^2 + (0.4-0.5)^2 = 0.02
  expect_equal(evaluate(loss, c(0.6, 0.4)), 0.02)
})

