test_that("target() appends a margin spec", {
  m <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.48, 0.52)))
  expect_length(m@targets, 1)
  expect_equal(m@targets[[1]]$tol, 0)
  expect_true(is.data.frame(m@targets[[1]]$df))
})

test_that("target() requires exactly one of p or n", {
  expect_error(target(raking(), data.frame(sex = "M")), "exactly one")
  expect_error(target(raking(), data.frame(sex = "M", p = 1, n = 1)),
               "exactly one")
})

test_that("calibrate() captures a formula and tolerance", {
  m <- raking() |> calibrate(mean(income) ~ 58000)
  expect_length(m@calibrations, 1)
  expect_true(inherits(m@calibrations[[1]]$formula, "formula"))
  expect_equal(m@calibrations[[1]]$tol, 0)

  m2 <- raking() |> calibrate(mean(income) ~ 58000, tol = 500)
  expect_equal(m2@calibrations[[1]]$tol, 500)
  expect_error(calibrate(raking(), mean(x) ~ 1, tol = -1), "tol")
})

test_that("penalty() sets the penalty and lambda", {
  m <- raking() |> penalty(entropy(), lambda = 2)
  expect_true(S7::S7_inherits(m@penalty, EntropyRegularizer))
  expect_equal(m@lambda, 2)
})

test_that("penalty() messages when replacing", {
  m <- raking() |> penalty(entropy())
  expect_message(penalty(m, sum_squares()), "replacing")
})

test_that("bounds() sets min and max", {
  m <- raking() |> bounds(min = 0.5, max = 2)
  expect_equal(m@bounds, list(min = 0.5, max = 2))
})

test_that("bounds() rejects min >= max", {
  expect_error(bounds(raking(), min = 2, max = 1), "must be <")
})

test_that("a full model chains together", {
  m <- raking(loss = "KL") |>
    target(data.frame(sex = c("M", "F"), p = c(0.48, 0.52)), tol = 0.01) |>
    calibrate(mean(income) ~ 58000) |>
    penalty(entropy()) |>
    bounds(min = 0.5, max = 2)
  expect_length(m@targets, 1)
  expect_length(m@calibrations, 1)
  expect_false(is.null(m@penalty))
  expect_equal(m@bounds$max, 2)
})
