test_that("rake_control returns the documented defaults", {
  ctrl <- rake_control()
  expect_equal(ctrl$margin_tol, 1e-4)
  expect_equal(ctrl$rho, 50)
  expect_equal(ctrl$maxiter, 5000)
})

test_that("rake_control nulls margin_tol when eps is given explicitly", {
  ctrl <- rake_control(eps_abs = 1e-7)
  expect_null(ctrl$margin_tol)
  expect_equal(ctrl$eps_abs, 1e-7)
})

test_that("resolve_control scales margin_tol by problem size", {
  ctrl <- resolve_control(rake_control(margin_tol = 1e-3), m = 4, n = 10)
  expect_equal(ctrl$eps_abs, 1e-3 / sqrt(4 + 2 * 10))
  expect_equal(ctrl$eps_abs, ctrl$eps_rel)
})

test_that("resolve_control honors raw eps from a plain list", {
  ctrl <- resolve_control(list(eps_abs = 1e-8, eps_rel = 1e-8), m = 4, n = 10)
  expect_equal(ctrl$eps_abs, 1e-8)
})

test_that("fit() accepts rake_control() and settings passed via ...", {
  data <- data.frame(sex = c("M", "F", "M", "F"))
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  res1 <- fit(model, data, control = rake_control(margin_tol = 1e-3))
  res2 <- fit(model, data, eps_abs = 1e-8, eps_rel = 1e-8)  # via ...
  expect_true(S7::S7_inherits(res1, raking_fit))
  expect_true(S7::S7_inherits(res2, raking_fit))
})
