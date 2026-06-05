# End-to-end fit() tests. tol = 0 (the default) is exact matching via
# EqualityLoss, so the global `loss` is irrelevant there; tol > 0 is a penalized
# band (the loss pulls toward the target, clipped into the band), where `loss`
# matters. Reference values were generated from a faithful Python port and
# cross-checked against rsw.

test_that("fit() matches margins and the mean exactly (rsw-verified)", {
  data <- data.frame(
    sex = c("M", "F", "M", "F", "M"),
    income = c(40000, 60000, 50000, 70000, 45000),
    stringsAsFactors = FALSE
  )
  model <- raking(loss = "KL") |>
    target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
    calibrate(mean(income) ~ 52000) |>
    penalty(entropy())
  res <- fit(model, data, control = list(eps_abs = 1e-8, eps_rel = 1e-8))

  expect_true(S7::S7_inherits(res, raking_fit))
  expect_equal(res@weights,
    c(1.553704729, 2.114690844, 0.283086240, 0.385309158, 0.663209029),
    tolerance = 1e-4)
  expect_equal(mean(res@weights), 1, tolerance = 1e-6)
  # sex margin and mean income are hit exactly
  expect_equal(res@balance$achieved[res@balance$level == "M"], 0.5,
               tolerance = 1e-4)
  expect_equal(res@balance$achieved[res@balance$variable == "income"], 52000,
               tolerance = 1)
  expect_true(all(res@balance$satisfied))
})

test_that("fit() errors on an empty model", {
  expect_error(fit(raking(), data.frame(x = 1)), "no targets or calibrations")
})

test_that("a fully-determined 1-DOF margin is matched exactly", {
  # 3 M / 2 F, target 0.5/0.5 -> M total 2.5 (0.8333 each), F total 2.5 (1.25).
  # At tol = 0 the match is exact, so the global loss is irrelevant; all agree.
  data <- data.frame(sex = c("M", "M", "M", "F", "F"), stringsAsFactors = FALSE)
  for (ls in c("KL", "L2", "L1", "chi_square")) {
    res <- fit(
      raking(loss = ls) |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
      data, control = list(eps_abs = 1e-8, eps_rel = 1e-8)
    )
    expect_equal(res@weights, c(5/6, 5/6, 5/6, 1.25, 1.25), tolerance = 1e-4)
  }
})

test_that("fit() prints a summary", {
  data <- data.frame(sex = c("M", "F", "M", "F"))
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  res <- fit(model, data)
  expect_output(print(res), "raking_fit")
})

test_that("fit() with tol > 0 enforces a range target", {
  data <- data.frame(sex = c("M", "M", "M", "F", "F"), stringsAsFactors = FALSE)
  model <- raking(loss = "KL") |>
    target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)), tol = 0.05) |>
    penalty(entropy())
  res <- fit(model, data, control = list(eps_abs = 1e-8, eps_rel = 1e-8))

  expect_equal(res@weights,
    c(0.9166666686, 0.9166666686, 0.9166666686, 1.124999997, 1.124999997),
    tolerance = 1e-4)
  # achieved M share is pulled to the +0.05 edge of the [0.45, 0.55] band
  expect_lte(abs(res@balance$achieved[1] - 0.5), 0.05 + 1e-6)
})

test_that("base weights shift weights within an under-determined cell", {
  # 3 M / 1 F, target M = 0.75: the 3 M's share a total of 3.0 with freedom.
  data <- data.frame(sex = c("M", "M", "M", "F"), stringsAsFactors = FALSE)
  tgt   <- data.frame(sex = c("M", "F"), p = c(0.75, 0.25))
  ctrl  <- list(eps_abs = 1e-9, eps_rel = 1e-9, maxiter = 20000)

  flat <- fit(raking() |> target(tgt) |> penalty(entropy()), data, control = ctrl)
  expect_equal(flat@weights, c(1, 1, 1, 1), tolerance = 1e-3)

  # base upweights the third respondent; entropy becomes KL toward the base
  wt <- fit(raking() |> target(tgt) |> penalty(entropy()),
            data, base = c(1, 1, 3, 1), control = ctrl)
  expect_equal(wt@weights, c(0.6, 0.6, 1.8, 1.0), tolerance = 1e-3)
})

test_that("fit() accepts base as a column name", {
  data <- data.frame(sex = c("M", "M", "M", "F"), dw = c(1, 1, 3, 1),
                     stringsAsFactors = FALSE)
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.75, 0.25))) |>
               penalty(entropy()),
             data, base = "dw", control = list(eps_abs = 1e-9, eps_rel = 1e-9,
                                               maxiter = 20000))
  expect_equal(res@weights[3], 1.8, tolerance = 1e-3)
})

test_that("fit() rejects non-positive base weights", {
  data <- data.frame(sex = c("M", "F"))
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data, base = c(1, 0)), "strictly positive")
})

test_that("fit() enforces weight bounds (feasible range + binding upper)", {
  # A tol band gives the slack needed for bounds to be feasible; the upper
  # bound binds on the heavy F weight.
  data <- data.frame(sex = c("M", "M", "M", "M", "F"), stringsAsFactors = FALSE)
  model <- raking(loss = "KL") |>
    target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)), tol = 0.2) |>
    penalty(entropy()) |>
    bounds(min = 0.7, max = 1.5)
  res <- fit(model, data, control = list(eps_abs = 1e-8, eps_rel = 1e-8))

  expect_true(all(res@weights >= 0.7 - 1e-6 & res@weights <= 1.5 + 1e-6))
  expect_equal(max(res@weights), 1.5, tolerance = 1e-3)
})

test_that("fit() handles a joint (sex x age) margin, matched exactly", {
  data <- data.frame(
    sex = c("M", "M", "F", "F", "M", "F"),
    age = c("Y", "O", "Y", "O", "Y", "O"),
    stringsAsFactors = FALSE
  )
  joint <- data.frame(
    sex = c("M", "M", "F", "F"),
    age = c("Y", "O", "Y", "O"),
    n   = c(1, 1, 1, 1)              # equal cells
  )
  res <- fit(raking() |> target(joint), data,
             control = list(eps_abs = 1e-8, eps_rel = 1e-8))
  # Each cell = 0.25 (weight 1.5); two-member cells split evenly to 0.75.
  # Rows: M:Y M:O F:Y F:O M:Y F:O
  expect_equal(res@weights, c(0.75, 1.5, 1.5, 0.75, 0.75, 0.75),
               tolerance = 1e-4)
})

test_that("calibrate() target can be a variable from the caller's scope", {
  # Regression: the formula's RHS must be evaluated in the formula's own
  # environment, not inside fit(). A literal (mean(income) ~ 58000) always
  # works because it needs no environment; a variable only resolves if fit()
  # uses environment(formula). Wrapping in a function puts `mu` in a scope
  # that fit()'s internals cannot see, so this fails until that is fixed.
  data <- data.frame(income = c(40000, 60000, 50000, 70000, 45000))
  fit_with <- function(mu) {
    model <- raking() |> calibrate(mean(income) ~ mu) |> penalty(entropy())
    fit(model, data, control = list(eps_abs = 1e-8, eps_rel = 1e-8))
  }
  res <- fit_with(52000)
  expect_equal(res@balance$achieved[res@balance$variable == "income"], 52000,
               tolerance = 1)
})
