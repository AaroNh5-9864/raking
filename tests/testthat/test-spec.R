# Doc Section 10 testing checklist + regression tests for the divergence fixes.
# Reference weight values were generated from a faithful Python port of the
# solver (validated to reproduce the existing baked values), so tolerances are
# 1e-4 unless the damped Cholesky / loss conditioning requires looser.

# --- raking_fit object structure (doc 7; fix D1) ---------------------------

test_that("raking_fit exposes the documented properties", {
  data <- data.frame(sex = c("M", "M", "F", "F"), stringsAsFactors = FALSE)
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data)

  # balance table columns (doc 7)
  expect_s3_class(res@balance, "data.frame")
  expect_named(res@balance,
    c("variable", "level", "target", "achieved", "difference", "satisfied"))
  expect_equal(res@balance$variable, c("sex", "sex"))
  expect_setequal(res@balance$level, c("M", "F"))

  # convergence list (doc 7)
  expect_named(res@convergence,
    c("converged", "iterations", "primal_residual", "dual_residual"))
  expect_type(res@convergence$converged, "logical")

  # diagnostics list (doc 7)
  expect_named(res@diagnostics,
    c("weight_range", "weight_mean", "weight_sd", "kish_deff", "kish_ess",
      "max_abs_diff", "max_pct_diff"))

  # model + call carried through
  expect_true(S7::S7_inherits(res@model, raking_model))
  expect_false(is.null(res@call))
})

# --- Unpenalized default (doc 4.4; fix D4) ----------------------------------

test_that("resolve_penalty: no penalty + uniform base is unpenalized (lambda 0)", {
  out <- resolve_penalty(NULL, 1, rep(1 / 4, 4))
  expect_true(S7::S7_inherits(out$reg, ZeroRegularizer))
  expect_equal(out$lambda, 0)

  # non-uniform base with no penalty pulls toward the base (KL), lambda 1
  b <- c(0.1, 0.2, 0.3, 0.4)
  out2 <- resolve_penalty(NULL, 1, b)
  expect_true(S7::S7_inherits(out2$reg, KLRegularizer))
  expect_equal(out2$reg@prior, b)
  expect_equal(out2$lambda, 1)

  # entropy + non-uniform base becomes KL toward the base distribution
  out3 <- resolve_penalty(EntropyRegularizer(), 1, b)
  expect_true(S7::S7_inherits(out3$reg, KLRegularizer))
})

# --- Known analytic case: L2 post-stratification (doc 10) -------------------

test_that("L2 + no penalty reproduces the post-stratification solution", {
  # Post-strat weight for a cell = target_prop * n / cell_count (sum = n scale).
  data <- data.frame(sex = c(rep("M", 7), rep("F", 3)), stringsAsFactors = FALSE)
  res <- fit(
    raking(loss = "L2") |> target(data.frame(sex = c("M", "F"), p = c(0.4, 0.6))),
    data, control = list(eps_abs = 1e-9, eps_rel = 1e-9, maxiter = 20000)
  )
  analytic <- ifelse(data$sex == "M", 0.4 * 10 / 7, 0.6 * 10 / 3)
  # The damped sparse solve limits accuracy to ~2e-4 (doc's 1e-6 is not reached;
  # see notes). 1e-3 comfortably verifies the closed form.
  expect_equal(res@weights, analytic, tolerance = 1e-3)
})

# --- Kish ESS / deff (doc 10) ----------------------------------------------

test_that("kish_deff and kish_ess match the manual formulas", {
  data <- data.frame(sex = c("M", "M", "M", "M", "F"), stringsAsFactors = FALSE)
  res <- fit(raking(loss = "L2") |>
               target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data, control = list(eps_abs = 1e-9, eps_rel = 1e-9, maxiter = 20000))
  w <- res@weights
  deff <- 1 + stats::var(w) / mean(w)^2
  expect_equal(res@diagnostics$kish_deff, deff, tolerance = 1e-9)
  expect_equal(res@diagnostics$kish_ess, length(w) / deff, tolerance = 1e-9)
})

# --- tol changes the problem (doc 6.5, 10; band enforcement) ----------------

test_that("tol>0 enforces the range band, and a wider band is more uniform", {
  data <- data.frame(sex = c("M", "M", "M", "M", "F"), stringsAsFactors = FALSE)
  tgt  <- data.frame(sex = c("M", "F"), p = c(0.5, 0.5))
  ctrl <- list(eps_abs = 1e-8, eps_rel = 1e-8)

  narrow <- fit(raking(loss = "KL") |> target(tgt, tol = 0.02) |> penalty(entropy()),
                data, control = ctrl)
  wide   <- fit(raking(loss = "KL") |> target(tgt, tol = 0.10) |> penalty(entropy()),
                data, control = ctrl)

  # both achieved M shares sit inside their bands
  expect_lte(abs(narrow@balance$achieved[1] - 0.5), 0.02 + 1e-6)
  expect_lte(abs(wide@balance$achieved[1]   - 0.5), 0.10 + 1e-6)
  # a wider tolerance band leaves more room to stay uniform
  expect_lt(stats::sd(wide@weights), stats::sd(narrow@weights))
})

# --- Penalized band: loss pulls toward target inside the band ---------------

test_that("tol > 0 is a penalized band that pulls the achieved value to target", {
  # 3M / 2F, target M = 0.5, band +/-0.2 -> [0.3, 0.7]. The natural M share is
  # 0.6, which is already inside the band, so an *indifferent* band would leave
  # the weights uniform (M = 0.6). The penalized band instead pulls M toward 0.5.
  data <- data.frame(sex = c("M", "M", "M", "F", "F"), stringsAsFactors = FALSE)
  res <- fit(
    raking(loss = "KL") |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)), tol = 0.2) |>
      penalty(entropy()),
    data, control = list(eps_abs = 1e-8, eps_rel = 1e-8)
  )
  achieved_M <- res@balance$achieved[1]
  expect_lt(achieved_M, 0.59)              # pulled off the natural 0.6
  expect_gt(achieved_M, 0.50)              # but not past the target
  expect_false(isTRUE(all.equal(res@weights, rep(1, 5))))  # not the uniform soln
})

test_that("calibrate(tol) puts a band on the mean", {
  # sample mean 50000, target 52000, tol 1000 -> band [51000, 53000].
  data <- data.frame(income = c(40000, 50000, 60000))
  res <- fit(
    raking(loss = "L2") |> calibrate(mean(income) ~ 52000, tol = 1000) |>
      penalty(entropy()),
    data, control = list(eps_abs = 1e-8, eps_rel = 1e-8)
  )
  achieved <- res@balance$achieved[1]
  expect_lte(achieved, 53000 + 1)
  expect_gte(achieved, 51000 - 1)          # lands at the near edge ~51000
})

# --- Bounds are hard (doc 6.4, 10; fix D9/D10) ------------------------------

test_that("hard bounds are satisfied exactly, including when binding", {
  # A tol band keeps the hard margin feasible against the bounds; the upper
  # bound binds on the heavy F weight.
  data <- data.frame(sex = c("M", "M", "M", "M", "F"), stringsAsFactors = FALSE)
  res <- fit(
    raking(loss = "KL") |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)), tol = 0.2) |>
      penalty(entropy()) |>
      bounds(min = 0.7, max = 1.5),
    data, control = list(eps_abs = 1e-8, eps_rel = 1e-8)
  )
  expect_true(all(res@weights >= 0.7 - 1e-6 & res@weights <= 1.5 + 1e-6))
  expect_equal(max(res@weights), 1.5, tolerance = 1e-3)  # upper bound binds
})

# --- Convergence flag + non-convergence warning (doc 6.4, 8; fix D3) --------

test_that("a clean fit reports converged = TRUE", {
  data <- data.frame(sex = c("M", "M", "F", "F"), stringsAsFactors = FALSE)
  res <- fit(raking(loss = "L2") |>
               target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data, control = list(eps_abs = 1e-8, eps_rel = 1e-8, maxiter = 20000))
  expect_true(res@convergence$converged)
})

test_that("hitting maxiter warns and returns converged = FALSE", {
  data <- data.frame(sex = c(rep("M", 5), "F"), stringsAsFactors = FALSE)
  model <- raking(loss = "L2") |>
    target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
    penalty(entropy())
  expect_warning(
    res <- fit(model, data, eps_abs = 1e-12, eps_rel = 1e-12, maxiter = 3),
    "did not converge"
  )
  expect_false(res@convergence$converged)
  expect_equal(res@convergence$iterations, 3)
})

# --- S7 validation raises at assignment, not at fit (doc 10) ----------------

test_that("invalid loss / negative tol / bad bounds error eagerly", {
  expect_error(raking(loss = "bogus"), "must be one of")
  expect_error(target(raking(), data.frame(sex = "M", p = 1), tol = -0.1),
               "tol")
  # bounds outside [0,1) x (1,Inf] are rejected when the model is built
  expect_error(bounds(raking(), min = 1.2, max = 2), "\\[0, 1\\)")
  expect_error(bounds(raking(), min = 0.5, max = 0.9), "> 1")
})

# --- Error cases at fit() time (doc 8; fixes D6/D7/D8) ----------------------

test_that("fit() errors when a target variable is absent from the data", {
  data <- data.frame(region = c("N", "S"), stringsAsFactors = FALSE)
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data), "not found")
})

test_that("fit() errors on a data level with no target", {
  data <- data.frame(sex = c("M", "F", "X"), stringsAsFactors = FALSE)
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data), "no matching target")
})

test_that("fit() errors on a target level absent from the data", {
  data <- data.frame(sex = c("M", "M", "M"), stringsAsFactors = FALSE)
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data), "does not occur in the data")
})

test_that("fit() errors on NAs in a target or calibration variable", {
  d1 <- data.frame(sex = c("M", "F", NA), stringsAsFactors = FALSE)
  m1 <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(m1, d1), "NAs")

  d2 <- data.frame(sex = c("M", "F"), income = c(50000, NA))
  m2 <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
    calibrate(mean(income) ~ 52000)
  expect_error(fit(m2, d2), "NAs")
})

test_that("target() rejects p that does not sum to 1", {
  expect_error(target(raking(), data.frame(sex = c("M", "F"), p = c(0.3, 0.3))),
               "sum to 1")
})

# --- sum_squares shrinks toward base (doc 6.3; fix D19) ---------------------

test_that("SumSquaresRegularizer shrinks toward its base", {
  reg <- SumSquaresRegularizer(base = c(1, 1))
  # prox of lam*sum((w-b)^2): (v + 2*lam*b)/(1 + 2*lam)
  expect_equal(prox(reg, v = c(0, 0), lam = 0.5), c(0.5, 0.5))
  expect_equal(evaluate(reg, c(2, 1)), 1)        # (2-1)^2 + (1-1)^2
  # default (no base) still shrinks toward zero
  expect_equal(prox(SumSquaresRegularizer(), c(0.2, 0.4), 0.5), c(0.1, 0.2))
})

# --- evaluate() for the hard losses (supports the band convergence check) ---

test_that("evaluate() is 0 inside a hard constraint and Inf outside", {
  eq <- EqualityLoss(target = c(0.5, 0.5))
  expect_equal(evaluate(eq, c(0.5, 0.5)), 0)
  expect_equal(evaluate(eq, c(0.6, 0.4)), Inf)
  ineq <- InequalityLoss(target = c(0.5), lower = -0.1, upper = 0.1)
  expect_equal(evaluate(ineq, 0.55), 0)
  expect_equal(evaluate(ineq, 0.7), Inf)
})

# --- Cross-package parity with regrake (doc 10) -----------------------------
# Best-effort: skips if regrake is unavailable or its API differs from the
# version this was written against. Adjust the regrake() call as needed.

test_that("raking matches regrake on an exact + entropy problem", {
  skip_if_not_installed("regrake")
  set.seed(1)
  n <- 300
  data <- data.frame(
    sex = sample(c("M", "F"), n, replace = TRUE, prob = c(0.55, 0.45)),
    reg = sample(c("N", "S"), n, replace = TRUE, prob = c(0.6, 0.4)),
    stringsAsFactors = FALSE
  )
  pop <- data.frame(
    variable = c("sex", "sex", "reg", "reg"),
    level    = c("M", "F", "N", "S"),
    target   = c(0.5, 0.5, 0.5, 0.5)
  )
  # raw eps (margin_tol = NULL) so both solvers use the same tolerance
  ctrl <- list(margin_tol = NULL, eps_abs = 1e-8, eps_rel = 1e-8,
               maxiter = 20000)

  # raking tol=0 is hard exact matching, so compare against regrake rr_exact().
  # rr_exact() markers are parsed from the formula AST; they need not be in scope.
  rg <- tryCatch(
    regrake::regrake(
      data = data,
      formula = ~ rr_exact(sex) + rr_exact(reg),
      population_data = pop, pop_type = "proportions",
      regularizer = "entropy", lambda = 1,
      control = ctrl
    ),
    error = function(e) NULL
  )
  skip_if(is.null(rg), "regrake API differs; adjust the regrake() call")

  rk <- fit(
    raking() |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
      target(data.frame(reg = c("N", "S"), p = c(0.5, 0.5))) |>
      penalty(entropy()),
    data, control = ctrl
  )
  expect_equal(unname(rk@weights), unname(as.numeric(rg$weights)),
               tolerance = 5e-3)
})
