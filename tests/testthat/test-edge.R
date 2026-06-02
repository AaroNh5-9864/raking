# Edge-case and robustness tests. Deterministic values are forced by exact
# constraints or were verified on the Python port; the rest assert structure,
# ordering, and error behavior.

ctrl <- list(eps_abs = 1e-8, eps_rel = 1e-8)

# --- Degenerate sizes -------------------------------------------------------

test_that("single-row data with a single-level margin gives weight 1", {
  res <- fit(raking() |> target(data.frame(sex = "M", p = 1)),
             data.frame(sex = "M"), control = ctrl)
  expect_equal(res@weights, 1, tolerance = 1e-5)
})

test_that("a single-level margin (p = 1) leaves weights uniform", {
  data <- data.frame(sex = c("M", "M", "M"), stringsAsFactors = FALSE)
  res <- fit(raking() |> target(data.frame(sex = "M", p = 1)), data,
             control = ctrl)
  expect_equal(res@weights, c(1, 1, 1), tolerance = 1e-4)
})

# --- Variable types and data hygiene ----------------------------------------

test_that("factor variables in data and targets work", {
  data <- data.frame(sex = factor(c("M", "F", "M", "F")))
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data, control = ctrl)
  expect_equal(sum(res@weights), 4, tolerance = 1e-6)
  expect_true(all(res@balance$satisfied))
})

test_that("unused columns and NAs outside the constraints are ignored", {
  data <- data.frame(
    sex = c("M", "F", "M", "F"),
    junk = c(1, 2, 3, 4),
    notes = c("a", NA, "c", NA),     # NA in a non-constraint column is fine
    stringsAsFactors = FALSE
  )
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data, control = ctrl)
  expect_length(res@weights, 4)
})

# --- p / n handling ---------------------------------------------------------

test_that("n counts are normalized to proportions", {
  data <- data.frame(sex = c("M", "M", "F", "F", "F"), stringsAsFactors = FALSE)
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), n = c(2, 3))),
             data, control = ctrl)
  expect_equal(res@balance$target[res@balance$level == "M"], 0.4, tolerance = 1e-9)
  expect_equal(res@balance$target[res@balance$level == "F"], 0.6, tolerance = 1e-9)
})

test_that("p within floating tolerance of 1 is accepted; far off errors", {
  expect_silent(
    target(raking(), data.frame(sex = c("M", "F"), p = c(0.5, 0.5 + 1e-8)))
  )
  expect_error(
    target(raking(), data.frame(sex = c("M", "F"), p = c(0.5, 0.55))),
    "sum to 1"
  )
})

# --- Multiple constraints ---------------------------------------------------

test_that("two margins are both matched exactly", {
  data <- data.frame(
    sex = c("M", "M", "F", "F", "M", "F"),
    reg = c("N", "S", "N", "S", "N", "S"),
    stringsAsFactors = FALSE
  )
  res <- fit(
    raking() |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
      target(data.frame(reg = c("N", "S"), p = c(0.5, 0.5))),
    data, control = ctrl
  )
  expect_true(all(res@balance$satisfied))
  expect_equal(sum(res@weights), 6, tolerance = 1e-6)
})

test_that("a margin and a mean calibration combine", {
  data <- data.frame(
    sex = c("M", "F", "M", "F"),
    income = c(40000, 60000, 50000, 70000),
    stringsAsFactors = FALSE
  )
  res <- fit(
    raking() |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
      calibrate(mean(income) ~ 55000),
    data, control = ctrl
  )
  expect_equal(res@balance$achieved[res@balance$variable == "income"], 55000,
               tolerance = 1)
  expect_equal(sum(res@weights), 4, tolerance = 1e-6)
})

# --- Model reuse across datasets --------------------------------------------

test_that("the same model fits multiple datasets independently", {
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  d1 <- data.frame(sex = c("M", "M", "F"))
  d2 <- data.frame(sex = c("M", "F", "F", "F"))
  r1 <- fit(model, d1, control = ctrl)
  r2 <- fit(model, d2, control = ctrl)
  expect_length(r1@weights, 3)
  expect_length(r2@weights, 4)
  expect_length(model@targets, 1)   # model unchanged by fitting
})

# --- Weights: ordering and invariants ---------------------------------------

test_that("weights are in input row order and respect group structure", {
  # 3 M / 1 F, target 0.5/0.5 -> each M = 2/3, F = 2.0.
  data <- data.frame(sex = c("M", "M", "M", "F"), stringsAsFactors = FALSE)
  tgt  <- data.frame(sex = c("M", "F"), p = c(0.5, 0.5))
  res  <- fit(raking() |> target(tgt), data, control = ctrl)
  expect_equal(res@weights, c(2/3, 2/3, 2/3, 2.0), tolerance = 1e-4)

  # reorder the rows; weights must follow the rows
  data2 <- data[c(4, 1, 2, 3), , drop = FALSE]
  res2  <- fit(raking() |> target(tgt), data2, control = ctrl)
  expect_equal(res2@weights, c(2.0, 2/3, 2/3, 2/3), tolerance = 1e-4)
})

# --- Penalties --------------------------------------------------------------

test_that("kl() and sum_squares() penalties run end to end", {
  data <- data.frame(sex = c("M", "M", "M", "F", "F"), stringsAsFactors = FALSE)
  tgt  <- data.frame(sex = c("M", "F"), p = c(0.5, 0.5))
  r_kl <- fit(raking() |> target(tgt) |> penalty(kl(prior = rep(1/5, 5))),
              data, control = ctrl)
  r_ss <- fit(raking() |> target(tgt) |> penalty(sum_squares()),
              data, control = ctrl)
  expect_equal(sum(r_kl@weights), 5, tolerance = 1e-6)
  expect_equal(sum(r_ss@weights), 5, tolerance = 1e-6)
})

test_that("lambda = 0 is accepted (explicitly unpenalized)", {
  data <- data.frame(sex = c("M", "M", "F"))
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
               penalty(entropy(), lambda = 0), data, control = ctrl)
  expect_equal(sum(res@weights), 3, tolerance = 1e-6)
})

# --- Bounds edge behavior ---------------------------------------------------

test_that("default bounds leave only the non-negativity requirement", {
  data <- data.frame(sex = c("M", "M", "M", "M", "F"), stringsAsFactors = FALSE)
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
               bounds(), data, control = ctrl)
  expect_true(all(res@weights >= 0))
  expect_equal(sum(res@weights), 5, tolerance = 1e-6)
})

# --- tol behavior at the extremes -------------------------------------------

test_that("a non-binding band gives the pure-soft solution regardless of width", {
  # natural M share 0.6 is inside both bands, so neither clips: identical result.
  data <- data.frame(sex = c("M", "M", "M", "F", "F"), stringsAsFactors = FALSE)
  tgt  <- data.frame(sex = c("M", "F"), p = c(0.5, 0.5))
  r2 <- fit(raking(loss = "KL") |> target(tgt, tol = 0.2) |> penalty(entropy()),
            data, control = ctrl)
  r4 <- fit(raking(loss = "KL") |> target(tgt, tol = 0.4) |> penalty(entropy()),
            data, control = ctrl)
  expect_equal(r2@weights, r4@weights, tolerance = 1e-5)
})

# --- chi_square requires positive targets (latent NaN guard) ----------------

test_that("chi_square with a zero target errors instead of producing NaN", {
  data <- data.frame(sex = c("M", "M", "F"), stringsAsFactors = FALSE)
  expect_error(
    fit(raking(loss = "chi_square") |>
          target(data.frame(sex = c("M", "F"), p = c(1, 0)), tol = 0.1),
        data, control = ctrl),
    "positive"
  )
})

# --- S7 validation at construction ------------------------------------------

test_that("invalid model/loss/regularizer arguments error eagerly", {
  expect_error(raking(name = c("a", "b")), "single string")
  expect_error(KLRegularizer(prior = c(0, 1)), "positive")
  expect_error(BandedLoss(inner = L1Loss(target = 0.5), target = 0.5,
                          lower = 0.6, upper = 0.4), "lower")
  expect_error(bounds(raking(), min = 0.5, max = 0.9), "> 1")
  expect_error(bounds(raking(), min = 1.2, max = 2), "\\[0, 1\\)")
  expect_error(bounds(raking(), min = -0.1, max = 2), "\\[0, 1\\)")
})

# --- Graceful failure on degenerate / malformed input -----------------------

test_that("non-positive maxiter and rho error clearly", {
  data <- data.frame(sex = c("M", "M", "F"), stringsAsFactors = FALSE)
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data, maxiter = 0), "maxiter")
  expect_error(fit(model, data, rho = 0), "rho")
})

test_that("zero-row data errors clearly", {
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data.frame(sex = character(0))), "no rows")
})

test_that("calibrate target of zero, and a non-mean LHS, error", {
  data <- data.frame(income = c(40000, 60000))
  expect_error(
    fit(raking() |> calibrate(mean(income) ~ 0), data),
    "non-zero"
  )
  expect_error(
    fit(raking() |> calibrate(var(income) ~ 1), data),
    "mean"
  )
})

test_that("malformed base weights error", {
  data <- data.frame(sex = c("M", "F", "M", "F"), grp = c("a", "b", "a", "b"),
                     stringsAsFactors = FALSE)
  model <- raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5)))
  expect_error(fit(model, data, base = c(1, 2, 3)), "length")        # wrong length
  expect_error(fit(model, data, base = "grp"), "numeric")            # non-numeric column
})

test_that("duplicate or missing target levels error", {
  # duplicate level in the target df -> a data row matches two cells
  expect_error(
    fit(raking() |> target(data.frame(sex = c("M", "M"), p = c(0.5, 0.5))),
        data.frame(sex = c("M", "M"))),
    "overlap"
  )
  # joint margin cell absent from the data
  joint <- data.frame(sex = c("M", "M", "F", "F"), age = c("Y", "O", "Y", "O"),
                      n = c(1, 1, 1, 1))
  data <- data.frame(sex = c("M", "M", "F"), age = c("Y", "O", "Y"),
                     stringsAsFactors = FALSE)  # no F:O
  expect_error(fit(raking() |> target(joint), data), "does not occur")
})
