# Larger problems, extreme magnitudes, and alternative data-frame inputs.
# These assert invariants (convergence, sum = n, satisfied margins, positivity)
# rather than baked values, so they hold for any feasible random draw.

test_that("end-to-end at n = 500: two binary margins + a mean (doc Section 10)", {
  set.seed(42)
  n <- 500
  data <- data.frame(
    sex    = sample(c("M", "F"), n, replace = TRUE, prob = c(0.55, 0.45)),
    reg    = sample(c("N", "S"), n, replace = TRUE, prob = c(0.60, 0.40)),
    income = rnorm(n, 50000, 12000),
    stringsAsFactors = FALSE
  )
  res <- fit(
    raking() |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
      target(data.frame(reg = c("N", "S"), p = c(0.5, 0.5))) |>
      calibrate(mean(income) ~ 52000),
    data   # default control (margin_tol = 1e-4); converges in ~85 iterations
  )
  expect_true(res@convergence$converged)
  expect_length(res@weights, n)
  expect_equal(sum(res@weights), n, tolerance = 1e-6)
  expect_true(all(res@weights > 0))
  expect_true(all(res@balance$satisfied))
  # exact constraints match to ~eps, well inside 1e-4 for the proportions
  expect_equal(res@balance$achieved[res@balance$level == "M"], 0.5,
               tolerance = 1e-4)
  expect_equal(res@balance$achieved[res@balance$level == "N"], 0.5,
               tolerance = 1e-4)
})

test_that("large-magnitude calibration is handled by normalization", {
  set.seed(7)
  n <- 300
  data <- data.frame(
    sex    = sample(c("M", "F"), n, replace = TRUE, prob = c(0.55, 0.45)),
    income = rnorm(n, 1e6, 2.5e5)
  )
  res <- fit(
    raking() |>
      target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
      calibrate(mean(income) ~ 1.05e6),
    data
  )
  expect_true(res@convergence$converged)
  expect_equal(sum(res@weights), n, tolerance = 1e-6)
  expect_true(all(res@balance$satisfied))
  achieved <- res@balance$achieved[res@balance$variable == "income"]
  expect_lt(abs(achieved - 1.05e6) / 1.05e6, 1e-4)   # relative error tiny
})

test_that("a penalized band at n = 500 stays feasible and binding-aware", {
  set.seed(99)
  n <- 500
  data <- data.frame(grp = sample(c("A", "B", "C"), n, replace = TRUE,
                                  prob = c(0.5, 0.3, 0.2)),
                     stringsAsFactors = FALSE)
  res <- fit(
    raking(loss = "KL") |>
      target(data.frame(grp = c("A", "B", "C"), p = c(1/3, 1/3, 1/3)),
             tol = 0.02) |>
      penalty(entropy()),
    data, control = list(eps_abs = 1e-7, eps_rel = 1e-7)
  )
  expect_true(res@convergence$converged)
  expect_equal(sum(res@weights), n, tolerance = 1e-6)
  expect_true(all(res@balance$satisfied))
  # every achieved proportion lands within the +/- 0.02 band (up to a tiny
  # convergence-residual slack)
  expect_true(all(abs(res@balance$difference) <= 0.02 + 1e-3))
})

# --- Alternative data-frame inputs ------------------------------------------

test_that("fit() works on a tibble", {
  skip_if_not_installed("tibble")
  data <- tibble::tibble(sex = c("M", "F", "M", "F"), income = c(40, 60, 50, 70))
  res <- fit(
    raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
      calibrate(mean(income) ~ 55),
    data, control = list(eps_abs = 1e-8, eps_rel = 1e-8)
  )
  expect_length(res@weights, 4)
  expect_equal(sum(res@weights), 4, tolerance = 1e-6)
  expect_true(all(res@balance$satisfied))
})

test_that("fit() works on a data.table", {
  skip_if_not_installed("data.table")
  data <- data.table::data.table(sex = c("M", "F", "M", "F"))
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data, control = list(eps_abs = 1e-8, eps_rel = 1e-8))
  expect_length(res@weights, 4)
  expect_equal(sum(res@weights), 4, tolerance = 1e-6)
})
