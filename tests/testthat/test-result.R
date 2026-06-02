test_that("weights() returns the fit's weights", {
  data <- data.frame(sex = c("M", "F", "M", "F"))
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data)
  expect_equal(weights(res), res@weights)
  expect_length(weights(res), 4)
})

test_that("summary() prints a balance table and returns the fit invisibly", {
  data <- data.frame(sex = c("M", "M", "F"))
  res <- fit(raking() |> target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))),
             data)
  expect_output(summary(res), "Balance")
  capture.output(out <- summary(res))
  expect_true(S7::S7_inherits(out, raking_fit))
})

test_that("kish_ess equals n for uniform weights", {
  expect_equal(kish_ess(rep(1, 10)), 10)
})
