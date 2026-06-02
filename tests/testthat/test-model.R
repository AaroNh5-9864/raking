test_that("raking() creates a model with defaults", {
  m <- raking()
  expect_true(S7::S7_inherits(m, raking_model))
  expect_equal(m@loss, "KL")
  expect_equal(m@name, "wt")
  expect_equal(m@targets, list())
  expect_null(m@penalty)
  expect_null(m@bounds)
  expect_equal(m@lambda, 1)
})

test_that("raking() accepts loss and name", {
  m <- raking(loss = "L2", name = "weight")
  expect_equal(m@loss, "L2")
  expect_equal(m@name, "weight")
})

test_that("raking() rejects an unknown loss", {
  expect_error(raking(loss = "bogus"), "must be one of")
})

test_that("raking_model rejects bad bounds", {
  expect_error(
    raking_model(loss = "KL", name = "wt", bounds = list(min = 2, max = 1)),
    "must be <"
  )
})

test_that("raking_model prints a summary", {
  expect_output(print(raking()), "raking_model")
})
