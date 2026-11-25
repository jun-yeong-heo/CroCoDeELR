test_that("floats are formatted like Python", {
  f <- CroCoDeELR:::.py_float_str
  expect_equal(f(c(1, 0.87, 0.1 + 0.2, 1e-5, 100, 0.00025)),
               c("1.0", "0.87", "0.30000000000000004", "1e-05", "100.0",
                 "0.00025"))
  expect_equal(CroCoDeELR:::.fmt_rate(c(0.741474, 0.0096)),
               c("7.41e-01", "9.60e-03"))
  expect_equal(CroCoDeELR:::.normalize_exponent("1.00e-003"), "1.00e-03")
})

test_that("rates are rounded to three significant digits", {
  r <- CroCoDeELR:::.round_conta_rate
  expect_equal(r(74.1474), 74.1)
  expect_equal(r(0.96354), 0.964)
  expect_error(r(0))
})

test_that("the caller's random number stream is left untouched", {
  set.seed(123)
  expected <- runif(3)
  set.seed(123)
  CroCoDeELR:::.with_seed(42, runif(10))
  expect_equal(runif(3), expected)
})
