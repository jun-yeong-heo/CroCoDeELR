test_that("RANSAC recovers a unit-slope line with outliers", {
  x <- seq(-6, -2, length.out = 40)
  y <- x + 1.5
  y[c(3, 10, 25)] <- y[c(3, 10, 25)] + c(1, -2, 3)
  fit <- CroCoDeELR:::.fit_unit_slope_ransac(x, y, 30L, 0.2, 42L)
  expect_equal(fit$intercept, 1.5)
  expect_equal(sum(fit$inlier_mask), 37L)
  expect_lte(fit$num_trials, 30L)
})

test_that("RANSAC subsets are reproducible", {
  s1 <- CroCoDeELR:::.candidate_subsets(57L, 30L, 2L, 999L)
  rm(list = "57_30_2_999", envir = CroCoDeELR:::.ransac_cache)
  s2 <- CroCoDeELR:::.candidate_subsets(57L, 30L, 2L, 999L)
  expect_identical(s1, s2)
  expect_equal(dim(s1), c(30L, 2L))
  expect_true(all(s1[, 1] != s1[, 2]))
})

test_that("the stopping rule matches sklearn", {
  f <- CroCoDeELR:::.dynamic_max_trials
  expect_equal(f(10, 10, 2, 0.99), 1)
  expect_equal(f(0, 10, 2, 0.99), Inf)
  expect_equal(f(5, 10, 2, 0.99), 17)
  expect_equal(f(5, 10, 2, 0), 0)
})
