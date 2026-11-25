test_that("the toy table is read", {
  ab <- read_ab_table(toy_file())
  expect_true(is.matrix(ab))
  expect_equal(dim(ab), c(1990L, 14L))
  expect_equal(colnames(ab)[1:2], c("sample1", "sample2"))
  expect_equal(rownames(ab)[1:2], c("1", "2"))
})

test_that("invalid tables are rejected", {
  prep <- CroCoDeELR:::.prepare_ab_table
  m <- matrix(c(1, 2, 0, 3), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_error(prep(replace(m, 1, -1)), "negative", class = "crocodeel_input_error")
  expect_error(prep(replace(m, 1, NA)), "missing", class = "crocodeel_input_error")
  expect_error(prep(cbind(m, s3 = 0)), "no non-zero", class = "crocodeel_input_error")
  expect_error(prep(m[, c(1, 1)]), "duplicated", class = "crocodeel_input_error")

  f <- tempfile(fileext = ".tsv")
  writeLines(c("id\ts1\ts2", "a\t1\tx", "b\t2\t3"), f)
  expect_error(read_ab_table(f), "not numeric", class = "crocodeel_input_error")
})

test_that("tables are normalized and log-transformed", {
  m <- matrix(c(1, 3, 0, 2, 2, 4), 3,
              dimnames = list(c("a", "b", "c"), c("s1", "s2")))
  out <- CroCoDeELR:::.prepare_ab_table(m)
  expect_equal(out[, "s1"], log10(c(a = 0.25, b = 0.75, c = 0)))
  expect_equal(sum(10^out[, "s2"]), 1)
})

test_that("data.frames with a species column are accepted", {
  df <- data.frame(species = c("a", "b"), s1 = c(1, 2), s2 = c(3, 4))
  out <- CroCoDeELR:::.prepare_ab_table(df)
  expect_equal(rownames(out), c("a", "b"))
  expect_equal(colnames(out), c("s1", "s2"))
})

test_that("low-abundance filtering follows --filter-low-ab", {
  m <- matrix(c(1, 5, 30, 100), 4, dimnames = list(letters[1:4], "s1"))
  out <- CroCoDeELR:::.filter_low_ab(m, 5)
  expect_equal(unname(out[, 1]), c(0, 0, 30, 100))
})

test_that("two tables are joined like pandas (sorted union, -Inf gaps)", {
  a <- log10(matrix(c(0.5, 0.5), 2, dimnames = list(c("b", "a"), "s1")))
  b <- log10(matrix(c(0.5, 0.5), 2, dimnames = list(c("c", "B"), "s2")))
  j <- CroCoDeELR:::.join_ab_tables(a, b)
  expect_equal(rownames(j), c("B", "a", "b", "c"))
  expect_equal(unname(j[, "s1"]), c(-Inf, log10(0.5), log10(0.5), -Inf))
})
