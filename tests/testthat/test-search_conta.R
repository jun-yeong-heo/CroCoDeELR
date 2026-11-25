test_that("with the Python RANSAC draws, results equal CroCoDeEL 1.3.1", {
  local_numpy_subsets()
  ab <- read_ab_table(toy_file())
  events <- search_conta(ab, verbose = FALSE)
  out <- tempfile(fileext = ".tsv")
  write_conta_events(events, out)
  upstream <- system.file("extdata", "contamination_events_upstream.tsv",
                          package = "CroCoDeELR")
  expect_identical(readLines(out), readLines(upstream))
})

test_that("two-table mode equals CroCoDeEL 1.3.1 with the same draws", {
  local_numpy_subsets()
  t1 <- read_ab_table(test_path("fixtures", "toy_table1.tsv"))
  t2 <- read_ab_table(test_path("fixtures", "toy_table2.tsv"))
  events <- suppressWarnings(search_conta(t1, t2, verbose = FALSE))
  out <- tempfile(fileext = ".tsv")
  write_conta_events(events, out)
  expect_identical(readLines(out),
                   readLines(test_path("fixtures", "two_tables_upstream.tsv")))
  expect_warning(search_conta(t1, t2, verbose = FALSE), "in common")
  expect_error(search_conta(t1, t1, verbose = FALSE), "both abundance tables",
               class = "crocodeel_input_error")
})

test_that("the R port reproduces its own expected toy results", {
  expect_true(suppressMessages(self_test()))
})

test_that("cutoffs filter events", {
  ab <- read_ab_table(toy_file())
  all_pairs <- search_conta(ab, probability_cutoff = 0, verbose = FALSE)
  default <- search_conta(ab, verbose = FALSE)
  expect_gt(nrow(all_pairs), nrow(default))
  high_rate <- search_conta(ab, rate_cutoff = 0.1, verbose = FALSE)
  expect_true(all(high_rate$rate > 0.1))
  expect_true(all(diff(all_pairs$rate) <= 0))
})

test_that("parallel search gives the same result", {
  skip_on_cran()
  ab <- read_ab_table(toy_file())
  expect_equal(search_conta(ab, nproc = 2, verbose = FALSE),
               search_conta(ab, nproc = 1, verbose = FALSE))
})
