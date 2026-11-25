test_that("event files round-trip and upstream files are read", {
  upstream <- read_conta_events(
    system.file("extdata", "contamination_events_upstream.tsv",
                package = "CroCoDeELR"))
  expect_equal(nrow(upstream), 7L)
  f <- tempfile(fileext = ".tsv")
  write_conta_events(upstream, f, header = "comment")
  expect_equal(read_conta_events(f), upstream)
  expect_equal(readLines(f)[1], "# comment")
})

test_that("malformed event files are rejected", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("source\ttarget\trate\tprobability", "a\tb\t0.1\t0.9"), f)
  expect_error(read_conta_events(f), "missing required",
               class = "crocodeel_input_error")
  writeLines(c("source\ttarget\trate\tprobability\tcontamination_specific_species",
               "a\tb\t2\t0.9\tx,y"), f)
  expect_error(read_conta_events(f), "at most 1",
               class = "crocodeel_input_error")
})

test_that("a PDF report is generated", {
  ab <- read_ab_table(toy_file())
  events <- search_conta(ab, verbose = FALSE)
  pdf_file <- tempfile(fileext = ".pdf")
  plot_conta(ab, events, file = pdf_file, nrow = 2, ncol = 2,
             color_conta_species = TRUE, verbose = FALSE)
  expect_true(file.exists(pdf_file))
  expect_gt(file.size(pdf_file), 1000)
  bad <- events
  bad$source[1] <- "unknown"
  expect_error(plot_conta(ab, bad, file = pdf_file, verbose = FALSE),
               "unknown source", class = "crocodeel_input_error")
})

test_that("easy_wf writes both outputs", {
  out <- tempfile()
  events <- suppressMessages(easy_wf(toy_file(),
                                     events_file = paste0(out, ".tsv"),
                                     pdf_file = paste0(out, ".pdf"),
                                     verbose = FALSE))
  lines <- readLines(paste0(out, ".tsv"))
  expect_match(lines[1], "^# CroCoDeELR version: ")
  expect_equal(read_conta_events(paste0(out, ".tsv")), events,
               tolerance = 1e-2)
  expect_true(file.exists(paste0(out, ".pdf")))
})
