# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/self_test.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

#' Run an end-to-end test of CroCoDeELR
#'
#' Runs [search_conta()] and [plot_conta()] on the toy dataset of CroCoDeEL
#' and compares the events with the results expected for this R port
#' (`contamination_events_R.tsv`). Sources, targets, species lists and
#' formatted rates must be identical; probabilities may differ by at most
#' `1e-9`.
#'
#' The expected results of the R port differ from those of the Python
#' CroCoDeEL only through the RANSAC random draws; see the README.
#'
#' @param keep_results If `TRUE`, the TSV and PDF files are kept and their
#'   paths are printed.
#' @return `TRUE`, invisibly, when the test succeeds. An error otherwise.
#' @examples
#' \donttest{
#' self_test()
#' }
#' @export
self_test <- function(keep_results = FALSE) {
  ab_file <- system.file("extdata", "mgs_profiles_test.tsv",
                         package = "CroCoDeELR")
  expected_file <- system.file("extdata", "contamination_events_R.tsv",
                               package = "CroCoDeELR")
  tsv <- tempfile("crocodeelr_test_", fileext = ".tsv")
  pdf <- tempfile("crocodeelr_test_", fileext = ".pdf")
  if (!keep_results) {
    on.exit(unlink(c(tsv, pdf)), add = TRUE)
  }

  message("Running tests on the toy dataset")
  ab <- read_ab_table(ab_file)
  events <- search_conta(ab, verbose = FALSE)
  write_conta_events(events, tsv)
  events <- read_conta_events(tsv)
  plot_conta(ab, events, file = pdf, verbose = FALSE)

  expected <- read_conta_events(expected_file)
  same_rows <- nrow(events) == nrow(expected) &&
    identical(events$source, expected$source) &&
    identical(events$target, expected$target) &&
    identical(.fmt_rate(events$rate), .fmt_rate(expected$rate)) &&
    identical(events$contamination_specific_species,
              expected$contamination_specific_species) &&
    all(abs(events$probability - expected$probability) <= 1e-9)
  if (!same_rows) {
    stop("Contamination events found are not those expected.", call. = FALSE)
  }
  message("All contamination events with expected rates and probabilities found")

  if (!file.exists(pdf) || file.size(pdf) < 1000) {
    stop("The PDF report was not generated correctly.", call. = FALSE)
  }
  message("PDF report generated")

  if (keep_results) {
    message("Results kept in ", tsv, " and ", pdf)
  }
  message("Tests completed successfully")
  invisible(TRUE)
}
