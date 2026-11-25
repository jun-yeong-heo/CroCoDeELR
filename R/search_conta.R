# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/search_conta.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

.SEARCH_CHUNKSIZE <- 100L

# Classify one batch of sample pairs (columns `source`, `target`).
.classify_sample_pairs <- function(pairs, ab, model) {
  feats <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    if (pairs$source[i] != pairs$target[i]) {
      feats[i] <- list(.extract_features(ab, pairs$source[i], pairs$target[i]))
    }
  }
  keep <- !vapply(feats, is.null, logical(1))
  if (!any(keep)) {
    return(NULL)
  }
  feats <- feats[keep]
  pairs <- pairs[keep, , drop = FALSE]

  feature_matrix <- matrix(
    unlist(lapply(feats, `[[`, "values"), use.names = FALSE),
    ncol = .NUM_FEATURES, byrow = TRUE
  )
  probability <- .predict_conta_proba(model, feature_matrix)
  offsets <- vapply(feats, `[[`, numeric(1), "conta_line_offset")

  data.frame(
    source = pairs$source,
    target = pairs$target,
    rate = 10^(-offsets),
    probability = probability,
    contamination_specific_species = vapply(
      feats, function(f) paste(f$conta_line_species, collapse = ","),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

.empty_events <- function() {
  data.frame(
    source = character(0),
    target = character(0),
    rate = numeric(0),
    probability = numeric(0),
    contamination_specific_species = character(0),
    stringsAsFactors = FALSE
  )
}

# Sort events by decreasing rate, then source and target (code point order).
.sort_events <- function(events) {
  if (nrow(events) == 0L) {
    return(events)
  }
  events <- events[.c_order(-events$rate, events$source, events$target), ,
                   drop = FALSE]
  rownames(events) <- NULL
  events
}

.check_cutoff <- function(x, name, upper = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 ||
      x > upper) {
    .input_error("`", name, "` must be a finite number between 0 and ",
                 upper, ".")
  }
}

#' Search for cross-sample contamination
#'
#' Port of `crocodeel search_conta`. Every ordered pair of samples is tested:
#' a contamination line is searched in the scatterplot of species abundances
#' (target on x, source on y), its features are computed and the Random
#' Forest model gives the probability of contamination.
#'
#' Unlike the original CroCoDeEL, RANSAC subsets are drawn with R's random
#' number generator. Results are therefore close to, but not always identical
#' with, those of the Python version. See the README for the measured
#' differences.
#'
#' @param ab_table Species abundance table: a numeric matrix or data.frame
#'   with species in rows (row names = species names) and samples in columns,
#'   for example the output of [read_ab_table()]. Abundances may be counts or
#'   relative abundances: each sample is normalized to sum to one.
#' @param ab_table_2 Optional second table. When given, only contamination
#'   from samples of `ab_table` (sources) to samples of `ab_table_2`
#'   (targets) is searched.
#' @param filter_low_ab Optional factor. In each sample, abundances lower than
#'   or equal to `filter_low_ab` times the smallest non-zero abundance are set
#'   to zero (option `--filter-low-ab` of CroCoDeEL).
#' @param model A `crocodeel_rf` model. `NULL` uses the default model.
#' @param probability_cutoff Minimum probability to report an event.
#' @param rate_cutoff Events with a rate lower than or equal to this value are
#'   not reported.
#' @param nproc Number of R processes. Values above 1 use a PSOCK cluster,
#'   which also works on Windows.
#' @param verbose Print progress messages.
#' @return A data.frame with one row per contamination event and columns
#'   `source`, `target`, `rate`, `probability` and
#'   `contamination_specific_species` (comma-separated species names), sorted
#'   by decreasing rate.
#' @examples
#' ab <- read_ab_table(system.file("extdata", "mgs_profiles_test.tsv",
#'                                 package = "CroCoDeELR"))
#' events <- search_conta(ab)
#' events[, c("source", "target", "rate", "probability")]
#' @export
search_conta <- function(ab_table, ab_table_2 = NULL, filter_low_ab = NULL,
                         model = NULL, probability_cutoff = 0.5,
                         rate_cutoff = 0, nproc = 1L,
                         verbose = interactive()) {
  .check_cutoff(probability_cutoff, "probability_cutoff")
  .check_cutoff(rate_cutoff, "rate_cutoff")
  if (!is.numeric(nproc) || length(nproc) != 1L || nproc < 1) {
    .input_error("`nproc` must be a positive integer.")
  }
  model <- if (is.null(model)) load_model() else .validate_model(model)

  tables <- .prepare_tables(ab_table, ab_table_2, filter_low_ab)
  ab <- tables$table

  if (is.null(tables$samples_2)) {
    sources <- tables$samples_1
    targets <- tables$samples_1
  } else {
    sources <- tables$samples_1
    targets <- tables$samples_2
  }

  # itertools.product(sources, targets): source varies slowest.
  pairs <- data.frame(
    source = rep(sources, each = length(targets)),
    target = rep(targets, times = length(sources)),
    stringsAsFactors = FALSE
  )

  .log(verbose, "Searching contamination in ", nrow(pairs),
       " sample pairs using ", nproc, " process",
       if (nproc == 1) "" else "es", "...")
  start <- proc.time()[["elapsed"]]

  batches <- lapply(.chunk(seq_len(nrow(pairs)), .SEARCH_CHUNKSIZE),
                    function(idx) pairs[idx, , drop = FALSE])

  results <- .par_lapply(batches, .classify_sample_pairs, nproc = nproc,
                         progress = verbose && nproc <= 1, ab = ab,
                         model = model)
  results <- results[!vapply(results, is.null, logical(1))]
  events <- if (length(results)) do.call(rbind, results) else .empty_events()

  events <- events[events$probability >= probability_cutoff &
                     events$rate > rate_cutoff, , drop = FALSE]
  events <- .sort_events(events)

  .log(verbose, "Search completed in ",
       sprintf("%.1f", proc.time()[["elapsed"]] - start), " seconds")
  .log(verbose, nrow(events), " contamination event",
       if (nrow(events) == 1L) "" else "s", " detected")
  .log(verbose, length(unique(events$target)), " sample",
       if (length(unique(events$target)) == 1L) "" else "s", " contaminated")

  events
}
