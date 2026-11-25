# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/conta_event.py and crocodeel/execution_description.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

.EVENT_COLUMNS <- c("source", "target", "rate", "probability",
                    "contamination_specific_species")

# Round a contamination rate (percentage) to three significant digits.
.round_conta_rate <- function(rate, significant_digits = 3L) {
  if (rate <= 0) {
    stop("Contamination rate must be strictly positive.")
  }
  magnitude <- floor(log10(rate))
  decimals <- as.integer(significant_digits - (magnitude + 1))
  .np_round(rate, decimals)
}

#' Write contamination events to a TSV file
#'
#' Writes events in the format of `crocodeel search_conta`, so that the file
#' can also be read by the Python CroCoDeEL and its interpretation interface.
#'
#' @param events Data.frame returned by [search_conta()].
#' @param file Output path.
#' @param header Optional comment line written before the column names
#'   (without the leading `#`).
#' @return `file`, invisibly.
#' @export
write_conta_events <- function(events, file, header = NULL) {
  missing <- setdiff(.EVENT_COLUMNS, names(events))
  if (length(missing) > 0L) {
    .input_error("`events` is missing column(s): ",
                 paste(missing, collapse = ", "))
  }
  lines <- character(0)
  if (!is.null(header)) {
    lines <- paste0("# ", header)
  }
  lines <- c(lines, paste(.EVENT_COLUMNS, collapse = "\t"))
  if (nrow(events) > 0L) {
    lines <- c(lines, paste(
      events$source,
      events$target,
      .fmt_rate(events$rate),
      .py_float_str(events$probability),
      events$contamination_specific_species,
      sep = "\t"
    ))
  }
  con <- file(file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
  invisible(file)
}

#' Read contamination events from a TSV file
#'
#' Reads and validates a file written by [write_conta_events()] or by the
#' Python CroCoDeEL. Comment lines before the header are skipped.
#'
#' @param file Path to the TSV file.
#' @return A data.frame with the same columns as [search_conta()].
#' @export
read_conta_events <- function(file) {
  lines <- readLines(file, encoding = "UTF-8", warn = FALSE)
  is_comment <- startsWith(lines, "#")
  first <- which(!is_comment)[1]
  if (is.na(first)) {
    .input_error("Contamination events file is empty or has no header.")
  }
  header <- strsplit(lines[first], "\t", fixed = TRUE)[[1]]
  missing <- setdiff(.EVENT_COLUMNS, header)
  if (length(missing) > 0L) {
    .input_error("Contamination events file is missing required column(s): ",
                 paste(missing, collapse = ", "))
  }

  body <- lines[-seq_len(first)]
  line_numbers <- first + seq_along(body)
  keep <- nzchar(body)
  body <- body[keep]
  line_numbers <- line_numbers[keep]

  fields <- strsplit(body, "\t", fixed = TRUE)
  events <- .empty_events()
  if (length(fields) == 0L) {
    return(events)
  }

  rows <- lapply(seq_along(fields), function(i) {
    f <- fields[[i]]
    ln <- line_numbers[i]
    # strsplit() drops a trailing empty field; restore it before checking.
    if (endsWith(body[i], "\t")) {
      f <- c(f, "")
    }
    if (length(f) != length(header)) {
      .input_error("Invalid number of fields on line ", ln,
                   " of the contamination events file.")
    }
    names(f) <- header
    rate <- suppressWarnings(as.numeric(f[["rate"]]))
    probability <- suppressWarnings(as.numeric(f[["probability"]]))
    if (is.na(rate) || !is.finite(rate)) {
      .input_error("Invalid rate '", f[["rate"]], "' on line ", ln,
                   " of the contamination events file.")
    }
    if (is.na(probability) || !is.finite(probability)) {
      .input_error("Invalid probability '", f[["probability"]], "' on line ",
                   ln, " of the contamination events file.")
    }
    if (!nzchar(f[["source"]])) {
      .input_error("Missing source sample on line ", ln,
                   " of the contamination events file.")
    }
    if (!nzchar(f[["target"]])) {
      .input_error("Missing target sample on line ", ln,
                   " of the contamination events file.")
    }
    if (!(rate > 0 && rate <= 1)) {
      .input_error("Invalid contamination rate '", f[["rate"]], "' on line ",
                   ln, ": value must be greater than 0 and at most 1.")
    }
    if (!(probability >= 0 && probability <= 1)) {
      .input_error("Invalid contamination probability '", f[["probability"]],
                   "' on line ", ln, ": value must be between 0 and 1.")
    }
    if (!nzchar(f[["contamination_specific_species"]])) {
      .input_error("Missing contamination-specific species on line ", ln,
                   " of the contamination events file.")
    }
    data.frame(
      source = f[["source"]],
      target = f[["target"]],
      rate = rate,
      probability = probability,
      contamination_specific_species = f[["contamination_specific_species"]],
      stringsAsFactors = FALSE
    )
  })
  events <- do.call(rbind, rows)
  rownames(events) <- NULL
  events
}

# Comment line describing the run (execution_description.py).
.execution_description <- function(ab_file, ab_file_2, model_file,
                                   filter_low_ab, probability_cutoff,
                                   rate_cutoff) {
  fields <- c(
    "CroCoDeELR version" = as.character(utils::packageVersion("CroCoDeELR")),
    "rf_model" = model_file,
    "hostname" = Sys.info()[["nodename"]],
    "username" = Sys.info()[["user"]],
    "datetime" = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    "species_ab_table" = ab_file
  )
  if (!is.null(ab_file_2)) {
    fields <- c(fields, "species_ab_table_2" = ab_file_2)
  }
  fields <- c(
    fields,
    "filtering_ab_thr_factor" = if (is.null(filter_low_ab)) "None" else
      .py_float_str(filter_low_ab),
    "probability_cutoff" = .py_float_str(probability_cutoff),
    "rate_cutoff" = .py_float_str(rate_cutoff)
  )
  paste(paste0(names(fields), ": ", fields), collapse = " | ")
}
