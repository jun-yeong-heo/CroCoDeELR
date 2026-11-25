# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of the easy_wf command of crocodeel/crocodeel.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

.decision_support_warning <- function() {
  message(
    "CroCoDeEL is a decision-support tool and should not be considered a ",
    "definitive contamination classification system. Reported contamination ",
    "events may include false positives, particularly for samples with ",
    "similar species abundance profiles (e.g. longitudinal samples). ",
    "Manually review the scatterplots of each predicted event. The CroCoDeEL ",
    "Interpretation Interface can help: ",
    "https://metagenopolis.github.io/CroCoDeEL_interpreter/"
  )
}

#' Run the complete CroCoDeEL workflow on files
#'
#' Port of `crocodeel easy_wf`: reads the abundance table(s), searches for
#' contamination, writes the events to a TSV file and the scatterplots to a
#' PDF report.
#'
#' @inheritParams search_conta
#' @inheritParams plot_conta
#' @param ab_file Path to the species abundance table (TSV).
#' @param events_file Path of the output TSV file.
#' @param pdf_file Path of the output PDF report.
#' @param ab_file_2 Optional path to a second abundance table.
#' @param model_file Optional path to a model saved by [train_model()] or
#'   [import_sklearn_model()]. `NULL` uses the default model.
#' @return The events data.frame, invisibly.
#' @examples
#' \donttest{
#' out <- tempfile()
#' easy_wf(system.file("extdata", "mgs_profiles_test.tsv",
#'                     package = "CroCoDeELR"),
#'         events_file = paste0(out, ".tsv"), pdf_file = paste0(out, ".pdf"))
#' }
#' @export
easy_wf <- function(ab_file, events_file = "contamination_events.tsv",
                    pdf_file = "contamination_events.pdf", ab_file_2 = NULL,
                    filter_low_ab = NULL, model_file = NULL,
                    probability_cutoff = 0.5, rate_cutoff = 0, nproc = 1L,
                    nrow = 4L, ncol = 4L, conta_line = TRUE,
                    color_conta_species = FALSE, verbose = TRUE) {
  ab <- read_ab_table(ab_file)
  ab2 <- NULL
  if (!is.null(ab_file_2)) {
    if (normalizePath(ab_file_2, mustWork = FALSE) ==
        normalizePath(ab_file, mustWork = FALSE)) {
      warning("Both abundance tables refer to the same file; the second one ",
              "is ignored and contamination is searched within a single ",
              "table.", call. = FALSE)
      ab_file_2 <- NULL
    } else {
      ab2 <- read_ab_table(ab_file_2)
    }
  }

  model <- load_model(model_file)
  model_label <- if (is.null(model_file)) {
    system.file("extdata", .DEFAULT_MODEL_FILE, package = "CroCoDeELR")
  } else {
    normalizePath(model_file)
  }

  events <- search_conta(ab, ab2, filter_low_ab = filter_low_ab,
                         model = model,
                         probability_cutoff = probability_cutoff,
                         rate_cutoff = rate_cutoff, nproc = nproc,
                         verbose = verbose)

  header <- .execution_description(
    normalizePath(ab_file),
    if (is.null(ab_file_2)) NULL else normalizePath(ab_file_2),
    model_label, filter_low_ab, probability_cutoff, rate_cutoff
  )
  write_conta_events(events, events_file, header = header)
  .log(verbose, "Contamination events saved in ", events_file)

  if (nrow(events) > 0L) {
    plot_conta(ab, events, file = pdf_file, ab_table_2 = ab2,
               filter_low_ab = filter_low_ab, nrow = nrow, ncol = ncol,
               conta_line = conta_line,
               color_conta_species = color_conta_species, verbose = verbose)
    .log(verbose, "PDF report saved in ", pdf_file)
    .decision_support_warning()
  } else {
    .log(verbose, "No PDF report generated.")
  }

  invisible(events)
}
