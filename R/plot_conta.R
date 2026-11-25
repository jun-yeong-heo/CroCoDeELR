# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/plot_conta.py (base R graphics instead of matplotlib)
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

.check_events_samples <- function(ab, events) {
  samples <- colnames(ab)
  unknown_sources <- sort(unique(events$source[!events$source %in% samples]),
                          method = "radix")
  unknown_targets <- sort(unique(events$target[!events$target %in% samples]),
                          method = "radix")
  msgs <- character(0)
  if (length(unknown_sources)) {
    msgs <- c(msgs, paste0("unknown source samples: ",
                           .format_sample_names(unknown_sources)))
  }
  if (length(unknown_targets)) {
    msgs <- c(msgs, paste0("unknown target samples: ",
                           .format_sample_names(unknown_targets)))
  }
  if (length(msgs)) {
    .input_error(paste(msgs, collapse = "; "))
  }
}

# Draw the scatterplot of one contamination event.
.plot_event <- function(ab, pseudo_zero, event, conta_line,
                        color_conta_species) {
  target <- ab[, event$target]
  source <- ab[, event$source]
  to_plot <- target > pseudo_zero | source > pseudo_zero

  edge <- rep("black", sum(to_plot))
  if (color_conta_species) {
    conta_species <- strsplit(event$contamination_specific_species, ",",
                              fixed = TRUE)[[1]]
    edge[rownames(ab)[to_plot] %in% conta_species] <- "orange"
  }

  ticks <- seq(pseudo_zero, 0)
  lims <- range(c(ticks, target[to_plot], source[to_plot]))

  graphics::plot(
    target[to_plot], source[to_plot],
    xlim = lims, ylim = lims,
    pch = 1, cex = 0.6, col = edge,
    axes = FALSE, xlab = event$target, ylab = event$source,
    main = paste0("prob = ", .py_float_str(event$probability),
                  ", rate = ",
                  .py_float_str(.round_conta_rate(100 * event$rate)), "%"),
    cex.main = 1, font.main = 1
  )

  labels <- lapply(ticks, function(t) bquote(10^.(t)))
  labels[[1]] <- "0"
  labels[[length(labels)]] <- "1"
  labels <- as.expression(labels)
  graphics::axis(1, at = ticks, labels = labels, tcl = 0.3, cex.axis = 0.8)
  graphics::axis(2, at = ticks, labels = labels, tcl = 0.3, cex.axis = 0.8,
                 las = 1)
  graphics::box()

  graphics::abline(a = 0, b = 1, col = "grey", lwd = 0.5)
  if (conta_line) {
    graphics::abline(a = -log10(event$rate), b = 1,
                     col = grDevices::adjustcolor("red", alpha.f = 0.5),
                     lwd = 0.5)
  }
}

#' Plot contamination events to a PDF report
#'
#' Port of `crocodeel plot_conta`. Draws, for each event, the scatterplot of
#' species abundances in the target (x) and the source (y) on a log scale,
#' with the identity line (grey) and the contamination line (red). Zeros are
#' drawn at a pseudo-zero below the smallest abundance.
#'
#' @inheritParams search_conta
#' @param events Data.frame returned by [search_conta()] or
#'   [read_conta_events()].
#' @param file Path of the PDF file.
#' @param nrow,ncol Number of plots per column and per row on each page
#'   (1 to 11).
#' @param conta_line Draw the contamination line.
#' @param color_conta_species Draw species of the contamination line in
#'   orange.
#' @return `file`, invisibly. No file is written when `events` is empty.
#' @examples
#' ab <- read_ab_table(system.file("extdata", "mgs_profiles_test.tsv",
#'                                 package = "CroCoDeELR"))
#' events <- search_conta(ab)
#' pdf_file <- tempfile(fileext = ".pdf")
#' plot_conta(ab, events, file = pdf_file)
#' @export
plot_conta <- function(ab_table, events, file = "contamination_events.pdf",
                       ab_table_2 = NULL, filter_low_ab = NULL, nrow = 4L,
                       ncol = 4L, conta_line = TRUE,
                       color_conta_species = FALSE, verbose = interactive()) {
  for (arg in c("nrow", "ncol")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1L || v < 1 || v > 11 || v != round(v)) {
      .input_error("`", arg, "` must be an integer between 1 and 11.")
    }
  }
  if (nrow(events) == 0L) {
    .log(verbose, "No PDF report generated.")
    return(invisible(file))
  }

  ab <- .prepare_tables(ab_table, ab_table_2, filter_low_ab)$table
  .check_events_samples(ab, events)

  pseudo_zero <- floor(min(ab[is.finite(ab)]))
  ab[ab == -Inf] <- pseudo_zero

  start <- proc.time()[["elapsed"]]
  grDevices::pdf(file, width = 4 * ncol, height = 4 * nrow)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(mfrow = c(nrow, ncol))
  # mfrow shrinks text for large grids; restore the default size.
  graphics::par(cex = 1, mar = c(4, 4, 2, 1), mgp = c(2.2, 0.5, 0))
  on.exit(graphics::par(old_par), add = TRUE, after = FALSE)

  for (i in seq_len(nrow(events))) {
    .plot_event(ab, pseudo_zero, events[i, , drop = FALSE], conta_line,
                color_conta_species)
  }

  .log(verbose, "PDF report generated in ",
       sprintf("%.1f", proc.time()[["elapsed"]] - start), " seconds")
  invisible(file)
}
