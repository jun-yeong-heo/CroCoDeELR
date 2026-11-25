# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/ab_table_utils.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

# Strings that pandas.read_csv() treats as missing values by default.
.PANDAS_NA_STRINGS <- c(
  "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
  "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a",
  "nan", "null"
)

#' Read a species abundance table
#'
#' Reads a tab-separated species abundance table. The first column holds
#' species names and every other column holds the abundances of one sample.
#' Lines and line ends starting with `#` are ignored, as in the original
#' CroCoDeEL. Abundances are returned as they are in the file; normalization
#' and log-transformation are done by [search_conta()], [plot_conta()] and
#' [train_model()].
#'
#' @param file Path to the TSV file.
#' @return A numeric matrix with species in rows and samples in columns.
#' @examples
#' ab <- read_ab_table(system.file("extdata", "mgs_profiles_test.tsv",
#'                                 package = "CroCoDeELR"))
#' dim(ab)
#' @export
read_ab_table <- function(file) {
  if (!file.exists(file)) {
    .input_error("Species abundance table '", file, "' does not exist.")
  }

  read_with <- function(col_classes, na_strings) {
    utils::read.delim(
      file,
      header = TRUE,
      sep = "\t",
      quote = "\"",
      comment.char = "#",
      colClasses = col_classes,
      na.strings = na_strings,
      check.names = FALSE,
      row.names = NULL,
      strip.white = FALSE,
      blank.lines.skip = TRUE
    )
  }

  # Fast path: numeric abundance columns. Reading everything as character
  # first needs several times more memory on large tables.
  header <- tryCatch(
    scan(file, what = "", sep = "\t", quote = "\"", comment.char = "#",
         nlines = 1L, quiet = TRUE, blank.lines.skip = TRUE),
    error = function(e) character(0)
  )
  fast <- NULL
  if (length(header) >= 2L) {
    fast <- tryCatch(
      suppressWarnings(read_with(
        c("character", rep("numeric", length(header) - 1L)),
        .PANDAS_NA_STRINGS
      )),
      error = function(e) NULL
    )
  }
  if (!is.null(fast) && ncol(fast) >= 1L &&
      all(vapply(fast[-1L], is.numeric, logical(1)))) {
    ab <- as.matrix(fast[-1L])
    dimnames(ab) <- list(fast[[1L]], names(fast)[-1L])
    return(ab)
  }

  # Slow path: read as character to report non-numeric samples precisely.
  raw <- tryCatch(
    read_with("character", character(0)),
    error = function(e) {
      if (grepl("no lines available", conditionMessage(e), fixed = TRUE)) {
        .input_error("Species abundance table '", file, "' is empty.")
      }
      .input_error("Invalid species abundance table '", file, "': ",
                   conditionMessage(e))
    }
  )

  if (ncol(raw) < 1L) {
    .input_error("Species abundance table '", file, "' is empty.")
  }

  species_names <- raw[[1L]]
  values <- raw[-1L]
  sample_names <- names(values)

  num_values <- lapply(values, function(column) {
    column[column %in% .PANDAS_NA_STRINGS] <- NA_character_
    suppressWarnings(as.numeric(column))
  })

  bad_samples <- sample_names[vapply(seq_along(values), function(i) {
    any(!is.na(values[[i]]) & !(values[[i]] %in% .PANDAS_NA_STRINGS) &
          is.na(num_values[[i]]))
  }, logical(1))]

  if (length(bad_samples) > 0L) {
    .input_error("Species abundance in the following samples is not numeric: ",
                 .format_sample_names(bad_samples))
  }

  ab <- matrix(
    unlist(num_values, use.names = FALSE),
    nrow = length(species_names),
    ncol = length(sample_names),
    dimnames = list(species_names, sample_names)
  )

  ab
}

# Convert user input (matrix or data.frame, species x samples) to a matrix.
.as_ab_matrix <- function(x, arg = "ab_table") {
  if (is.data.frame(x)) {
    if (is.null(rownames(x)) ||
        identical(rownames(x), as.character(seq_len(nrow(x))))) {
      if (ncol(x) >= 1L && !is.numeric(x[[1L]])) {
        rn <- as.character(x[[1L]])
        x <- x[-1L]
        rownames_x <- rn
      } else {
        rownames_x <- rownames(x)
      }
    } else {
      rownames_x <- rownames(x)
    }
    non_numeric <- names(x)[!vapply(x, is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      .input_error("Species abundance in the following samples is not numeric: ",
                   .format_sample_names(non_numeric))
    }
    x <- as.matrix(x)
    rownames(x) <- rownames_x
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    .input_error("`", arg, "` must be a numeric matrix or a data.frame ",
                 "with species in rows and samples in columns.")
  }

  storage.mode(x) <- "double"

  if (is.null(rownames(x))) {
    .input_error("`", arg, "` has no species names (row names).")
  }
  if (is.null(colnames(x))) {
    .input_error("`", arg, "` has no sample names (column names).")
  }
  if (anyNA(rownames(x))) {
    .input_error("Species names must not be missing.")
  }

  dup_samples <- unique(colnames(x)[duplicated(colnames(x))])
  if (length(dup_samples) > 0L) {
    .input_error("The following sample names are duplicated: ",
                 .format_sample_names(dup_samples))
  }

  x
}

# Remove species with low abundances from each sample.
.filter_low_ab <- function(ab, filtering_ab_thr_factor) {
  min_positive <- apply(ab, 2L, function(v) {
    v <- v[v > 0]
    if (length(v) == 0L) Inf else min(v)
  })
  thr <- matrix(filtering_ab_thr_factor * min_positive,
                nrow = nrow(ab), ncol = ncol(ab), byrow = TRUE)
  ab[ab <= thr] <- 0
  ab
}

.check_non_empty_samples <- function(ab) {
  empty <- colnames(ab)[colSums(ab) <= 0]
  if (length(empty) > 0L) {
    .input_error("The following samples have no non-zero species abundances: ",
                 .format_sample_names(empty))
  }
}

# Validate, optionally filter, normalize and log10-transform a table.
# Port of ab_table_utils.read_filter_normalize() after the reading step.
.prepare_ab_table <- function(x, filter_low_ab = NULL, arg = "ab_table") {
  ab <- .as_ab_matrix(x, arg)

  missing_samples <- colnames(ab)[colSums(is.na(ab)) > 0]
  if (length(missing_samples) > 0L) {
    .input_error("The following samples contain missing species abundances: ",
                 .format_sample_names(missing_samples))
  }

  negative_samples <- colnames(ab)[colSums(ab < 0) > 0]
  if (length(negative_samples) > 0L) {
    .input_error("The following samples contain negative species abundances: ",
                 .format_sample_names(negative_samples))
  }

  .check_non_empty_samples(ab)

  if (!is.null(filter_low_ab)) {
    if (!is.numeric(filter_low_ab) || length(filter_low_ab) != 1L ||
        !is.finite(filter_low_ab) || filter_low_ab <= 0) {
      .input_error("`filter_low_ab` must be a finite number greater than 0.")
    }
    ab <- .filter_low_ab(ab, filter_low_ab)
    .check_non_empty_samples(ab)
  }

  # Relative abundances, then log10 (zeros become -Inf).
  ab <- sweep(ab, 2L, colSums(ab), "/")
  log10(ab)
}

# Check that no sample name appears in both tables.
.check_distinct_sample_names <- function(ab, ab2) {
  shared <- colnames(ab)[colnames(ab) %in% colnames(ab2)]
  if (length(shared) > 0L) {
    .input_error("The following sample names are present in both abundance ",
                 "tables: ", .format_sample_names(shared), ". ",
                 "Sample names must be unique across the two tables.")
  }
}

# Warn when the two tables do not share the same species.
.compare_species_names <- function(ab, ab2) {
  s1 <- unique(rownames(ab))
  s2 <- unique(rownames(ab2))
  if (!setequal(s1, s2)) {
    warning(
      "Abundance tables have only ", length(intersect(s1, s2)),
      " species names in common. Make sure the abundance tables were ",
      "generated with the same tool and database. Missing abundance values ",
      "will be filled with zeros for non-shared species.",
      call. = FALSE
    )
  }
}

# Outer join of two log10 tables on species names (pandas join how="outer"):
# the union of species is sorted by code point and gaps are filled with -Inf.
.join_ab_tables <- function(ab, ab2) {
  if (anyDuplicated(rownames(ab)) || anyDuplicated(rownames(ab2))) {
    .input_error("Species names must be unique in each table when two ",
                 "abundance tables are compared.")
  }
  species <- union(rownames(ab), rownames(ab2))
  species <- species[.c_order(species)]
  out <- matrix(-Inf, nrow = length(species), ncol = ncol(ab) + ncol(ab2),
                dimnames = list(species, c(colnames(ab), colnames(ab2))))
  out[rownames(ab), colnames(ab)] <- ab
  out[rownames(ab2), colnames(ab2)] <- ab2
  out
}

# Prepare one or two tables for searching or plotting.
.prepare_tables <- function(ab_table, ab_table_2, filter_low_ab) {
  ab <- .prepare_ab_table(ab_table, filter_low_ab, "ab_table")
  if (is.null(ab_table_2)) {
    return(list(table = ab, samples_1 = colnames(ab), samples_2 = NULL))
  }
  ab2 <- .prepare_ab_table(ab_table_2, filter_low_ab, "ab_table_2")
  .check_distinct_sample_names(ab, ab2)
  .compare_species_names(ab, ab2)
  list(
    table = .join_ab_tables(ab, ab2),
    samples_1 = colnames(ab),
    samples_2 = colnames(ab2)
  )
}
