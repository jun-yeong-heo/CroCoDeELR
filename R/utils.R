# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

# Signal an input error with a dedicated condition class.
.input_error <- function(...) {
  msg <- paste0(...)
  cond <- structure(
    class = c("crocodeel_input_error", "error", "condition"),
    list(message = msg, call = NULL)
  )
  stop(cond)
}

# Format sample names, truncating long lists (utils.format_sample_names).
.format_sample_names <- function(sample_names, max_names = 5L) {
  if (length(sample_names) <= max_names) {
    return(paste(sample_names, collapse = ", "))
  }
  paste0(
    length(sample_names), " sample names (showing ", max_names, "): ",
    paste(sample_names[seq_len(max_names)], collapse = ", ")
  )
}

# Order strings by code point, like Python's sort, independent of locale.
.c_order <- function(...) {
  order(..., method = "radix")
}

# Format a double like Python's str()/repr(): shortest round-trip digits.
.py_float_str <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) {
      return("nan")
    }
    if (is.infinite(v)) {
      return(if (v > 0) "inf" else "-inf")
    }
    if (v == 0) {
      return(if (1 / v < 0) "-0.0" else "0.0")
    }
    # Smallest number of significant digits that round-trips.
    d <- 1L
    while (d < 17L && as.numeric(sprintf("%.*e", d - 1L, v)) != v) {
      d <- d + 1L
    }
    sci <- sprintf("%.*e", d - 1L, v)
    exponent <- as.integer(sub("^.*e", "", sci))
    # Python uses positional notation for 1e-4 <= |x| < 1e16.
    if (exponent >= -4L && exponent < 16L) {
      s <- sprintf("%.*f", max(0L, d - 1L - exponent), v)
      if (!grepl(".", s, fixed = TRUE)) {
        s <- paste0(s, ".0")
      }
      return(s)
    }
    mantissa <- sub("e.*$", "", sci)
    paste0(mantissa, "e", if (exponent < 0) "-" else "+",
           formatC(abs(exponent), width = 2, flag = "0"))
  }, character(1), USE.NAMES = FALSE)
}

# Normalize C exponents to at least two digits ("1e-005" -> "1e-05").
.normalize_exponent <- function(s) {
  sub("e([+-])0*([0-9]{2,})$", "e\\1\\2", sub("e([+-])([0-9])$", "e\\10\\2", s))
}

# Format with two decimals in exponent notation, as Python's f"{x:.2e}".
.fmt_rate <- function(x) {
  .normalize_exponent(sprintf("%.2e", x))
}

# numpy.round(): round half to even after scaling by a power of ten.
.np_round <- function(x, decimals) {
  if (decimals >= 0) {
    scale <- 10^decimals
    round(x * scale) / scale
  } else {
    scale <- 10^(-decimals)
    round(x / scale) * scale
  }
}

# Run code with a fixed seed and RNG kind, then restore the caller's RNG.
.with_seed <- function(seed, code) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  }
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(seed)
  force(code)
}

# Log an informational message unless quiet.
.log <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " :: INFO :: ", ...)
  }
}

# Split a vector into consecutive chunks of at most `size` elements.
.chunk <- function(x, size) {
  if (length(x) == 0L) {
    return(list())
  }
  split(x, ceiling(seq_along(x) / size))
}

# Apply FUN over a list, optionally on a PSOCK cluster (works on Windows).
.par_lapply <- function(X, FUN, nproc, progress = FALSE, ...) {
  nproc <- as.integer(nproc)
  if (length(X) == 0L) {
    return(list())
  }
  if (nproc <= 1L || length(X) == 1L) {
    if (progress) {
      pb <- utils::txtProgressBar(min = 0, max = length(X), style = 3)
      on.exit(close(pb), add = TRUE)
    }
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
      out[[i]] <- FUN(X[[i]], ...)
      if (progress) utils::setTxtProgressBar(pb, i)
    }
    return(out)
  }
  cl <- parallel::makePSOCKcluster(min(nproc, length(X)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function() {
    suppressPackageStartupMessages(requireNamespace("CroCoDeELR", quietly = TRUE))
    NULL
  })
  parallel::parLapplyLB(cl, X, FUN, ...)
}
