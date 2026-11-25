# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Random Forest model storage and prediction.
# Copyright (C) the CroCoDeEL authors (original Python code and default model)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.
#
# A model is stored as a list of decision trees in the same layout as
# sklearn's tree arrays. Prediction reproduces sklearn's
# RandomForestClassifier.predict_proba(): features are cast to single
# precision (for models imported from sklearn), each tree returns the class
# fractions of the reached leaf, and the fractions are averaged over trees.

.DEFAULT_MODEL_FILE <- "crocodeel_rf_Feb2026_2.rds"

# Build a crocodeel_rf object from per-tree arrays.
.new_crocodeel_rf <- function(trees, input_precision, origin) {
  structure(
    list(
      trees = trees,
      n_features = .NUM_FEATURES,
      feature_names = .FEATURE_NAMES,
      input_precision = input_precision,
      origin = origin
    ),
    class = "crocodeel_rf"
  )
}

#' @export
print.crocodeel_rf <- function(x, ...) {
  cat("CroCoDeEL Random Forest model\n")
  cat("  trees:           ", length(x$trees), "\n")
  cat("  features:        ", x$n_features, "\n")
  cat("  input precision: ", x$input_precision, "\n")
  for (nm in names(x$origin)) {
    cat("  ", nm, ": ", paste(format(x$origin[[nm]]), collapse = ", "), "\n",
        sep = "")
  }
  invisible(x)
}

.validate_model <- function(model, what = "model") {
  if (!inherits(model, "crocodeel_rf")) {
    .input_error("The ", what, " is not a CroCoDeELR Random Forest model ",
                 "(class 'crocodeel_rf').")
  }
  if (!identical(as.integer(model$n_features), .NUM_FEATURES)) {
    .input_error("The Random Forest ", what, " expects ", model$n_features,
                 " features, but ", .NUM_FEATURES, " are required.")
  }
  model
}

#' Load a CroCoDeELR Random Forest model
#'
#' With `file = NULL`, loads the default model of CroCoDeEL
#' (`crocodeel_rf_Feb2026_2`), converted from the original scikit-learn file
#' without retraining. Otherwise reads a model saved by [train_model()] or by
#' [import_sklearn_model()].
#'
#' @param file Path to an `.rds` model file, or `NULL` for the default model.
#' @return An object of class `crocodeel_rf`.
#' @examples
#' load_model()
#' @export
load_model <- function(file = NULL) {
  if (is.null(file)) {
    file <- system.file("extdata", .DEFAULT_MODEL_FILE, package = "CroCoDeELR")
    if (!nzchar(file)) {
      stop("The default model file is missing from the CroCoDeELR installation.")
    }
  } else if (!file.exists(file)) {
    .input_error("The model file '", file, "' does not exist.")
  }
  model <- tryCatch(
    readRDS(file),
    error = function(e) {
      .input_error("The model file '", file, "' could not be read: ",
                   conditionMessage(e))
    }
  )
  .validate_model(model, paste0("model file ", file))
}

#' Import a Random Forest model trained with the Python CroCoDeEL
#'
#' Converts a model trained by `crocodeel train_model` (scikit-learn) into a
#' `crocodeel_rf` object. First export the `.joblib` file to JSON with the
#' Python script shipped in this package:
#'
#' ```
#' python <path>/export_sklearn_model.py model.joblib model.json
#' ```
#'
#' where `<path>` is `system.file("python", package = "CroCoDeELR")`.
#' Predictions of the imported model are the same as those of the
#' scikit-learn model.
#'
#' @param json_file Path to the JSON file written by `export_sklearn_model.py`.
#' @param output Optional path of an `.rds` file where the model is saved.
#' @return An object of class `crocodeel_rf` (invisibly when `output` is set).
#' @export
import_sklearn_model <- function(json_file, output = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to import scikit-learn models.")
  }
  payload <- jsonlite::fromJSON(json_file, simplifyVector = TRUE)
  if (!identical(payload$format, "crocodeel_rf_sklearn_v1")) {
    .input_error("'", json_file, "' was not written by export_sklearn_model.py.")
  }
  if (!identical(as.logical(payload$classes), c(FALSE, TRUE))) {
    .input_error("The model must be a binary classifier with classes ",
                 "False and True.")
  }
  if (!identical(as.integer(payload$n_features), .NUM_FEATURES)) {
    .input_error("The Random Forest model expects ", payload$n_features,
                 " features, but ", .NUM_FEATURES, " are required.")
  }

  trees_df <- payload$trees
  trees <- lapply(seq_len(nrow(trees_df)), function(i) {
    left <- as.integer(trees_df$children_left[[i]])
    right <- as.integer(trees_df$children_right[[i]])
    value <- trees_df$value[[i]]
    if (!is.matrix(value)) {
      value <- matrix(value, ncol = 2L)
    }
    # sklearn's DecisionTreeClassifier.predict_proba() normalization.
    normalizer <- value[, 1] + value[, 2]
    normalizer[normalizer == 0] <- 1
    list(
      # 1-based children; 0 marks a leaf (sklearn uses -1).
      left = ifelse(left < 0L, 0L, left + 1L),
      right = ifelse(right < 0L, 0L, right + 1L),
      feature = as.integer(trees_df$feature[[i]]) + 1L,
      threshold = as.numeric(trees_df$threshold[[i]]),
      prob = value[, 2] / normalizer
    )
  })

  model <- .new_crocodeel_rf(
    trees,
    input_precision = "float32",
    origin = list(
      engine = "scikit-learn RandomForestClassifier",
      source_file = basename(json_file)
    )
  )

  if (!is.null(output)) {
    saveRDS(model, output)
    return(invisible(model))
  }
  model
}

# Round doubles to single precision (as numpy's astype(float32)).
.to_float32 <- function(x) {
  out <- readBin(writeBin(as.double(x), raw(), size = 4L),
                 what = "double", size = 4L, n = length(x))
  dim(out) <- dim(x)
  out
}

# Probability of contamination for each row of a feature matrix.
.predict_conta_proba <- function(model, features) {
  features <- matrix(as.double(features), ncol = .NUM_FEATURES)
  n <- nrow(features)
  if (n == 0L) {
    return(numeric(0))
  }
  if (identical(model$input_precision, "float32")) {
    features <- .to_float32(features)
  }

  rows <- seq_len(n)
  proba <- numeric(n)

  for (tree in model$trees) {
    node <- rep.int(1L, n)
    active <- tree$left[node] > 0L
    while (any(active)) {
      idx <- rows[active]
      current <- node[idx]
      go_left <- features[cbind(idx, tree$feature[current])] <=
        tree$threshold[current]
      node[idx] <- ifelse(go_left, tree$left[current], tree$right[current])
      active <- tree$left[node] > 0L
    }
    proba <- proba + tree$prob[node]
  }

  proba / length(model$trees)
}
