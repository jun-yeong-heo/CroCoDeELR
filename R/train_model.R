# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/train_model.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.
#
# Differences from the original: the forest is trained with 'ranger' instead
# of scikit-learn, and the train/test split uses R's random number generator.
# The ranger settings mirror sklearn's defaults for RandomForestClassifier
# (bootstrap sampling, mtry = floor(sqrt(6)) = 2, fully grown trees, Gini).

.SAMPLE_NAME_PATTERN <- "^(conta_|non_conta_)(source_|target_)case_[0-9]+$"

# Reconstruct source-target pairs from the sample names of a training table.
.reconstruct_sample_pairs <- function(sample_names) {
  invalid <- sample_names[!grepl(.SAMPLE_NAME_PATTERN, sample_names)]
  if (length(invalid) > 0L) {
    .input_error("The following sample names do not match the expected ",
                 "pattern '", .SAMPLE_NAME_PATTERN, "' (",
                 .format_sample_names(invalid), ").")
  }
  sources <- sample_names[grepl("_source_", sample_names, fixed = TRUE)]
  targets <- sample_names[grepl("_target_", sample_names, fixed = TRUE)]

  without_targets <- sources[!sub("_source_", "_target_", sources,
                                  fixed = TRUE) %in% targets]
  without_sources <- targets[!sub("_target_", "_source_", targets,
                                  fixed = TRUE) %in% sources]
  msgs <- character(0)
  if (length(without_targets)) {
    msgs <- c(msgs, paste0("source samples without corresponding targets: ",
                           .format_sample_names(without_targets)))
  }
  if (length(without_sources)) {
    msgs <- c(msgs, paste0("target samples without corresponding sources: ",
                           .format_sample_names(without_sources)))
  }
  if (length(msgs)) {
    .input_error(paste(msgs, collapse = "; "))
  }

  data.frame(
    source = sources,
    target = sub("_source_", "_target_", sources, fixed = TRUE),
    stringsAsFactors = FALSE
  )
}

# Features of a batch of pairs; NA rows where extraction failed.
.features_for_pairs <- function(pairs, ab) {
  out <- matrix(NA_real_, nrow = nrow(pairs), ncol = .NUM_FEATURES)
  for (i in seq_len(nrow(pairs))) {
    f <- .extract_features(ab, pairs$source[i], pairs$target[i])
    if (!is.null(f)) {
      out[i, ] <- f$values
    }
  }
  out
}

# sklearn.metrics.classification_report(output_dict = TRUE).
.classification_report <- function(truth, pred) {
  labels <- sort(unique(c(truth, pred)))
  per_class <- lapply(labels, function(l) {
    tp <- sum(pred == l & truth == l)
    fp <- sum(pred == l & truth != l)
    fn <- sum(pred != l & truth == l)
    precision <- if (tp + fp > 0) tp / (tp + fp) else 0
    recall <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1 <- if (precision + recall > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      0
    }
    list(precision = precision, recall = recall, `f1-score` = f1,
         support = sum(truth == l))
  })
  names(per_class) <- ifelse(labels, "True", "False")
  support <- vapply(per_class, `[[`, numeric(1), "support")
  avg <- function(metric, weights) {
    v <- vapply(per_class, `[[`, numeric(1), metric)
    sum(v * weights) / sum(weights)
  }
  c(
    per_class,
    list(
      accuracy = mean(truth == pred),
      `macro avg` = list(
        precision = avg("precision", rep(1, length(labels))),
        recall = avg("recall", rep(1, length(labels))),
        `f1-score` = avg("f1-score", rep(1, length(labels))),
        support = sum(support)
      ),
      `weighted avg` = list(
        precision = avg("precision", support),
        recall = avg("recall", support),
        `f1-score` = avg("f1-score", support),
        support = sum(support)
      )
    )
  )
}

# Convert a ranger probability forest to the crocodeel_rf tree layout.
.ranger_to_crocodeel_rf <- function(fit, seed) {
  trees <- lapply(seq_len(fit$num.trees), function(k) {
    info <- ranger::treeInfo(fit, tree = k)
    info <- info[order(info$nodeID), , drop = FALSE]
    n <- nrow(info)
    if (!identical(info$nodeID, seq_len(n) - 1L)) {
      stop("Unexpected node numbering in ranger tree ", k, ".")
    }
    terminal <- info$terminal
    feature <- match(info$splitvarName, .FEATURE_NAMES)
    list(
      left = ifelse(terminal, 0L, as.integer(info$leftChild) + 1L),
      right = ifelse(terminal, 0L, as.integer(info$rightChild) + 1L),
      feature = ifelse(terminal, 0L, feature),
      threshold = ifelse(terminal, 0, info$splitval),
      prob = ifelse(terminal, info[["pred.TRUE"]], 0)
    )
  })
  .new_crocodeel_rf(
    trees,
    input_precision = "float64",
    origin = list(
      engine = paste0("ranger ", utils::packageVersion("ranger")),
      trained_with = paste0("CroCoDeELR ", utils::packageVersion("CroCoDeELR")),
      seed = seed,
      date = format(Sys.Date())
    )
  )
}

#' Train a Random Forest model for contamination detection
#'
#' Port of `crocodeel train_model`. The abundance table must contain labeled
#' sample pairs named `conta_source_case_<i>` / `conta_target_case_<i>`
#' (contaminated) and `non_conta_source_case_<i>` /
#' `non_conta_target_case_<i>` (not contaminated), as in the training dataset
#' of CroCoDeEL (<https://doi.org/10.57745/N6JSHQ>).
#'
#' Features are computed for each pair, pairs without a contamination line
#' are removed, and the remaining pairs are split into a training and a test
#' set. The forest is trained with 'ranger' using settings equivalent to
#' scikit-learn's defaults in the original tool.
#'
#' @inheritParams search_conta
#' @param test_size Fraction of sample pairs used for testing.
#' @param ntrees Number of trees.
#' @param seed Random seed for the train/test split and the forest.
#' @param model_file Optional `.rds` path where the model is saved.
#' @param report_file Optional `.json` path for the performance report
#'   (requires 'jsonlite').
#' @return A list, invisibly, with the `model` (class `crocodeel_rf`, usable
#'   in [search_conta()]), the performance `report` (sklearn-style
#'   classification reports on the train and test sets), `data` (one row per
#'   sample pair kept for training: label, `train`/`test` set and the six
#'   features) and the underlying `ranger` object.
#' @export
train_model <- function(ab_table, filter_low_ab = NULL, test_size = 0.3,
                        ntrees = 100L, seed = 0L, nproc = 1L,
                        model_file = NULL, report_file = NULL,
                        verbose = interactive()) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required to train a model.")
  }
  if (!is.null(report_file) && !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write the performance report.")
  }
  if (!is.numeric(test_size) || length(test_size) != 1L ||
      !(test_size > 0 && test_size < 1)) {
    .input_error("`test_size` must be a number strictly between 0 and 1.")
  }

  ab <- .prepare_ab_table(ab_table, filter_low_ab)
  pairs <- .reconstruct_sample_pairs(colnames(ab))
  is_contaminated <- startsWith(pairs$source, "conta")

  .log(verbose, "Abundance table contains ", nrow(pairs), " sample pairs: ",
       sum(is_contaminated), " contaminated and ", sum(!is_contaminated),
       " non-contaminated")

  start <- proc.time()[["elapsed"]]
  batches <- lapply(.chunk(seq_len(nrow(pairs)), 50L),
                    function(idx) pairs[idx, , drop = FALSE])
  features <- do.call(rbind, .par_lapply(batches, .features_for_pairs,
                                         nproc = nproc,
                                         progress = verbose && nproc <= 1,
                                         ab = ab))
  .log(verbose, "Feature computation completed in ",
       sprintf("%.1f", proc.time()[["elapsed"]] - start), " seconds")

  invalid <- rowSums(is.na(features)) > 0
  if (sum(is_contaminated & invalid)) {
    .log(verbose, "Feature computation failed for ",
         sum(is_contaminated & invalid), " contaminated sample pairs")
  }
  if (sum(!is_contaminated & invalid)) {
    .log(verbose, "Feature computation failed for ",
         sum(!is_contaminated & invalid), " non-contaminated sample pairs")
  }
  features <- features[!invalid, , drop = FALSE]
  is_contaminated <- is_contaminated[!invalid]
  valid_pairs <- pairs[!invalid, , drop = FALSE]
  colnames(features) <- .FEATURE_NAMES

  n <- nrow(features)
  n_test <- ceiling(test_size * n)
  if (n_test < 1L || n - n_test < 1L) {
    .input_error("Not enough valid sample pairs to split into training and ",
                 "test sets.")
  }
  .log(verbose, n, " sample pairs remain: ", sum(is_contaminated),
       " contaminated, ", sum(!is_contaminated), " non-contaminated")

  perm <- .with_seed(seed, sample.int(n))
  test_idx <- perm[seq_len(n_test)]
  train_idx <- perm[-seq_len(n_test)]

  x_train <- as.data.frame(features[train_idx, , drop = FALSE])
  y_train <- factor(is_contaminated[train_idx], levels = c(FALSE, TRUE))

  .log(verbose, "Training Random Forest model with ", ntrees, " trees...")
  start <- proc.time()[["elapsed"]]
  fit <- ranger::ranger(
    x = x_train,
    y = y_train,
    num.trees = ntrees,
    mtry = floor(sqrt(.NUM_FEATURES)),
    min.node.size = 1L,
    replace = TRUE,
    sample.fraction = 1,
    splitrule = "gini",
    probability = TRUE,
    num.threads = as.integer(nproc),
    seed = seed,
    verbose = FALSE
  )
  .log(verbose, "Training completed in ",
       sprintf("%.1f", proc.time()[["elapsed"]] - start), " seconds")

  model <- .ranger_to_crocodeel_rf(fit, seed)

  predict_class <- function(idx) {
    .predict_conta_proba(model, features[idx, , drop = FALSE]) > 0.5
  }
  report <- list(
    train = list(classification_report = .classification_report(
      is_contaminated[train_idx], predict_class(train_idx))),
    test = list(classification_report = .classification_report(
      is_contaminated[test_idx], predict_class(test_idx)))
  )

  if (!is.null(model_file)) {
    saveRDS(model, model_file)
    .log(verbose, "Model saved to ", model_file)
  }
  if (!is.null(report_file)) {
    jsonlite::write_json(report, report_file, auto_unbox = TRUE, digits = NA,
                         pretty = 4)
    .log(verbose, "Model performance report saved to ", report_file)
  }

  set <- rep("train", n)
  set[test_idx] <- "test"
  data <- data.frame(valid_pairs, is_contaminated = is_contaminated,
                     set = set, features, stringsAsFactors = FALSE)
  rownames(data) <- NULL

  invisible(list(model = model, report = report, data = data, ranger = fit))
}
