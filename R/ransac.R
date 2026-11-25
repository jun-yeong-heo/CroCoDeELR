# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/ransac.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.
#
# Difference from the original: candidate subsets are drawn with R's
# Mersenne-Twister generator (sample.int, seed 42) instead of NumPy's
# generator through sklearn. Model selection, tie-breaking and the dynamic
# stopping rule are the same as in the original code.

# Cache of candidate subsets, keyed by number of observations and settings.
# Subsets do not depend on the data, only on these settings.
.ransac_cache <- new.env(parent = emptyenv())

.ransac_key <- function(num_samples, max_trials, min_samples, random_state) {
  paste(num_samples, max_trials, min_samples, random_state, sep = "_")
}

# Return a max_trials x min_samples matrix of 1-based observation indices.
.candidate_subsets <- function(num_samples, max_trials, min_samples,
                               random_state) {
  key <- .ransac_key(num_samples, max_trials, min_samples, random_state)
  subsets <- .ransac_cache[[key]]
  if (is.null(subsets)) {
    subsets <- .with_seed(random_state, {
      matrix(
        unlist(lapply(seq_len(max_trials), function(i) {
          sample.int(num_samples, min_samples)
        })),
        nrow = max_trials, ncol = min_samples, byrow = TRUE
      )
    })
    assign(key, subsets, envir = .ransac_cache)
  }
  subsets
}

# Number of trials required by RANSAC's stopping criterion (sklearn).
.dynamic_max_trials <- function(num_inliers, num_samples, min_samples,
                                probability) {
  eps <- .Machine$double.eps
  inlier_ratio <- num_inliers / num_samples
  numerator <- max(eps, 1 - probability)
  denominator <- max(eps, 1 - inlier_ratio^min_samples)
  if (numerator == 1) {
    return(0)
  }
  if (denominator == 1) {
    return(Inf)
  }
  abs(ceiling(log(numerator) / log(denominator)))
}

# Fit ordinate = abscissa + intercept with RANSAC.
# Returns NULL when no valid consensus set is found, otherwise a list with
# the inlier mask, the fitted intercept and the number of trials performed.
.fit_unit_slope_ransac <- function(abscissa, ordinate, max_trials,
                                   residual_threshold, random_state,
                                   min_samples = 2L, stop_probability = 0.99) {
  residuals <- ordinate - abscissa
  num_samples <- length(residuals)

  subsets <- .candidate_subsets(num_samples, max_trials, min_samples,
                                random_state)

  # Intercept of each candidate = mean residual of its sampled observations.
  candidate_offsets <- rowSums(
    matrix(residuals[subsets], nrow = max_trials)
  ) / min_samples

  candidate_inliers <- abs(
    matrix(residuals, nrow = max_trials, ncol = num_samples, byrow = TRUE) -
      candidate_offsets
  ) <= residual_threshold
  candidate_num_inliers <- rowSums(candidate_inliers)

  best_trial <- -1L
  best_num_inliers <- 1
  best_score <- -Inf
  dynamic_max_trials <- max_trials
  num_trials <- 0L

  while (num_trials < dynamic_max_trials) {
    num_trials <- num_trials + 1L
    trial <- num_trials

    if (candidate_num_inliers[trial] < best_num_inliers) {
      next
    }

    inlier_residuals <- residuals[candidate_inliers[trial, ]]
    deviations <- (inlier_residuals - candidate_offsets[trial])^2
    score <- -(sum(deviations) / length(deviations))

    if (candidate_num_inliers[trial] == best_num_inliers &&
        score < best_score) {
      next
    }

    best_trial <- trial
    best_num_inliers <- candidate_num_inliers[trial]
    best_score <- score

    dynamic_max_trials <- min(
      dynamic_max_trials,
      .dynamic_max_trials(best_num_inliers, num_samples, min_samples,
                          stop_probability)
    )
  }

  if (best_trial < 0L) {
    return(NULL)
  }

  inlier_mask <- candidate_inliers[best_trial, ]
  inlier_residuals <- residuals[inlier_mask]

  list(
    inlier_mask = inlier_mask,
    intercept = sum(inlier_residuals) / length(inlier_residuals),
    num_trials = num_trials
  )
}
