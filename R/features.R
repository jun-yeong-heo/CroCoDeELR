# CroCoDeELR: unofficial R port of CroCoDeEL 1.3.1 (modified: translated
# from Python to R by Jun-yeong Heo, November 2026)
# Port of crocodeel/conta_features.py
# Copyright (C) the CroCoDeEL authors (original Python code)
# Copyright (C) Jun-yeong Heo (R port)
# License: GPL (>= 3). See the COPYING file.

.FEATURE_NAMES <- c(
  "ransac_num_trials",
  "num_species_conta_line",
  "mean_distance_to_conta_line",
  "mean_distance_to_nearest_neighbors",
  "diff_mean_ab_top10_source_species_vs_ab_cutoff1",
  "diff_mean_ab_top10_source_species_vs_ab_cutoff2"
)
.NUM_FEATURES <- 6L
.CONTA_LINE_MIN_NUM_SPECIES <- 8L
.UPPER_LEFT_QUADRANT_MAX_NUM_SPECIES <- 2L
.RANSAC_MAX_TRIALS <- 30L
.RANSAC_RANDOM_STATE <- 42L
.RANSAC_RESIDUAL_THRESHOLD <- 0.2

# numpy.quantile(x, q) with the default "linear" method, including NumPy's
# two-sided interpolation formula.
.np_quantile <- function(x, q) {
  x <- sort(x)
  n <- length(x)
  vapply(q, function(p) {
    virtual_index <- (n - 1) * p
    previous_index <- floor(virtual_index)
    next_index <- min(previous_index + 1, n - 1)
    gamma <- virtual_index - previous_index
    a <- x[previous_index + 1]
    b <- x[next_index + 1]
    diff_b_a <- b - a
    if (gamma >= 0.5) b - diff_b_a * (1 - gamma) else a + diff_b_a * gamma
  }, numeric(1))
}

# Mean distance of each point to its k nearest neighbours (the point itself
# included, as in sklearn's kneighbors() on the training data).
.mean_distance_to_nearest_neighbors <- function(x, y, num_neighbors = 5L) {
  n <- length(x)
  k <- min(num_neighbors, n)
  dx <- outer(x, x, "-")
  dy <- outer(y, y, "-")
  d <- sqrt(dx * dx + dy * dy)
  nearest <- apply(d, 1L, function(row) sort.int(row, partial = seq_len(k))[seq_len(k)])
  mean(nearest)
}

# Select candidate species for a contamination line.
# Returns a logical mask over the species of the upper triangle.
.conta_line_candidates <- function(x, y) {
  is_left <- outer(x, x, "<=")
  is_upper <- outer(y, y, ">=")
  upper_left_count <- colSums(is_left & is_upper) - 1
  upper_left_count <- as.integer(upper_left_count)
  upper_left_count <= .UPPER_LEFT_QUADRANT_MAX_NUM_SPECIES
}

# Extract the contamination features for one source-target pair.
# `ab` is a log10 relative abundance matrix (zeros as -Inf).
# Returns NULL when no contamination line can be identified.
.extract_features <- function(ab, source, target) {
  target_ab <- ab[, target]
  source_ab <- ab[, source]

  # Step 1: species present in the target and at least as abundant in the
  # source, then those with few species in their upper-left quadrant.
  upper_triangle <- (source_ab >= target_ab) & (target_ab != -Inf)
  x_upper <- target_ab[upper_triangle]
  y_upper <- source_ab[upper_triangle]

  if (length(x_upper) < .CONTA_LINE_MIN_NUM_SPECIES) {
    return(NULL)
  }

  candidates <- .conta_line_candidates(x_upper, y_upper)
  cand_x <- x_upper[candidates]
  cand_y <- y_upper[candidates]

  if (length(cand_x) < .CONTA_LINE_MIN_NUM_SPECIES) {
    return(NULL)
  }

  # Step 2: contamination line offset with RANSAC (slope fixed to one).
  fit <- .fit_unit_slope_ransac(
    cand_x, cand_y,
    max_trials = .RANSAC_MAX_TRIALS,
    residual_threshold = .RANSAC_RESIDUAL_THRESHOLD,
    random_state = .RANSAC_RANDOM_STATE
  )
  if (is.null(fit)) {
    return(NULL)
  }

  line_x <- cand_x[fit$inlier_mask]
  line_y <- cand_y[fit$inlier_mask]
  offset <- fit$intercept

  if (length(line_x) < .CONTA_LINE_MIN_NUM_SPECIES) {
    return(NULL)
  }

  # Step 3: features describing the contamination line.
  source_specific <- (target_ab == -Inf) & (source_ab != -Inf)

  num_species_line <- length(line_x)
  mean_nn <- .mean_distance_to_nearest_neighbors(line_x, line_y)
  distances <- abs(line_y - line_x - offset) / sqrt(2)
  mean_distance_line <- sum(distances) / length(distances)

  if (!any(source_specific)) {
    diff_cutoff1 <- 0
    diff_cutoff2 <- 0
  } else {
    # 'm' in the paper: mean source abundance of the 10 most abundant
    # source-specific species.
    values <- source_ab[source_specific]
    if (length(values) > 10L) {
      values <- sort.int(values, decreasing = TRUE)[1:10]
    }
    m <- sum(values) / length(values)

    # 'c1': pseudo zero shifted by the contamination line offset.
    pseudo_zero <- min(target_ab[target_ab != -Inf]) - 1
    diff_cutoff1 <- abs(m - (pseudo_zero + offset))

    # 'c2': mean source abundance of the 10% least abundant line species.
    num_low <- max(as.integer(0.1 * length(line_y)), 1L)
    low_values <- sort.int(line_y)[seq_len(num_low)]
    diff_cutoff2 <- abs(m - sum(low_values) / length(low_values))
  }

  values <- c(
    fit$num_trials,
    num_species_line,
    mean_distance_line,
    mean_nn,
    diff_cutoff1,
    diff_cutoff2
  )

  list(
    values = values,
    conta_line_offset = offset,
    conta_line_species = .refine_conta_line(
      rownames(ab), target_ab, source_ab, upper_triangle, line_x, line_y
    )
  )
}

# Species belonging to the refined contamination line.
.refine_conta_line <- function(species, target_ab, source_ab, upper_triangle,
                               line_x, line_y) {
  ratio <- line_y - line_x
  q <- .np_quantile(ratio, c(0.25, 0.75))
  iqr <- q[2] - q[1]
  lower <- max(q[1] - 1.5 * iqr, min(ratio))
  upper <- min(q[2] + 1.5 * iqr, max(ratio))

  upper_ratio <- source_ab[upper_triangle] - target_ab[upper_triangle]
  keep <- upper_ratio >= lower & upper_ratio <= upper
  species[upper_triangle][keep]
}
