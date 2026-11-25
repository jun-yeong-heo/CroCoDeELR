simulate_training_table <- function(ab, n_cases, seed) {
  rel <- sweep(ab, 2, colSums(ab), "/")
  set.seed(seed)
  out <- list()
  for (i in seq_len(n_cases)) {
    st <- sample(ncol(rel), 2)
    r <- 10^runif(1, -3, -1)
    out[[paste0("conta_source_case_", i)]] <- rel[, st[1]]
    out[[paste0("conta_target_case_", i)]] <- (1 - r) * rel[, st[2]] + r * rel[, st[1]]
    st <- sample(ncol(rel), 2)
    out[[paste0("non_conta_source_case_", i)]] <- rel[, st[1]]
    out[[paste0("non_conta_target_case_", i)]] <- rel[, st[2]]
  }
  m <- do.call(cbind, out)
  rownames(m) <- rownames(ab)
  m
}

test_that("a model can be trained, saved and used", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("jsonlite")
  training <- simulate_training_table(read_ab_table(toy_file()), 40, 1)
  model_file <- tempfile(fileext = ".rds")
  report_file <- tempfile(fileext = ".json")
  res <- train_model(training, ntrees = 20, model_file = model_file,
                     report_file = report_file, verbose = FALSE)
  expect_s3_class(res$model, "crocodeel_rf")
  expect_length(res$model$trees, 20L)
  expect_true(file.exists(report_file))
  rep <- jsonlite::read_json(report_file)
  expect_true(all(c("train", "test") %in% names(rep)))
  expect_true("accuracy" %in% names(rep$test$classification_report))

  # Converted trees reproduce ranger's own predictions.
  x <- matrix(runif(600, 0, 3), ncol = 6)
  x[, 1] <- sample(1:30, 100, TRUE)
  x[, 2] <- sample(8:200, 100, TRUE)
  colnames(x) <- CroCoDeELR:::.FEATURE_NAMES
  p_ranger <- predict(res$ranger, as.data.frame(x))$predictions[, "TRUE"]
  expect_equal(CroCoDeELR:::.predict_conta_proba(res$model, x), p_ranger)

  model <- load_model(model_file)
  events <- search_conta(read_ab_table(toy_file()), model = model,
                         verbose = FALSE)
  expect_s3_class(events, "data.frame")
})

test_that("training tables must follow the naming scheme", {
  ab <- read_ab_table(toy_file())
  expect_error(train_model(ab, verbose = FALSE), "do not match",
               class = "crocodeel_input_error")
  pairs <- c("conta_source_case_1", "conta_target_case_1",
             "conta_source_case_2")
  expect_error(CroCoDeELR:::.reconstruct_sample_pairs(pairs),
               "without corresponding targets", class = "crocodeel_input_error")
})

test_that("classification reports follow sklearn", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
  pred <- c(TRUE, FALSE, FALSE, TRUE, TRUE)
  r <- CroCoDeELR:::.classification_report(truth, pred)
  expect_equal(r$True$precision, 2 / 3)
  expect_equal(r$True$recall, 2 / 3)
  expect_equal(r$False$support, 2)
  expect_equal(r$accuracy, 0.6)
  expect_equal(r$`weighted avg`$recall, 0.6)
})
