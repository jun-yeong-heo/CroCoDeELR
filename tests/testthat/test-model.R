test_that("the default model gives the same probabilities as scikit-learn", {
  ref <- utils::read.delim(test_path("fixtures", "sklearn_predictions.tsv"))
  model <- load_model()
  expect_s3_class(model, "crocodeel_rf")
  expect_length(model$trees, 100L)
  p <- CroCoDeELR:::.predict_conta_proba(model, as.matrix(ref[, 1:6]))
  # Includes rows whose features sit exactly on split thresholds.
  expect_identical(p, ref$proba)
})

test_that("invalid models are rejected", {
  expect_error(search_conta(read_ab_table(toy_file()), model = list()),
               class = "crocodeel_input_error")
  f <- tempfile(fileext = ".rds")
  saveRDS(list(a = 1), f)
  expect_error(load_model(f), class = "crocodeel_input_error")
})

test_that("an sklearn JSON export can be imported", {
  skip_if_not_installed("jsonlite")
  json <- tempfile(fileext = ".json")
  tree <- list(children_left = c(1, -1, -1), children_right = c(2, -1, -1),
               feature = c(0, -2, -2), threshold = c(5, -2, -2),
               value = list(c(0.5, 0.5), c(1, 0), c(0.25, 0.75)))
  jsonlite::write_json(list(format = "crocodeel_rf_sklearn_v1",
                            classes = c(FALSE, TRUE), n_features = 6,
                            trees = list(tree, tree)),
                       json, auto_unbox = TRUE, digits = NA)
  m <- import_sklearn_model(json)
  x <- rbind(c(5, 0, 0, 0, 0, 0), c(6, 0, 0, 0, 0, 0))
  expect_equal(CroCoDeELR:::.predict_conta_proba(m, x), c(0, 0.75))
})
