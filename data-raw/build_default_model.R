# Build inst/extdata/crocodeel_rf_Feb2026_2.rds from the CroCoDeEL 1.3.1
# default model, without retraining.
#
# 1. In a Python environment with CroCoDeEL 1.3.1 (scikit-learn >= 1.8):
#      python inst/python/export_sklearn_model.py \
#        <crocodeel>/models/crocodeel_rf_Feb2026_2.joblib default_model.json
# 2. In R, from the package root:
devtools::load_all()
m <- import_sklearn_model("default_model.json")
m$origin <- list(
  engine = "scikit-learn RandomForestClassifier",
  sklearn_version = "1.8.0",
  upstream = "CroCoDeEL 1.3.1 (crocodeel/models/crocodeel_rf_Feb2026_2.joblib)",
  converted_with = "inst/python/export_sklearn_model.py (no retraining)"
)
saveRDS(m, "inst/extdata/crocodeel_rf_Feb2026_2.rds", compress = "xz")
