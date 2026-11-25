# Apply models trained on the official CroCoDeEL training dataset by
# `crocodeel train_model` (scikit-learn, imported exactly) and by
# train_model() (ranger) to real datasets, with identical features
# (R RANSAC draws). The shipped default model is shown for reference.
library(CroCoDeELR)
dir <- "/home/claude/data"
models <- list(
  sklearn_retrained = import_sklearn_model(file.path(dir, "train", "py_model.json")),
  ranger_retrained = load_model(file.path(dir, "train", "r_model.rds")),
  default = load_model()
)
datasets <- c("PRJEB83730", "PRJNA698986_P2", "PRJNA698986_P3", "PRJEB12449",
              "PRJEB10878", "PRJNA352475", "PRJEB6337")
res <- list()
for (d in datasets) {
  ab <- CroCoDeELR:::.prepare_ab_table(read_ab_table(file.path(dir, paste0(d, ".meteor.tab"))))
  f <- list()
  for (s in colnames(ab)) for (t in colnames(ab)) {
    if (s == t) next
    x <- CroCoDeELR:::.extract_features(ab, s, t)
    if (!is.null(x)) f[[length(f) + 1]] <- x$values
  }
  f <- do.call(rbind, f)
  p <- lapply(models, function(m) CroCoDeELR:::.predict_conta_proba(m, f))
  ev <- lapply(p, function(x) x >= 0.5)
  jac <- function(a, b) if (any(a | b)) sum(a & b) / sum(a | b) else 1
  res[[d]] <- data.frame(
    dataset = d, pairs_with_line = nrow(f),
    events_sklearn = sum(ev$sklearn_retrained), events_ranger = sum(ev$ranger_retrained),
    events_default = sum(ev$default),
    jaccard_ranger_vs_sklearn = jac(ev$ranger_retrained, ev$sklearn_retrained),
    jaccard_sklearn_vs_default = jac(ev$sklearn_retrained, ev$default),
    jaccard_ranger_vs_default = jac(ev$ranger_retrained, ev$default),
    cor_prob_ranger_sklearn = cor(p$ranger_retrained, p$sklearn_retrained)
  )
}
out <- do.call(rbind, res); rownames(out) <- NULL
print(out, digits = 3)
write.table(out, file.path(dir, "comparison_trained_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
