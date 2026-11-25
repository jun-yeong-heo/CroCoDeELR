# train_model() on the official training table; keep per-pair results.
library(CroCoDeELR)
dir <- "/home/claude/data/train"
ab <- read_ab_table(file.path(dir, "training_dataset.meteor.tsv"))
t0 <- proc.time()[["elapsed"]]
res <- train_model(ab, nproc = 2, model_file = file.path(dir, "r_model.rds"),
                   report_file = file.path(dir, "r_report.json"), verbose = TRUE)
cat("train_model elapsed:", proc.time()[["elapsed"]] - t0, "s\n")
d <- res$data
d$prob_retrained <- CroCoDeELR:::.predict_conta_proba(res$model, as.matrix(d[, CroCoDeELR:::.FEATURE_NAMES]))
d$prob_default <- CroCoDeELR:::.predict_conta_proba(load_model(), as.matrix(d[, CroCoDeELR:::.FEATURE_NAMES]))
write.table(d, file.path(dir, "r_train_eval.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
