# Compare train_model() (R, ranger) with `crocodeel train_model` (Python,
# scikit-learn) on the official CroCoDeEL training dataset.
library(CroCoDeELR)
dir <- "/home/claude/data/train"
py <- read.delim(file.path(dir, "py_train_eval.tsv"))
r <- read.delim(file.path(dir, "r_train_eval.tsv"))
py$is_contaminated <- py$is_contaminated == "True"

metrics <- function(truth, prob) {
  pred <- prob > 0.5
  tp <- sum(pred & truth); fp <- sum(pred & !truth); fn <- sum(!pred & truth)
  pos <- prob[truth]; neg <- prob[!truth]
  auc <- (sum(rank(c(pos, neg))[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
  c(n = length(truth), n_conta = sum(truth), accuracy = mean(pred == truth),
    precision = tp / (tp + fp), recall = tp / (tp + fn),
    f1 = 2 * tp / (2 * tp + fp + fn), auc = auc)
}

cat("Python retrained model vs shipped default model, max |prob diff| on Python features:",
    max(abs(py$prob_retrained - py$prob_default)), "\n")

# 1. Each tool on its own test split.
own <- rbind(
  python_own_test = metrics(py$is_contaminated[py$set == "test"], py$prob_retrained[py$set == "test"]),
  r_own_test = metrics(r$is_contaminated[r$set == "test"], r$prob_retrained[r$set == "test"])
)

# 2. Pairs held out by both tools.
common <- intersect(py$source[py$set == "test"], r$source[r$set == "test"])
pyc <- py[match(common, py$source), ]
rc <- r[match(common, r$source), ]
stopifnot(identical(pyc$is_contaminated, rc$is_contaminated))
shared_test <- rbind(
  python_pipeline = metrics(pyc$is_contaminated, pyc$prob_retrained),
  r_pipeline = metrics(rc$is_contaminated, rc$prob_retrained)
)
same_call <- mean((pyc$prob_retrained > 0.5) == (rc$prob_retrained > 0.5))

out <- rbind(own, shared_test)
print(round(out, 4))
cat("Pairs held out by both:", length(common), " same call:", round(same_call, 4), "\n")
cat("Pairs kept: Python", nrow(py), " R", nrow(r), " shared",
    length(intersect(py$source, r$source)), "\n")
write.table(round(out, 4), file.path(dir, "comparison_training.tsv"), sep = "\t", quote = FALSE)

# 3. Rate estimation on the contaminated training pairs.
meta <- read.delim(file.path(dir, "training_dataset.metadata.tsv"))
meta <- meta[startsWith(meta$source_name, "conta"), ]
ab <- CroCoDeELR:::.prepare_ab_table(read_ab_table(file.path(dir, "training_dataset.meteor.tsv")))
r_rate <- vapply(seq_len(nrow(meta)), function(i) {
  f <- CroCoDeELR:::.extract_features(ab, meta$source_name[i], meta$target_name[i])
  if (is.null(f)) NA_real_ else 10^(-f$conta_line_offset)
}, numeric(1))
ok <- !is.na(r_rate) & !is.na(meta$crocodeel_conta_rate) & meta$crocodeel_conta_rate > 0
lerr <- function(a, b) abs(log10(a) - log10(b))
rate_tab <- rbind(
  r_vs_theoretical = summary(lerr(r_rate[ok], meta$theoretical_conta_rate[ok])),
  python_vs_theoretical = summary(lerr(meta$crocodeel_conta_rate[ok], meta$theoretical_conta_rate[ok])),
  r_vs_python = summary(lerr(r_rate[ok], meta$crocodeel_conta_rate[ok]))
)
cat("\nAbsolute log10 rate error, contaminated pairs (n =", sum(ok), ")\n")
print(round(rate_tab, 4))
write.table(round(rate_tab, 4), file.path(dir, "comparison_training_rates.tsv"), sep = "\t", quote = FALSE)
