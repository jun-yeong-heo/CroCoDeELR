# Event-level comparison on the MetaPhlAn4 dataset (filter_low_ab = 20).
library(CroCoDeELR)
dir <- "/home/claude/data"
rd <- function(f) { d <- read_conta_events(f); d$key <- paste(d$source, d$target, sep = "\r"); d }
cmp <- function(a, b) {
  s <- intersect(a$key, b$key)
  m <- merge(a[a$key %in% s, ], b[b$key %in% s, ], by = "key", suffixes = c(".a", ".b"))
  dr <- abs(m$rate.b - m$rate.a) / m$rate.a
  c(events_a = nrow(a), events_b = nrow(b), shared = length(s),
    jaccard = length(s) / length(union(a$key, b$key)),
    max_dprob = max(abs(m$probability.b - m$probability.a)),
    median_drate = median(dr), p95_drate = unname(quantile(dr, 0.95)), max_drate = max(dr))
}
py <- rd(file.path(dir, "py_out", "metaphlan.tsv"))
out <- rbind(r_vs_py = cmp(py, rd(file.path(dir, "r_out", "metaphlan.r_rng.tsv"))),
             pyseed1_vs_py = cmp(py, rd(file.path(dir, "py_out", "metaphlan.seed1.tsv"))))
print(round(out, 4))
write.table(round(out, 4), file.path(dir, "comparison_metaphlan.tsv"), sep = "\t", quote = FALSE)
