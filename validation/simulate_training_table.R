# Build a labeled training table by in-silico contamination of real profiles.
# Used only to check that train_model() works like `crocodeel train_model`;
# the official training dataset (doi:10.57745/N6JSHQ) could not be downloaded
# in the validation environment.
# Usage: Rscript simulate_training_table.R input.tab output.tsv n_cases seed log10_rate_min log10_rate_max
args <- commandArgs(TRUE)
ab <- CroCoDeELR::read_ab_table(args[1])
n_cases <- as.integer(args[3]); set.seed(as.integer(args[4]))
rel <- sweep(ab, 2, colSums(ab), "/")
out <- list()
for (i in seq_len(n_cases)) {
  st <- sample(ncol(rel), 2)
  s <- rel[, st[1]]; t <- rel[, st[2]]
  r <- 10^runif(1, as.numeric(args[5]), as.numeric(args[6]))
  out[[paste0("conta_source_case_", i)]] <- s
  out[[paste0("conta_target_case_", i)]] <- (1 - r) * t + r * s
  st <- sample(ncol(rel), 2)
  out[[paste0("non_conta_source_case_", i)]] <- rel[, st[1]]
  out[[paste0("non_conta_target_case_", i)]] <- rel[, st[2]]
}
m <- do.call(cbind, out); rownames(m) <- rownames(ab)
write.table(cbind(id = rownames(m), as.data.frame(m)), args[2], sep = "\t",
            quote = FALSE, row.names = FALSE)
