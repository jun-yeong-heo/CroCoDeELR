# Run CroCoDeELR on the validation datasets.
# mode "r_rng": package as installed (R random numbers), all pairs (cutoff 0).
# mode "numpy": RANSAC subsets of the Python version injected, single process.
args <- commandArgs(TRUE)
mode <- args[1]; datasets <- args[-1]
library(CroCoDeELR)
source("/home/claude/CroCoDeELR/validation/helpers.R")
if (mode == "numpy") {
  inject_numpy_subsets(NULL, "/home/claude/data/numpy_subsets.tsv", CroCoDeELR:::.ransac_cache)
}
dir.create("/home/claude/data/r_out", showWarnings = FALSE)
for (d in datasets) {
  file <- if (d == "toy") system.file("extdata", "mgs_profiles_test.tsv", package = "CroCoDeELR") else
    file.path("/home/claude/data", paste0(d, ".meteor.tab"))
  ab <- read_ab_table(file)
  t0 <- proc.time()[["elapsed"]]
  ev <- search_conta(ab, probability_cutoff = 0, nproc = if (mode == "numpy") 1 else 2, verbose = FALSE)
  el <- proc.time()[["elapsed"]] - t0
  write_conta_events(ev, file.path("/home/claude/data/r_out", paste0(d, ".", mode, ".all.tsv")))
  cat(d, mode, ncol(ab), "samples", sprintf("%.1f", el), "s\n", file = "/home/claude/data/r_out/times.txt", append = TRUE)
}
