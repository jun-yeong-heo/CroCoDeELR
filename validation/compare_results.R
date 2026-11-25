# Compare CroCoDeELR results with those of the Python CroCoDeEL 1.3.1.
#
# Inputs are event files written with probability cutoff 0, i.e. every sample
# pair with a contamination line:
#   py_out/<d>.all.tsv          Python, RANSAC seed 42 (reference)
#   py_out/<d>.seed<s>.all.tsv  Python, RANSAC seed s (yardstick)
#   r_out/<d>.numpy.all.tsv     R port with NumPy's RANSAC subsets injected
#   r_out/<d>.r_rng.all.tsv     R port as installed (R random numbers)
library(CroCoDeELR)
dir <- "/home/claude/data"
datasets <- c("toy", "PRJEB83730", "PRJNA698986_P2", "PRJNA698986_P3",
              "PRJEB12449", "PRJEB10878", "PRJNA352475", "PRJEB6337")
seeds <- 1:3

read_all <- function(f) {
  d <- read_conta_events(f)
  d$key <- paste(d$source, d$target, sep = "\r")
  d
}

num_samples <- function(d) {
  f <- if (d == "toy") {
    system.file("extdata", "mgs_profiles_test.tsv", package = "CroCoDeELR")
  } else {
    file.path(dir, paste0(d, ".meteor.tab"))
  }
  ncol(read_ab_table(f))
}

# Event-level comparison of run `b` against reference run `a`.
compare_runs <- function(a, b, cutoff = 0.5) {
  a_ev <- a[a$probability >= cutoff, ]
  b_ev <- b[b$probability >= cutoff, ]
  shared <- intersect(a_ev$key, b_ev$key)
  m <- merge(a_ev[a_ev$key %in% shared, ], b_ev[b_ev$key %in% shared, ],
             by = "key", suffixes = c(".a", ".b"))
  dprob <- abs(m$probability.b - m$probability.a)
  drate <- abs(m$rate.b - m$rate.a) / m$rate.a
  n_union <- length(union(a_ev$key, b_ev$key))
  c(
    events_a = nrow(a_ev),
    events_b = nrow(b_ev),
    shared = length(shared),
    only_a = nrow(a_ev) - length(shared),
    only_b = nrow(b_ev) - length(shared),
    jaccard = if (n_union) length(shared) / n_union else 1,
    max_dprob = if (length(dprob)) max(dprob) else 0,
    median_drate = if (length(drate)) median(drate) else 0,
    p95_drate = if (length(drate)) unname(quantile(drate, 0.95)) else 0,
    max_drate = if (length(drate)) max(drate) else 0
  )
}

rows <- list()
for (d in datasets) {
  py <- read_all(file.path(dir, "py_out", paste0(d, ".all.tsv")))
  rn <- read_all(file.path(dir, "r_out", paste0(d, ".numpy.all.tsv")))
  rr <- read_all(file.path(dir, "r_out", paste0(d, ".r_rng.all.tsv")))

  # 1. NumPy subsets injected: every pair must be identical to Python.
  m1 <- merge(py, rn, by = "key", suffixes = c(".py", ".r"))
  numpy_identical <- nrow(py) == nrow(rn) && nrow(m1) == nrow(py) &&
    all(m1$contamination_specific_species.py ==
          m1$contamination_specific_species.r) &&
    all(CroCoDeELR:::.fmt_rate(m1$rate.py) ==
          CroCoDeELR:::.fmt_rate(m1$rate.r)) &&
    all(m1$probability.py == m1$probability.r)

  # 2. R random numbers vs Python (seed 42).
  r_vs_py <- compare_runs(py, rr)

  # 3. Python with other seeds vs Python (seed 42).
  seed_cmp <- sapply(seeds, function(s) {
    compare_runs(py, read_all(file.path(dir, "py_out",
                                        paste0(d, ".seed", s, ".all.tsv"))))
  })

  rows[[d]] <- data.frame(
    dataset = d,
    samples = num_samples(d),
    pairs_with_line = nrow(py),
    numpy_identical = numpy_identical,
    events_py = r_vs_py[["events_a"]],
    events_r = r_vs_py[["events_b"]],
    shared = r_vs_py[["shared"]],
    only_py = r_vs_py[["only_a"]],
    only_r = r_vs_py[["only_b"]],
    jaccard_r = r_vs_py[["jaccard"]],
    jaccard_pyseeds_min = min(seed_cmp["jaccard", ]),
    jaccard_pyseeds_mean = mean(seed_cmp["jaccard", ]),
    max_dprob_r = r_vs_py[["max_dprob"]],
    max_dprob_pyseeds = max(seed_cmp["max_dprob", ]),
    median_drate_r = r_vs_py[["median_drate"]],
    median_drate_pyseeds = mean(seed_cmp["median_drate", ]),
    p95_drate_r = r_vs_py[["p95_drate"]],
    p95_drate_pyseeds = mean(seed_cmp["p95_drate", ]),
    max_drate_r = r_vs_py[["max_drate"]],
    max_drate_pyseeds = max(seed_cmp["max_drate", ]),
    stringsAsFactors = FALSE
  )
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
write.table(res, file.path(dir, "comparison.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
print(res, digits = 3)
