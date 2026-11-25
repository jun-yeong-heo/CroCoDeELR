# Figures of the validation results for the README (man/figures/).
# Base R graphics. Run after the validation scripts listed in README.md.
library(CroCoDeELR)
dir <- "/home/claude/data"
out_dir <- "/home/claude/CroCoDeELR/man/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Palette (validated categorical slots 1-2, neutral text inks).
col_r <- "#257de0"      # R port
col_py <- "#f47230"     # Python, other seed / other model
col_grey <- "#9c9c9c"
col_text <- "#191919"
col_text2 <- "#666666"
alpha <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)

datasets <- c("toy", "PRJEB83730", "PRJNA698986_P2", "PRJNA698986_P3",
              "PRJEB12449", "PRJEB10878", "PRJNA352475", "PRJEB6337")
labels <- c("toy (14)", "PRJEB83730 (17)", "PRJNA698986_P2 (78)",
            "PRJNA698986_P3 (91)", "PRJEB12449 (110)", "PRJEB10878 (128)",
            "PRJNA352475 (182)", "PRJEB6337 (237)")

rd <- function(f) {
  d <- read_conta_events(f)
  d$key <- paste(d$source, d$target, sep = "\r")
  d
}
setup_par <- function() {
  graphics::par(family = "sans", col.axis = col_text2, col.lab = col_text,
                col.main = col_text, fg = col_text2, las = 1, tcl = -0.25,
                mgp = c(2.4, 0.5, 0), cex.axis = 0.85, cex.lab = 0.95,
                font.main = 1, cex.main = 1, bty = "l")
}
panel_title <- function(letter, text) {
  graphics::mtext(paste0(letter, "  ", text), side = 3, line = 0.8, adj = 0,
                  cex = 0.95, col = col_text, font = 2)
}
grid_lines <- function(at_x = NULL, at_y = NULL) {
  if (!is.null(at_x)) graphics::abline(v = at_x, col = "#ececec", lwd = 1)
  if (!is.null(at_y)) graphics::abline(h = at_y, col = "#ececec", lwd = 1)
}
identity_panel <- function(x, y, xlab, ylab, letter, title, note, log = FALSE) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  lim <- range(c(x, y))
  graphics::plot(NA, xlim = lim, ylim = lim, xlab = xlab, ylab = ylab,
                 axes = FALSE)
  ticks <- if (log) seq(floor(lim[1]), ceiling(lim[2])) else pretty(lim)
  grid_lines(ticks, ticks)
  graphics::abline(0, 1, col = col_grey, lwd = 1)
  graphics::points(x, y, pch = 16, cex = 0.45, col = alpha(col_r, 0.35))
  if (log) {
    lab <- parse(text = paste0("10^", ticks))
    graphics::axis(1, at = ticks, labels = lab)
    graphics::axis(2, at = ticks, labels = lab)
  } else {
    graphics::axis(1, at = ticks)
    graphics::axis(2, at = ticks)
  }
  panel_title(letter, title)
  graphics::legend("topleft", legend = note, bty = "n", text.col = col_text,
                   cex = 0.85)
}

# ---------------------------------------------------------------------------
# Figure 1. Exactness without the random draws.
model_sk <- scan(file.path(dir, "model_check_sklearn.txt"), quiet = TRUE)
model_r <- scan(file.path(dir, "model_check_r.txt"), quiet = TRUE)

py_all <- do.call(rbind, lapply(datasets, function(d)
  rd(file.path(dir, "py_out", paste0(d, ".all.tsv")))))
np_all <- do.call(rbind, lapply(datasets, function(d)
  rd(file.path(dir, "r_out", paste0(d, ".numpy.all.tsv")))))
m1 <- merge(py_all, np_all, by = c("key"), suffixes = c(".py", ".r"))

grDevices::png(file.path(out_dir, "validation-exactness.png"), width = 2400,
               height = 820, res = 200)
setup_par()
graphics::par(mfrow = c(1, 3), mar = c(4, 4.2, 3, 1))
identity_panel(model_sk, model_r, "scikit-learn probability",
               "CroCoDeELR probability", "a", "Default model, 20,000 inputs",
               sprintf("max |difference| = %g", max(abs(model_sk - model_r))))
identity_panel(m1$probability.py, m1$probability.r, "CroCoDeEL probability",
               "CroCoDeELR probability", "b",
               sprintf("NumPy draws, %s pairs", format(nrow(m1), big.mark = ",")),
               sprintf("max |difference| = %g",
                       max(abs(m1$probability.py - m1$probability.r))))
identity_panel(m1$rate.py, m1$rate.r, "CroCoDeEL rate", "CroCoDeELR rate",
               "c", "NumPy draws, contamination rate",
               sprintf("identical formatted rates: %d / %d",
                       sum(CroCoDeELR:::.fmt_rate(m1$rate.py) ==
                             CroCoDeELR:::.fmt_rate(m1$rate.r)), nrow(m1)),
               log = TRUE)
grDevices::dev.off()

# ---------------------------------------------------------------------------
# Figure 2. R random draws vs the seed dependence of CroCoDeEL.
cmp <- read.delim(file.path(dir, "comparison.tsv"))
jac <- function(a, b, cutoff = 0.5) {
  ea <- a$key[a$probability >= cutoff]
  eb <- b$key[b$probability >= cutoff]
  u <- union(ea, eb)
  if (length(u)) length(intersect(ea, eb)) / length(u) else 1
}
seed_jac <- sapply(datasets, function(d) {
  py <- rd(file.path(dir, "py_out", paste0(d, ".all.tsv")))
  vapply(1:3, function(s) jac(py, rd(file.path(dir, "py_out",
                                              paste0(d, ".seed", s, ".all.tsv")))),
         numeric(1))
})
mp <- read.delim(file.path(dir, "comparison_metaphlan.tsv"))
jac_r <- c(cmp$jaccard_r, mp["r_vs_py", "jaccard"])
jac_seed <- cbind(seed_jac, metaphlan = c(mp["pyseed1_vs_py", "jaccard"], NA, NA))
all_labels <- c(labels, "PRJNA763023 MetaPhlAn4 (900)")

# Relative rate differences on shared events, pooled over Meteor datasets.
rate_diffs <- function(b_file_fun) {
  unlist(lapply(datasets, function(d) {
    a <- rd(file.path(dir, "py_out", paste0(d, ".all.tsv")))
    b <- rd(b_file_fun(d))
    a <- a[a$probability >= 0.5, ]
    b <- b[b$probability >= 0.5, ]
    m <- merge(a, b, by = "key", suffixes = c(".a", ".b"))
    abs(m$rate.b - m$rate.a) / m$rate.a
  }))
}
dr_r <- rate_diffs(function(d) file.path(dir, "r_out", paste0(d, ".r_rng.all.tsv")))
dr_seed <- unlist(lapply(1:3, function(s) rate_diffs(function(d)
  file.path(dir, "py_out", paste0(d, ".seed", s, ".all.tsv")))))

pair_probs <- function(b_file_fun) {
  do.call(rbind, lapply(datasets, function(d) {
    a <- rd(file.path(dir, "py_out", paste0(d, ".all.tsv")))
    b <- rd(b_file_fun(d))
    merge(a, b, by = "key", suffixes = c(".a", ".b"))[, c("probability.a", "probability.b")]
  }))
}
pp_r <- pair_probs(function(d) file.path(dir, "r_out", paste0(d, ".r_rng.all.tsv")))
pp_seed <- pair_probs(function(d) file.path(dir, "py_out", paste0(d, ".seed1.all.tsv")))

prob_panel <- function(pp, ylab, letter, title, col) {
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
                 xlab = "CroCoDeEL probability (seed 42)", ylab = ylab)
  grid_lines(seq(0, 1, 0.25), seq(0, 1, 0.25))
  graphics::abline(0, 1, col = col_grey)
  graphics::abline(v = 0.5, h = 0.5, col = col_grey, lty = 3)
  graphics::points(pp$probability.a, pp$probability.b, pch = 16, cex = 0.45,
                   col = alpha(col, 0.3))
  graphics::axis(1, at = seq(0, 1, 0.25))
  graphics::axis(2, at = seq(0, 1, 0.25))
  disc <- sum((pp$probability.a >= 0.5) != (pp$probability.b >= 0.5))
  graphics::legend("topleft", bty = "n", cex = 0.85, text.col = col_text,
                   legend = c(sprintf("%s pairs", format(nrow(pp), big.mark = ",")),
                              sprintf("different call: %d", disc),
                              sprintf("same probability: %.1f%%",
                                      100 * mean(pp$probability.a == pp$probability.b))))
  panel_title(letter, title)
}

grDevices::png(file.path(out_dir, "validation-random-draws.png"), width = 2400,
               height = 2000, res = 200)
setup_par()
graphics::layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1.3, 1))
# a. Jaccard per dataset
graphics::par(mar = c(4, 13, 3, 1))
n <- length(all_labels)
y <- rev(seq_len(n))
graphics::plot(NA, xlim = c(0.75, 1), ylim = c(0.5, n + 0.5), axes = FALSE,
               xlab = "Event Jaccard index vs CroCoDeEL (seed 42)", ylab = "")
grid_lines(seq(0.75, 1, 0.05), NULL)
for (i in seq_len(n)) {
  s <- jac_seed[, i]
  s <- s[!is.na(s)]
  graphics::segments(min(s, jac_r[i]), y[i], max(s, jac_r[i]), y[i],
                     col = "#dddddd", lwd = 2)
  graphics::points(s, rep(y[i], length(s)), pch = 21, cex = 1.3, lwd = 1.5,
                   col = col_py, bg = "white")
  graphics::points(jac_r[i], y[i], pch = 16, cex = 1.4, col = col_r)
}
graphics::axis(1, at = seq(0.75, 1, 0.05))
graphics::axis(2, at = y, labels = all_labels, tick = FALSE, cex.axis = 0.85)
panel_title("a", "Detected events (probability >= 0.5)")
graphics::legend("topleft", bty = "n", pch = c(16, 21), pt.cex = c(1.4, 1.3),
                 pt.lwd = c(1, 1.5), col = c(col_r, col_py), pt.bg = "white",
                 legend = c("CroCoDeELR (R draws)",
                            "CroCoDeEL seed 1, 2, 3"),
                 text.col = col_text, cex = 0.85)
# b. ECDF of rate differences
graphics::par(mar = c(4, 4.2, 3, 1))
graphics::plot(NA, xlim = c(0, 0.25), ylim = c(0, 1), axes = FALSE,
               xlab = "Relative rate difference vs seed 42",
               ylab = "Fraction of shared events")
grid_lines(seq(0, 0.25, 0.05), seq(0, 1, 0.25))
gx <- seq(0, 0.25, length.out = 500)
graphics::lines(gx, stats::ecdf(dr_seed)(gx), col = col_py, lwd = 2, type = "s")
graphics::lines(gx, stats::ecdf(dr_r)(gx), col = col_r, lwd = 2, type = "s")
graphics::axis(1, at = seq(0, 0.25, 0.05), labels = paste0(seq(0, 25, 5), "%"))
graphics::axis(2, at = seq(0, 1, 0.25))
graphics::legend("bottomright", bty = "n", lwd = 2, col = c(col_r, col_py),
                 legend = c("CroCoDeELR", "CroCoDeEL seeds 1-3"),
                 text.col = col_text, cex = 0.85)
panel_title("b", "Rates of shared events (Meteor)")
# c, d. pair-level probabilities
graphics::par(mar = c(4, 4.2, 3, 1))
prob_panel(pp_r, "CroCoDeELR probability", "c",
           "Every pair (Meteor), R draws", col_r)
prob_panel(pp_seed, "CroCoDeEL probability (seed 1)", "d",
           "Every pair (Meteor), seed 1", col_py)
grDevices::dev.off()

# ---------------------------------------------------------------------------
# Figure 3. Model training on the official CroCoDeEL training dataset.
tr <- file.path(dir, "train")
py_tr <- read.delim(file.path(tr, "py_train_eval.tsv"))
r_tr <- read.delim(file.path(tr, "r_train_eval.tsv"))
py_tr$is_contaminated <- py_tr$is_contaminated == "True"
common <- intersect(py_tr$source[py_tr$set == "test"], r_tr$source[r_tr$set == "test"])
pyc <- py_tr[match(common, py_tr$source), ]
rc <- r_tr[match(common, r_tr$source), ]

meta <- read.delim(file.path(tr, "training_dataset.metadata.tsv"))
meta <- meta[startsWith(meta$source_name, "conta"), ]
ab <- CroCoDeELR:::.prepare_ab_table(read_ab_table(file.path(tr, "training_dataset.meteor.tsv")))
r_rate <- vapply(seq_len(nrow(meta)), function(i) {
  f <- CroCoDeELR:::.extract_features(ab, meta$source_name[i], meta$target_name[i])
  if (is.null(f)) NA_real_ else 10^(-f$conta_line_offset)
}, numeric(1))
ok <- !is.na(r_rate) & !is.na(meta$crocodeel_conta_rate) & meta$crocodeel_conta_rate > 0

tm <- read.delim(file.path(dir, "comparison_trained_models.tsv"))

grDevices::png(file.path(out_dir, "validation-training.png"), width = 2400,
               height = 820, res = 200)
setup_par()
graphics::par(mfrow = c(1, 3), mar = c(4, 4.2, 3, 1))
# a. held-out probabilities
graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
               xlab = "CroCoDeEL (scikit-learn) probability",
               ylab = "CroCoDeELR (ranger) probability")
grid_lines(seq(0, 1, 0.25), seq(0, 1, 0.25))
graphics::abline(0, 1, col = col_grey)
graphics::abline(v = 0.5, h = 0.5, col = col_grey, lty = 3)
graphics::points(pyc$prob_retrained[!pyc$is_contaminated],
                 rc$prob_retrained[!rc$is_contaminated], pch = 21, cex = 0.8,
                 col = col_py, bg = "white")
graphics::points(pyc$prob_retrained[pyc$is_contaminated],
                 rc$prob_retrained[rc$is_contaminated], pch = 16, cex = 0.6,
                 col = alpha(col_r, 0.5))
graphics::axis(1, at = seq(0, 1, 0.25))
graphics::axis(2, at = seq(0, 1, 0.25))
graphics::legend("left", bty = "n", pch = c(16, 21), col = c(col_r, col_py),
                 legend = c(sprintf("contaminated (%d)", sum(pyc$is_contaminated)),
                            sprintf("not contaminated (%d)", sum(!pyc$is_contaminated))),
                 text.col = col_text, cex = 0.85)
panel_title("a", sprintf("Pairs held out by both (%d)", nrow(pyc)))
# b. rates on contaminated training pairs
identity_panel(meta$crocodeel_conta_rate[ok], r_rate[ok], "CroCoDeEL rate",
               "CroCoDeELR rate", "b",
               sprintf("Contaminated pairs (%s)", format(sum(ok), big.mark = ",")),
               sprintf("median |log10 difference| = %.4f",
                       stats::median(abs(log10(r_rate[ok]) -
                                           log10(meta$crocodeel_conta_rate[ok])))),
               log = TRUE)
# c. real data
graphics::par(mar = c(4, 9, 3, 1))
n <- nrow(tm)
y <- rev(seq_len(n))
graphics::plot(NA, xlim = c(0.7, 1), ylim = c(0.5, n + 0.5), axes = FALSE,
               xlab = "Event Jaccard index", ylab = "")
grid_lines(seq(0.7, 1, 0.1), NULL)
graphics::points(tm$jaccard_sklearn_vs_default, y, pch = 21, cex = 1.3,
                 lwd = 1.5, col = col_py, bg = "white")
graphics::points(tm$jaccard_ranger_vs_sklearn, y, pch = 16, cex = 1.4, col = col_r)
graphics::axis(1, at = seq(0.7, 1, 0.1))
graphics::axis(2, at = y, labels = tm$dataset, tick = FALSE, cex.axis = 0.85)
graphics::legend("topleft", bty = "n", pch = c(16, 21), pt.lwd = c(1, 1.5),
                 col = c(col_r, col_py), pt.bg = "white",
                 legend = c("ranger vs scikit-learn", "scikit-learn vs default"),
                 text.col = col_text, cex = 0.85)
panel_title("c", "Retrained models on real data")
grDevices::dev.off()
