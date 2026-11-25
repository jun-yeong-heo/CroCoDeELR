# Validation helpers (not part of the package).
inject_numpy_subsets <- function(ns, file, cache) {
  s <- utils::read.delim(file, header = FALSE, col.names = c("n", "i", "j"))
  for (n in unique(s$n)) {
    m <- as.matrix(s[s$n == n, c("i", "j")]); dimnames(m) <- NULL
    assign(paste(n, 30L, 2L, 42L, sep = "_"), m, envir = cache)
  }
  invisible(TRUE)
}
r_features_all_pairs <- function(ab_file, ns) {
  ab <- ns$.prepare_ab_table(ns$read_ab_table(ab_file))
  model <- ns$load_model()
  out <- list()
  for (s in colnames(ab)) for (t in colnames(ab)) {
    if (s == t) next
    f <- ns$.extract_features(ab, s, t)
    if (is.null(f)) next
    out[[length(out) + 1]] <- data.frame(source = s, target = t,
      f1 = f$values[1], f2 = f$values[2], f3 = f$values[3], f4 = f$values[4],
      f5 = f$values[5], f6 = f$values[6], offset = f$conta_line_offset,
      species = paste(f$conta_line_species, collapse = ","), stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  res$probability <- ns$.predict_conta_proba(model, as.matrix(res[, paste0("f", 1:6)]))
  res
}
read_py_features <- function(file) {
  d <- utils::read.delim(file, header = FALSE, colClasses = c("character", "character", rep("numeric", 8), "character"),
                         col.names = c("source", "target", paste0("f", 1:6), "offset", "probability", "species"))
  d
}
