toy_file <- function() {
  system.file("extdata", "mgs_profiles_test.tsv", package = "CroCoDeELR")
}

# Replace the R RANSAC draws by those of the Python CroCoDeEL (NumPy RNG)
# for the duration of the calling test.
local_numpy_subsets <- function(env = parent.frame()) {
  cache <- CroCoDeELR:::.ransac_cache
  saved <- as.list(cache)
  s <- utils::read.delim(test_path("fixtures", "numpy_ransac_subsets.tsv"),
                         header = FALSE, col.names = c("n", "i", "j"))
  rm(list = ls(cache), envir = cache)
  for (n in unique(s$n)) {
    m <- as.matrix(s[s$n == n, c("i", "j")])
    dimnames(m) <- NULL
    assign(paste(n, 30L, 2L, 42L, sep = "_"), m, envir = cache)
  }
  withr::defer({
    rm(list = ls(cache), envir = cache)
    for (k in names(saved)) assign(k, saved[[k]], envir = cache)
  }, envir = env)
  invisible(TRUE)
}
