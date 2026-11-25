# CroCoDeELR on the MetaPhlAn4 dataset (900 samples, --filter-low-ab 20).
library(CroCoDeELR)
ab <- read_ab_table("/home/claude/data/PRJNA763023_PRJDB4176.metaphlan4.tab")
t0 <- proc.time()[["elapsed"]]
ev <- search_conta(ab, filter_low_ab = 20, nproc = 2, verbose = FALSE)
el <- proc.time()[["elapsed"]] - t0
write_conta_events(ev, "/home/claude/data/r_out/metaphlan.r_rng.tsv")
cat("metaphlan r_rng", ncol(ab), "samples", sprintf("%.1f", el), "s\n",
    file = "/home/claude/data/r_out/times.txt", append = TRUE)
