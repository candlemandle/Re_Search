# B3: export full-sample loss matrices and bootstrap indices for an independent MCS check.
# Usage: Rscript docs/review_B/mcs_export.R <forecasts dir> <file suffix> <out dir>
# For every cell: losses (CSV), our p-values with the pipeline settings (B = 5000,
# block = 20, seed = cfg$seed + h) for 5 seeds, and for one seed with B = 1000 the
# exact bootstrap index matrix used by mcs_pvalues() (int32, row-major, B x n).
source("R/code_b_config.R"); source_code_b()
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; sfx <- if (length(args) > 1) args[2] else ""; out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cfg <- code_b_config("full"); cl <- code_b_classes()
cells <- data.frame(
  sample = c("dj30", "dj30", "dj30", "dj30", "dj30", "dj30", "mag7", "mag7"),
  table = c("Table 2", "Table 2", "Table 2", "Table 3", "Table 3", "Table 18", "Table 17", "Table 19"),
  class = c("mfbm", "mfbm", "mfbm", "vhar", "vhar", "vhar", "vhar", "vhar"),
  period = "full", metric = c("MSFE", "MSFE", "MSFE", "MSFE", "MSFE", "QLIKE", "MSFE", "QLIKE"),
  h = c(1L, 5L, 20L, 1L, 20L, 15L, 20L, 10L))
fcs <- list(dj30 = readRDS(file.path(dir, paste0("code_b_forecasts_dj30", sfx, ".rds"))),
            mag7 = readRDS(file.path(dir, paste0("code_b_forecasts_mag7", sfx, ".rds"))))
res <- list()
for (k in seq_len(nrow(cells))) {
  c0 <- cells[k, ]; id <- sprintf("cell%d", k)
  per <- (if (c0$sample == "dj30") cfg$dj30_periods else cfg$mag7_periods)[[c0$period]]
  Lm <- class_loss_matrix(fcs[[c0$sample]], cl[[c0$class]], per, c0$metric, c0$h)
  utils::write.csv(Lm, file.path(out, paste0(id, "_losses.csv")))
  pv <- sapply(0:4, function(r) { set.seed(cfg$seed + c0$h + 1000L * r); mcs_pvalues(Lm, 5000L, 20L) })
  set.seed(cfg$seed + c0$h); p_pipeline <- mcs_pvalues(Lm, 5000L, 20L)
  # exact resamples: replay the RNG stream of mcs_pvalues() with B = 1000
  set.seed(7L); idx <- t(vapply(seq_len(1000L), function(b) as.integer(block_bootstrap_index(nrow(Lm), 20L)), integer(nrow(Lm))))
  writeBin(as.integer(t(idx)), file.path(out, paste0(id, "_idx.bin")), size = 4L)
  set.seed(7L); p_exact <- mcs_pvalues(Lm, 1000L, 20L)
  res[[k]] <- data.frame(cell = id, c0[rep(1, ncol(Lm)), ], model = colnames(Lm), n = nrow(Lm),
                         p_pipeline_seed = p_pipeline, p_seed_mean = rowMeans(pv),
                         p_seed_min = apply(pv, 1, min), p_seed_max = apply(pv, 1, max),
                         p_B1000_seed7 = p_exact, row.names = NULL)
}
utils::write.csv(do.call(rbind, res), file.path(out, "r_pvalues.csv"), row.names = FALSE)
cat("exported", nrow(cells), "cells to", out, "\n")
