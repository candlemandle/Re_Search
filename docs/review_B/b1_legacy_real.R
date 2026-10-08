# B1 on real data: realised trading-day distance of legacy (common_obs) forecasts.
# Input: legacy full forecasts saved by Code B (inspected artifacts, not rerun here).
# Usage: Rscript docs/review_B/b1_legacy_real.R <dir with code_b_forecasts_{dj30,mag7}.rds>
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results/forecasts"
out <- list()
for (nm in c("dj30", "mag7")) {
  fc <- readRDS(file.path(dir, sprintf("code_b_forecasts_%s.rds", nm)))
  p <- read.csv(sprintf("data/processed/%s_logvol.csv", nm)); p$date <- as.Date(p$date)
  pos <- match(fc$target_date, p$date) - match(fc$origin_date, p$date)
  fc$trading_days <- pos
  x <- unique(fc[fc$model == "fBm", c("h", "origin_date", "target_date", "trading_days")])
  per_h <- do.call(rbind, lapply(split(x, x$h), function(z) data.frame(sample = nm, h = z$h[1],
    n_forecasts = nrow(z), n_shifted = sum(z$trading_days != z$h),
    max_extra_days = max(z$trading_days - z$h),
    first_target = min(z$target_date), last_target = max(z$target_date))))
  out[[nm]] <- per_h
  cat("\n==", nm, ": origins", length(unique(fc$origin_date)), " first origin", format(min(fc$origin_date)),
      " models", length(unique(fc$model)), "\n")
}
res <- do.call(rbind, out); rownames(res) <- NULL
print(res)
utils::write.csv(res, "docs/review_B/b1_legacy_shifted_targets.csv", row.names = FALSE)
