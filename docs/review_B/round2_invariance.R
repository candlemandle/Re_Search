# Second review: retained DJ30 VHAR4 forecasts (h = 15, 20) after the fix versus the
# saved full runs of all eight horizons. Usage:
#   Rscript --vanilla docs/review_B/round2_invariance.R <after rds> <trading full rds> <common_obs full rds>
a <- commandArgs(trailingOnly = TRUE)
after <- readRDS(a[1])
saved <- list(trading = readRDS(a[2]), common_obs = readRDS(a[3]))
cols <- c("model", "h", "origin_date", "target_date", "forecast", "actual", "train_start", "n_train_days")
key <- function(x) paste(x$h, x$origin_date, x$target_date)
for (cal in names(saved)) {
  x <- after[[cal]]
  y <- saved[[cal]]; y <- y[y$model == "VHAR4" & y$h %in% c(15L, 20L), ]
  x <- x[order(key(x)), ]; y <- y[order(key(y)), ]
  cc <- intersect(cols, intersect(names(x), names(y)))
  cat(sprintf("%-10s rows after %d, saved %d | identical keys %s | identical %s | max|forecast diff| %g\n",
    cal, nrow(x), nrow(y), identical(key(x), key(y)),
    paste(cc[mapply(function(u, v) identical(u, v), x[cc], y[cc])], collapse = ","),
    max(abs(x$forecast - y$forecast))))
}
