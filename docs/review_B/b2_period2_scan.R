# B2 EXPLORATORY (not used by the pipeline): which period-2 date set could give the
# paper's Table 2 period-2 numbers? Our h = 1 values equal the paper to 4 decimals for
# the full period and period 1, but not for period 2.
# One free boundary is fitted on ONE cell (fBm, h = 1); the other 39 period-2 cells
# are then out-of-sample checks of that hypothesis. Also checks April 11 vs 12.
# Usage: Rscript docs/review_B/b2_period2_scan.R <forecasts rds> <out csv>
source("R/code_b_config.R"); source_code_b()
args <- commandArgs(trailingOnly = TRUE)
fc <- readRDS(args[1]); out <- args[2]
pp <- read.csv("docs/paper_values/code_b_forecast_tables.csv")
pp <- pp[pp$table == "Table 2" & pp$panel == "period2", ]
mfbm <- code_b_classes()$mfbm
fc <- fc[fc$model %in% mfbm, ]
fc$l <- 100 * loss_se(fc$actual, fc$forecast)
cell_means <- function(start, end) {
  x <- fc[fc$target_date >= as.Date(start) & fc$target_date <= as.Date(end), ]
  a <- aggregate(l ~ model + h, x, mean)
  m <- merge(a, pp[, c("model", "h", "value")], by = c("model", "h"))
  m$rel_diff <- m$l / m$value - 1
  m
}
fit_cell <- function(m) m$rel_diff[m$model == "fBm" & m$h == 1]
summ <- function(label, start, end) {
  m <- cell_means(start, end)
  others <- m[!(m$model == "fBm" & m$h == 1), ]
  data.frame(label = label, start = start, end = end, fit_cell_rel_diff = fit_cell(m),
             max_abs_rel_diff_other_39 = max(abs(others$rel_diff)),
             median_abs_rel_diff_other_39 = stats::median(abs(others$rel_diff)),
             n_within_0.1pct_of_40 = sum(abs(m$rel_diff) < 0.001))
}
ends <- sort(unique(fc$target_date[fc$target_date >= as.Date("2018-01-01")]))
starts <- sort(unique(fc$target_date[fc$target_date >= as.Date("2016-01-01") & fc$target_date <= as.Date("2018-12-31")]))
x1 <- fc[fc$model == "fBm" & fc$h == 1, ]
x1 <- x1[order(x1$target_date), ]
target <- pp$value[pp$model == "fBm" & pp$h == 1]
mean_between <- function(s, e) mean(x1$l[x1$target_date >= s & x1$target_date <= e])
best_end <- ends[which.min(abs(vapply(ends, function(e) mean_between(as.Date("2017-04-12"), e), 0) - target))]
best_start <- starts[which.min(abs(vapply(starts, function(s) mean_between(s, as.Date("2021-07-30")), 0) - target))]
res <- rbind(
  summ("stated (Table 2 caption): 2017-04-12 to 2021-07-30", "2017-04-12", "2021-07-30"),
  summ("prose: 2017-04-11 to 2021-07-30", "2017-04-11", "2021-07-30"),
  summ("EXPLORATORY end fitted on fBm h=1", "2017-04-12", format(best_end)),
  summ("EXPLORATORY start fitted on fBm h=1", format(best_start), "2021-07-30"))
print(res, digits = 4, row.names = FALSE)
utils::write.csv(res, out, row.names = FALSE)
