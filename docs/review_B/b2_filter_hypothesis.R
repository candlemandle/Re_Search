# B2 EXPLORATORY (not used by the pipeline): could the paper's linear (V)HAR MSFE
# cells have been computed WITHOUT our positivity filter (raw OLS forecasts, RV <= 0
# allowed; MSFE is defined for them, QLIKE is not)?
# The hypothesis is tested on every linear HAR MSFE cell of Tables 3, 15 and 17, not
# fitted on one cell. Usage:
#   Rscript docs/review_B/b2_filter_hypothesis.R <forecasts dir> <suffix> <calendar> <out csv>
source("R/code_b_config.R"); source_code_b()
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; sfx <- args[2]; calendar <- args[3]; out <- args[4]
cfg <- code_b_config("full"); cl <- code_b_classes()
paper <- utils::read.csv("docs/paper_values/code_b_forecast_tables.csv")
jobs <- list(
  list(sample = "dj30", table = "Table 3", panel = "full", models = cl$vhar, type = "vhar"),
  list(sample = "dj30", table = "Table 15", panel = "full", models = cl$vharf[-1], type = "vharf"),
  list(sample = "mag7", table = "Table 17", panel = "vhar", models = cl$vhar, type = "vhar"),
  list(sample = "mag7", table = "Table 17", panel = "vharf", models = cl$vharf[-1], type = "vharf"))
rows <- list()
for (j in jobs) {
  fc <- readRDS(file.path(dir, paste0("code_b_forecasts_", j$sample, sfx, ".rds")))
  panel <- load_logvol_panel(j$sample)
  assets <- cfg[[paste0(j$sample, "_assets")]]
  per <- (if (j$sample == "dj30") cfg$dj30_periods else cfg$mag7_periods)$full
  sub <- if (calendar == "common_obs") panel[stats::complete.cases(panel[, assets]), ] else panel
  for (m in j$models) {
    k <- if (m == "HAR") 1L else as.integer(sub("\\D+", "", m))
    L <- as.matrix(sub[, assets[1:k]])
    x <- fc[fc$model == m, ]
    origins <- sort(unique(x$origin_date))
    raw <- do.call(rbind, parallel::mclapply(match(origins, sub$date), function(t) {
      f <- har_fit_window(L[1:t, , drop = FALSE], cfg$horizons, cfg$window, j$type,
                          na_avg = calendar == "trading")
      cbind(f, origin_date = sub$date[t])
    }, mc.cores = 4L))
    x <- merge(x, raw[, c("origin_date", "h", "forecast", "pred_raw", "filtered")],
               by = c("origin_date", "h"), suffixes = c("", "_refit"))
    stopifnot(isTRUE(all.equal(x$forecast, x$forecast_refit)))
    x <- x[x$target_date >= as.Date(per[1]) & x$target_date <= as.Date(per[2]), ]
    for (h in cfg$horizons) {
      z <- x[x$h == h, ]
      pv <- paper$value[paper$table == j$table & paper$panel == j$panel & paper$model == m & paper$h == h]
      rows[[length(rows) + 1]] <- data.frame(sample = j$sample, table = j$table, panel = j$panel,
        model = m, h = h, calendar = calendar, n = nrow(z), n_filtered = sum(z$filtered),
        paper = pv, pipeline = 100 * mean(loss_se(z$actual, z$forecast)),
        no_filter = 100 * mean(loss_se(z$actual, z$pred_raw)))
    }
  }
}
res <- do.call(rbind, rows)
res$rel_pipeline <- res$pipeline / res$paper - 1
res$rel_no_filter <- res$no_filter / res$paper - 1
utils::write.csv(res, out, row.names = FALSE)
act <- res$n_filtered > 0
cat("cells:", nrow(res), " with filter events:", sum(act), "\n")
cat("cells with filter events, |rel diff| < 1%: pipeline", sum(abs(res$rel_pipeline[act]) < .01),
    " no_filter", sum(abs(res$rel_no_filter[act]) < .01), "\n")
cat("median |rel diff| (cells with filter events): pipeline", round(median(abs(res$rel_pipeline[act])), 4),
    " no_filter", round(median(abs(res$rel_no_filter[act])), 4), "\n")
print(res[act, c("sample", "table", "model", "h", "n_filtered", "paper", "pipeline", "no_filter", "rel_pipeline", "rel_no_filter")],
      digits = 4, row.names = FALSE)
