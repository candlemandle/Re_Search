# B2: row-level diagnosis of the five large discrepancies named by the review.
# Usage: Rscript docs/review_B/b2_diagnose.R <forecasts dir> <suffix> <calendar> <out prefix>
#   e.g. Rscript docs/review_B/b2_diagnose.R results/forecasts "" trading docs/review_B/b2_trading
# Reads saved forecasts, refits only the HAR model of each cell to recover the raw
# OLS prediction (before the positivity filter). Writes <out prefix>_cells.csv,
# _top_losses.csv, _variants.csv. Variants change ONE factor at a time; the filter
# variants are an exploratory sensitivity analysis, not a tuning of the pipeline.
source("R/code_b_config.R"); source_code_b()
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; sfx <- args[2]; calendar <- args[3]; out <- args[4]
cfg <- code_b_config("full"); cl <- code_b_classes()
paper <- utils::read.csv("docs/paper_values/code_b_forecast_tables.csv")
cells <- data.frame(
  sample = c("mag7", "mag7", "dj30", "dj30", "dj30"),
  table = c("Table 17", "Table 19", "Table 18", "Table 3", "Table 2"),
  panel = c("vhar", "vhar", "vhar", "full", "period2"),
  class = c("vhar", "vhar", "vhar", "vhar", "mfbm"),
  period = c("full", "full", "full", "full", "period2"),
  metric = c("MSFE", "QLIKE", "QLIKE", "MSFE", "MSFE"),
  model = c("VHAR5", "VHAR3", "VHAR4", "VHAR4", "fBm"),
  h = c(20L, 10L, 15L, 20L, 1L), stringsAsFactors = FALSE)
fcs <- list(dj30 = readRDS(file.path(dir, paste0("code_b_forecasts_dj30", sfx, ".rds"))),
            mag7 = readRDS(file.path(dir, paste0("code_b_forecasts_mag7", sfx, ".rds"))))
panels <- list(dj30 = load_logvol_panel("dj30"), mag7 = load_logvol_panel("mag7"))
loss <- function(metric, a, f) if (metric == "MSFE") 100 * loss_se(a, f) else loss_qlike(a, f)

# raw HAR fits at the saved origins of one model and h, plus window statistics
refit_har <- function(panel, assets, model, h, origin_dates, window, extra = 0L) {
  k <- if (model %in% c("HAR", "LHAR")) 1L else as.integer(sub("\\D+", "", model))
  type <- if (grepl("^VHARF", model)) "vharf" else "vhar"
  log <- grepl("LHAR", model)
  sub <- if (calendar == "common_obs") panel[stats::complete.cases(panel[, assets]), ] else panel
  L <- as.matrix(sub[, assets[1:k]])
  t_all <- match(origin_dates, sub$date)
  do.call(rbind, lapply(t_all, function(t) {
    # the 500-pair variant needs h more rows; the first origins use what exists
    f <- har_fit_window(L[1:t, , drop = FALSE], h, min(window + extra, t), type, log, na_avg = calendar == "trading")
    y <- exp(L[(t - window + 1):t, 1])
    cbind(f, origin_date = sub$date[t], y_min = min(y, na.rm = TRUE), y_max = max(y, na.rm = TRUE),
          y_mean = mean(y, na.rm = TRUE))
  }))
}

cell_rows <- top_rows <- var_rows <- list()
for (k in seq_len(nrow(cells))) {
  c0 <- cells[k, ]
  assets <- cfg[[paste0(c0$sample, "_assets")]]
  per <- (if (c0$sample == "dj30") cfg$dj30_periods else cfg$mag7_periods)[[c0$period]]
  fc <- fcs[[c0$sample]]
  x <- fc[fc$model == c0$model & fc$h == c0$h & fc$target_date >= as.Date(per[1]) &
          fc$target_date <= as.Date(per[2]), ]
  x <- x[order(x$target_date), ]
  x$loss <- loss(c0$metric, x$actual, x$forecast)
  pv <- paper$value[paper$table == c0$table & paper$panel == c0$panel & paper$model == c0$model & paper$h == c0$h]
  ord <- order(-x$loss)
  share <- function(m) sum(x$loss[ord[seq_len(m)]]) / sum(x$loss)
  row <- data.frame(c0, calendar = calendar, n = nrow(x), first_target = min(x$target_date),
    last_target = max(x$target_date), ours = mean(x$loss), paper = pv, rel_diff = mean(x$loss) / pv - 1,
    top1_share = share(1), top10_share = share(10), top1pct_share = share(ceiling(nrow(x) / 100)),
    mean_without_top10 = mean(x$loss[-ord[1:10]]),
    min_forecast = min(x$forecast), n_forecast_below_0.05 = sum(x$forecast < 0.05),
    n_filtered = NA_integer_, n_below_window_min = NA_integer_, n_rank_deficient = NA_integer_)
  top <- head(x[ord, c("origin_date", "target_date", "forecast", "actual", "loss")], 10)
  top$share <- top$loss / sum(x$loss)
  if (c0$class != "mfbm") {
    raw <- refit_har(panels[[c0$sample]], assets, c0$model, c0$h, x$origin_date, cfg$window)
    stopifnot(all.equal(raw$forecast, x$forecast)) # the refit reproduces the saved forecasts
    row$n_filtered <- sum(raw$filtered)
    row$n_below_window_min <- sum(raw$pred_raw < raw$y_min)
    row$n_rank_deficient <- sum(raw$rank < raw$n_coef)
    top$pred_raw <- raw$pred_raw[ord[1:10]]
    top$filtered <- raw$filtered[ord[1:10]]
    top$window_min <- raw$y_min[ord[1:10]]
    p <- raw$pred_raw
    variants <- list(
      pipeline_filter_le0_to_window_min = x$forecast,
      exploratory_clip_to_window_range = pmin(pmax(p, raw$y_min), raw$y_max),
      exploratory_outside_range_to_window_mean = ifelse(p < raw$y_min | p > raw$y_max, raw$y_mean, p))
    if (c0$metric == "MSFE") variants$exploratory_no_filter <- p
    # window endpoint: 500 regression pairs instead of 500 - h (pairs reach h days further back)
    raw500 <- refit_har(panels[[c0$sample]], assets, c0$model, c0$h, x$origin_date, cfg$window, extra = c0$h)
    variants$endpoint_500_pairs <- raw500$forecast
    for (v in names(variants)) var_rows[[length(var_rows) + 1]] <- data.frame(c0, calendar = calendar,
      variant = v, value = mean(loss(c0$metric, x$actual, variants[[v]])), paper = pv,
      rel_diff = mean(loss(c0$metric, x$actual, variants[[v]])) / pv - 1)
  }
  cell_rows[[k]] <- row
  top_rows[[k]] <- cbind(c0[rep(1, nrow(top)), c("sample", "table", "model", "h")], top, row.names = NULL)
}
utils::write.csv(do.call(rbind, cell_rows), paste0(out, "_cells.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(top_rows, function(z) { z[setdiff(c("pred_raw", "filtered", "window_min"), names(z))] <- NA; z })),
                 paste0(out, "_top_losses.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, var_rows), paste0(out, "_variants.csv"), row.names = FALSE)
print(do.call(rbind, cell_rows)[, c("table", "model", "h", "period", "n", "ours", "paper", "rel_diff",
  "top1_share", "top10_share", "mean_without_top10", "min_forecast", "n_filtered", "n_below_window_min", "n_rank_deficient")], digits = 4)
print(do.call(rbind, var_rows)[, c("table", "model", "h", "variant", "value", "paper", "rel_diff")], digits = 4)
