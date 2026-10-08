# forecast losses and the model confidence set (appendix F.4, F.5).

# squared forecast error. the paper reports its mean times 100
loss_se <- function(actual, forecast) (actual - forecast)^2

# QLIKE, appendix F.5: RV / RV_hat - log(RV / RV_hat) - 1
loss_qlike <- function(actual, forecast) {
  r <- actual / forecast
  r - log(r) - 1
}

# indices of a circular moving block bootstrap sample of length n
block_bootstrap_index <- function(n, block) {
  starts <- sample.int(n, ceiling(n / block), replace = TRUE)
  idx <- as.vector(outer(0:(block - 1), starts, "+"))
  ((idx - 1) %% n + 1)[1:n]
}

# MCS p-values of hansen, lunde and nason (2011) with the T_max statistic.
#   L     : n x M matrix of losses (rows = dates, columns = models)
#   B     : bootstrap replications (5000 in the paper)
#   block : block length (20 in the paper)
# at every step we test "all models left are equally good". if rejected we drop
# the model with the largest t statistic. the p-value of a dropped model is the
# largest test p-value seen up to that step, the last model gets 1
# when block >= n every circular block resample is a rotation of the whole sample,
# so every bootstrap mean equals the sample mean and the test has no information
# (all p-values would be 1). those p-values are returned as NA, see mcs_status()
mcs_pvalues <- function(L, B = 5000L, block = 20L) {
  L <- as.matrix(L)
  if (!all(is.finite(L))) stop("mcs_pvalues: losses must be finite")
  n <- nrow(L)
  if (mcs_status(n, block) != "ok") return(rep(NA_real_, ncol(L)))
  # the same bootstrap samples are used in every step. we keep them as weights:
  # W[b, t] = (times day t is in sample b) / n, so a bootstrap mean is W %*% d
  W <- matrix(0, B, n)
  for (b in seq_len(B)) W[b, ] <- tabulate(block_bootstrap_index(n, block), n) / n

  alive <- seq_len(ncol(L))
  pval <- rep(NA_real_, ncol(L))
  p_max <- 0
  while (length(alive) > 1) {
    # loss of each model relative to the average of the models still alive
    d <- L[, alive, drop = FALSE] - rowMeans(L[, alive, drop = FALSE])
    d_bar <- colMeans(d)
    d_boot <- t(W %*% d) # models x B
    se <- sqrt(rowMeans((d_boot - d_bar)^2))
    se[se == 0] <- Inf # identical forecasts. they can never be the worst

    t_stat <- d_bar / se
    t_max_boot <- apply((d_boot - d_bar) / se, 2, max)
    p <- mean(t_max_boot >= max(t_stat))

    p_max <- max(p_max, p)
    worst <- alive[which.max(t_stat)]
    pval[worst] <- p_max
    alive <- setdiff(alive, worst)
  }
  pval[alive] <- 1
  pval
}

# is the block bootstrap MCS informative for n losses and this block length?
mcs_status <- function(n, block) {
  if (n < 2L) return("unavailable_n_lt_2")
  if (block >= n) return("unavailable_block_ge_n")
  "ok"
}

# loss table for one model class and one period: mean loss for every model and h
# plus the MCS p-values computed separately for every h.
#   fc      : forecasts from run_all_models()
#   models  : models of the class in the order of the table
#   metric  : "MSFE" (mean squared error x 100) or "QLIKE"
evaluate_class <- function(fc, models, period, metric = c("MSFE", "QLIKE"),
                           B = 5000L, block = 20L, seed = 1L) {
  metric <- match.arg(metric)
  scale <- if (metric == "MSFE") 100 else 1
  out <- list()
  in_period <- fc$target_date >= as.Date(period[1]) & fc$target_date <= as.Date(period[2])
  for (h in sort(unique(fc$h[in_period & fc$model %in% models]))) {
    Lm <- class_loss_matrix(fc, models, period, metric, h)
    set.seed(seed + h)
    out[[length(out) + 1]] <- data.frame(
      model = models, h = h, metric = metric,
      value = colMeans(Lm) * scale,
      mcs_p = mcs_pvalues(Lm, B, block),
      n_forecasts = nrow(Lm),
      first_target = min(as.Date(rownames(Lm))),
      last_target = max(as.Date(rownames(Lm))),
      mcs_status = mcs_status(nrow(Lm), block),
      mcs_B = B, mcs_block = block, mcs_seed = seed + h
    )
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

# target dates x models loss matrix of one class, period (by target date) and h.
# rows line up because all models share the dates; not scaled by 100
class_loss_matrix <- function(fc, models, period, metric = c("MSFE", "QLIKE"), h) {
  metric <- match.arg(metric)
  x <- fc[fc$target_date >= as.Date(period[1]) & fc$target_date <= as.Date(period[2]) &
          fc$model %in% models & fc$h == h, ]
  loss_fun <- if (metric == "MSFE") loss_se else loss_qlike
  x$loss <- loss_fun(x$actual, x$forecast)
  tapply(x$loss, list(as.character(x$target_date), x$model), sum)[, models, drop = FALSE]
}

# build several paper tables at once.
#   specs   : data.frame with columns table, panel, class, period, metric
#             (period is a name from `periods`, class a name from code_b_classes())
#   periods : named list of c(start, end) target dates
paper_tables <- function(fc, specs, periods, B, block, seed) {
  classes <- code_b_classes()
  out <- lapply(seq_len(nrow(specs)), function(k) {
    s <- specs[k, ]
    tab <- evaluate_class(fc, classes[[s$class]], periods[[s$period]], s$metric, B, block, seed)
    cbind(table = s$table, panel = s$panel, class = s$class, period = s$period, tab)
  })
  do.call(rbind, out)
}

# our tables next to the paper numbers (docs/paper_values/code_b_forecast_tables.csv)
compare_with_paper <- function(ours, paper) {
  paper <- paper[, c("table", "panel", "model", "h", "value", "mcs_p")]
  names(paper)[5:6] <- c("paper_value", "paper_mcs_p")
  cmp <- merge(ours, paper, by = c("table", "panel", "model", "h"), all.x = TRUE)
  cmp$diff <- cmp$value - cmp$paper_value
  cmp$rel_diff <- cmp$diff / cmp$paper_value
  cmp[order(cmp$table, cmp$panel, match(cmp$model, unlist(code_b_classes())), cmp$h), ]
}

# short "original vs ours": for every table, panel and h the loss of the largest
# model (5 assets) relative to the univariate model of the same class.
# ratio < 1 means the extra assets help. same_direction checks we agree with the paper
original_vs_ours <- function(cmp) {
  classes <- code_b_classes()
  keys <- unique(cmp[, c("table", "panel", "class", "h")])
  rows <- lapply(seq_len(nrow(keys)), function(k) {
    x <- merge(keys[k, ], cmp)
    uni <- classes[[keys$class[k]]][1]
    big <- classes[[keys$class[k]]][5]
    ours <- x$value[x$model == big] / x$value[x$model == uni]
    paper <- x$paper_value[x$model == big] / x$paper_value[x$model == uni]
    data.frame(keys[k, ], quantity = paste0(big, " / ", uni),
               paper = paper, ours = ours,
               same_direction = sign(paper - 1) == sign(ours - 1),
               best_paper = x$model[which.min(x$paper_value)],
               best_ours = x$model[which.min(x$value)])
  })
  do.call(rbind, rows)
}

# wide version for reading: one row per model, one column per h
to_wide <- function(tab, col = "value") {
  w <- stats::reshape(tab[, c("model", "h", col)], idvar = "model", timevar = "h",
                      direction = "wide")
  names(w) <- sub(paste0(col, "."), paste0(col, "_h"), names(w), fixed = TRUE)
  rownames(w) <- NULL
  w
}

# plot of original_vs_ours(): loss of the 5-asset model / univariate model by h,
# one colour per class, solid = ours, dotted = paper
plot_class_ratio <- function(ovo, file, main = "") {
  cls <- c(mfbm = "black", vhar = "red", vharf = "blue", vlhar = "darkgreen")
  lab <- c(mfbm = "mfBm5 / fBm", vhar = "VHAR5 / HAR", vharf = "VHARF5 / HAR", vlhar = "VLHAR5 / LHAR")
  grDevices::pdf(file, width = 6, height = 4.5)
  plot(NA, xlim = range(ovo$h), ylim = range(c(ovo$ours, ovo$paper, 1), na.rm = TRUE),
       xlab = "horizon h (days)", ylab = "loss of 5-asset model / univariate model", main = main)
  for (k in names(cls)) {
    x <- ovo[ovo$class == k, ]
    lines(x$h, x$ours, col = cls[k], type = "b", pch = 20)
    lines(x$h, x$paper, col = cls[k], lty = 3)
  }
  abline(h = 1, col = "grey")
  legend("topleft", c(lab, "ours", "paper"), col = c(cls, "grey30", "grey30"),
         lty = c(1, 1, 1, 1, 1, 3), bty = "n", cex = 0.8)
  grDevices::dev.off()
}

# --- calendar conventions (review B1) -------------------------------------------

# configuration of one calendar run; the list of conventions itself is not part of
# the run identity
code_b_calendar_config <- function(cfg, calendar) {
  cfg$calendar <- match.arg(calendar, code_b_calendars())
  cfg$calendars <- NULL
  cfg
}

# file suffix: the primary convention keeps the historical file names
code_b_calendar_suffix <- function(cfg, calendar) {
  if (identical(calendar, cfg$calendar)) "" else paste0("_", calendar)
}

# the same paper cells under both conventions, next to the paper value
calendar_sensitivity <- function(tabs, paper) {
  key <- c("table", "panel", "class", "period", "model", "h", "metric")
  cols <- c(key, "value", "mcs_p", "n_forecasts", "first_target", "last_target")
  m <- merge(tabs[["trading"]][, cols], tabs[["common_obs"]][, cols], by = key,
             suffixes = c("_trading", "_common_obs"))
  p <- paper[, c("table", "panel", "model", "h", "value")]
  names(p)[5] <- "paper_value"
  m <- merge(m, p, by = c("table", "panel", "model", "h"), all.x = TRUE)
  m$rel_diff_trading_vs_common_obs <- m$value_trading / m$value_common_obs - 1
  m$rel_diff_trading_vs_paper <- m$value_trading / m$paper_value - 1
  m$rel_diff_common_obs_vs_paper <- m$value_common_obs / m$paper_value - 1
  m[order(m$table, m$panel, match(m$model, unlist(code_b_classes())), m$h), ]
}

# forecast keys of the two conventions, per h (one model is enough: keys are shared).
#   shifted_targets  : common_obs keys whose target is not h trading days after the origin
#   same_origin_*    : origins used by both conventions
#   train_start_diff : of those, origins whose 500-row window starts on another day
calendar_key_audit <- function(fc_trading, fc_common, dates, model = "fBm") {
  a <- fc_trading[fc_trading$model == model, ]
  b <- fc_common[fc_common$model == model, ]
  key <- function(x) paste(x$h, x$origin_date, x$target_date)
  rows <- lapply(sort(unique(c(a$h, b$h))), function(h) {
    ah <- a[a$h == h, ]; bh <- b[b$h == h, ]
    dist <- match(bh$target_date, dates) - match(bh$origin_date, dates)
    so <- intersect(ah$origin_date, bh$origin_date)
    ts_a <- ah$train_start[match(so, ah$origin_date)]
    ts_b <- bh$train_start[match(so, bh$origin_date)]
    data.frame(h = h, n_trading = nrow(ah), n_common_obs = nrow(bh),
               n_identical_keys = length(intersect(key(ah), key(bh))),
               n_only_trading = length(setdiff(key(ah), key(bh))),
               n_only_common_obs = length(setdiff(key(bh), key(ah))),
               shifted_targets_common_obs = sum(dist != h),
               max_target_distance_common_obs = max(dist),
               same_origins = length(so),
               train_start_diff = if (is.null(ts_a) || is.null(ts_b)) NA_integer_ else sum(ts_a != ts_b))
  })
  do.call(rbind, rows)
}

# identity of one run: input checksum, code fingerprint and the full configuration
code_b_provenance <- function(script, inputs, cfg) {
  data.frame(script = script, time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             input = inputs, input_md5 = unname(tools::md5sum(inputs)),
             code_fingerprint = code_b_fingerprint(list()),
             config = paste(deparse(cfg, control = c("keepNA", "niceNames")), collapse = " "),
             r_version = R.version.string)
}
