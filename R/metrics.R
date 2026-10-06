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
mcs_pvalues <- function(L, B = 5000L, block = 20L) {
  L <- as.matrix(L)
  if (!all(is.finite(L))) stop("mcs_pvalues: losses must be finite")
  n <- nrow(L)
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

# loss table for one model class and one period: mean loss for every model and h
# plus the MCS p-values computed separately for every h.
#   fc      : forecasts from run_all_models()
#   models  : models of the class in the order of the table
#   metric  : "MSFE" (mean squared error x 100) or "QLIKE"
evaluate_class <- function(fc, models, period, metric = c("MSFE", "QLIKE"),
                           B = 5000L, block = 20L, seed = 1L) {
  metric <- match.arg(metric)
  in_period <- fc$target_date >= as.Date(period[1]) & fc$target_date <= as.Date(period[2])
  fc <- fc[in_period & fc$model %in% models, ]
  loss_fun <- if (metric == "MSFE") loss_se else loss_qlike
  scale <- if (metric == "MSFE") 100 else 1
  fc$loss <- loss_fun(fc$actual, fc$forecast)

  out <- list()
  for (h in sort(unique(fc$h))) {
    x <- fc[fc$h == h, ]
    # dates x models loss matrix. rows line up because all models share the dates
    Lm <- tapply(x$loss, list(as.character(x$target_date), x$model), sum)[, models, drop = FALSE]
    set.seed(seed + h)
    out[[length(out) + 1]] <- data.frame(
      model = models, h = h, metric = metric,
      value = colMeans(Lm) * scale,
      mcs_p = mcs_pvalues(Lm, B, block),
      n_forecasts = nrow(Lm)
    )
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
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
