# HAR family, appendix F.2: HAR (33), vector HAR (34), vector HAR with one
# common factor VHARF (35) and the log versions LHAR / VLHAR (F.2.4).
# all are direct h-step regressions estimated by OLS on the rolling window.

# HAR regressors of one series: today, mean of last 5 days, mean of last 22 days.
# value at row t only uses rows t-21 .. t so there is no look-ahead.
# na_avg (trading-calendar mode, review B1): rows are trading days and NA means
# not observed. the weekly / monthly terms are then the mean of the observed days
# among the last 5 / 22 trading days; the daily term stays NA on a missing day.
# without NA both versions are the same numbers.
har_features <- function(v, na_avg = FALSE) {
  if (na_avg && anyNA(v)) {
    return(cbind(day = v, week = rolling_mean_observed(v, 5), month = rolling_mean_observed(v, 22)))
  }
  week <- stats::filter(v, rep(1 / 5, 5), sides = 1)
  month <- stats::filter(v, rep(1 / 22, 22), sides = 1)
  cbind(day = v, week = as.numeric(week), month = as.numeric(month))
}

# mean of the non-missing values among v[t - k + 1], ..., v[t] (NA for t < k or
# when all k values are missing)
rolling_mean_observed <- function(v, k) {
  ok <- !is.na(v)
  cs <- cumsum(ifelse(ok, v, 0))
  cn <- cumsum(ok)
  lag_k <- function(x) c(rep(0, k), x[seq_len(length(x) - k)])
  n <- cn - lag_k(cn)
  out <- (cs - lag_k(cs)) / n
  out[seq_len(min(k - 1, length(v)))] <- NA
  out[n == 0] <- NA
  out
}

# design matrix for one model.
#   V    : matrix of RV (or log RV), column 1 is the series we forecast
#   type : "har" own regressors only
#          "vhar" own + every other series (34)
#          "vharf" own + the average of the other series (35)
# the common factor of a day is defined only if all other series are observed
har_design <- function(V, type = c("har", "vhar", "vharf"), na_avg = FALSE) {
  type <- match.arg(type)
  V <- as.matrix(V)
  own <- har_features(V[, 1], na_avg)
  if (type == "har" || ncol(V) == 1) return(own)
  others <- V[, -1, drop = FALSE]
  if (type == "vhar") {
    extra <- do.call(cbind, lapply(seq_len(ncol(others)), function(j) har_features(others[, j], na_avg)))
  } else {
    # common factor: features of the other series averaged with equal weights
    extra <- har_features(rowMeans(others), na_avg)
  }
  cbind(own, extra)
}

# forecast of RV of column 1 for every h, using only the rows we are given.
#   L      : log vol up to the forecast origin (last row = today), columns = assets
#   window : how many last days go into the regression
#   type   : "har", "vhar" or "vharf"
#   log    : FALSE regresses RV, TRUE regresses log RV and returns exp(fit + s2 / 2)
#   na_avg : trading-calendar mode, see har_features()
forecast_har_window <- function(L, horizons, window, type = "har", log = FALSE, na_avg = FALSE) {
  har_fit_window(L, horizons, window, type, log, na_avg)$forecast
}

# the same fit with diagnostics, one row per h:
#   pred_raw : OLS prediction before the positivity filter (log scale for LHAR)
#   filtered : TRUE when a linear HAR predicted RV <= 0 and the filter was used
#   n_train, n_coef, rank : regression size and numerical rank
har_fit_window <- function(L, horizons, window, type = "har", log = FALSE, na_avg = FALSE) {
  L <- as.matrix(L)
  V <- if (log) L else exp(L)
  X <- cbind(1, har_design(V, type, na_avg))
  y <- V[, 1]
  t <- nrow(V)
  first <- t - window + 1
  if (anyNA(X[t, ])) stop("HAR regressors are not observed at the origin")

  res <- vapply(horizons, function(h) {
    # training pairs (row s, target s + h) with s + h <= t. so the newest target is today.
    # a pair needs all regressors at s and the response at s + h
    s <- first:(t - h)
    s <- s[stats::complete.cases(X[s, , drop = FALSE]) & !is.na(y[s + h])]
    fit <- stats::lm.fit(X[s, , drop = FALSE], y[s + h])
    b <- fit$coefficients
    b[is.na(b)] <- 0 # in case of a collinear column
    pred <- sum(X[t, ] * b)
    filtered <- FALSE
    if (!log) {
      # a linear HAR can forecast RV <= 0 (rare, VHAR at long h). then QLIKE is not
      # defined so we use the smallest RV of the window instead ("insanity filter")
      filtered <- !(pred > 0)
      fc <- if (filtered) min(y[first:t], na.rm = TRUE) else pred
    } else {
      s2 <- sum(fit$residuals^2) / (length(s) - fit$rank)
      fc <- exp(pred + s2 / 2)
    }
    c(fc, pred, filtered, length(s), fit$rank)
  }, numeric(5))
  data.frame(h = horizons, forecast = res[1, ], pred_raw = res[2, ], filtered = res[3, ] == 1,
             n_train = res[4, ], n_coef = ncol(X), rank = res[5, ])
}
