# HAR family, appendix F.2: HAR (33), vector HAR (34), vector HAR with one
# common factor VHARF (35) and the log versions LHAR / VLHAR (F.2.4).
# all are direct h-step regressions estimated by OLS on the rolling window.

# HAR regressors of one series: today, mean of last 5 days, mean of last 22 days.
# value at row t only uses rows t-21 .. t so there is no look-ahead
har_features <- function(v) {
  week <- stats::filter(v, rep(1 / 5, 5), sides = 1)
  month <- stats::filter(v, rep(1 / 22, 22), sides = 1)
  cbind(day = v, week = as.numeric(week), month = as.numeric(month))
}

# design matrix for one model.
#   V    : matrix of RV (or log RV), column 1 is the series we forecast
#   type : "har" own regressors only
#          "vhar" own + every other series (34)
#          "vharf" own + the average of the other series (35)
har_design <- function(V, type = c("har", "vhar", "vharf")) {
  type <- match.arg(type)
  V <- as.matrix(V)
  own <- har_features(V[, 1])
  if (type == "har" || ncol(V) == 1) return(own)
  others <- V[, -1, drop = FALSE]
  if (type == "vhar") {
    extra <- do.call(cbind, lapply(seq_len(ncol(others)), function(j) har_features(others[, j])))
  } else {
    # common factor: features of the other series averaged with equal weights
    extra <- har_features(rowMeans(others))
  }
  cbind(own, extra)
}

# forecast of RV of column 1 for every h, using only the rows we are given.
#   L      : log vol up to the forecast origin (last row = today), columns = assets
#   window : how many last days go into the regression
#   type   : "har", "vhar" or "vharf"
#   log    : FALSE regresses RV, TRUE regresses log RV and returns exp(fit + s2 / 2)
forecast_har_window <- function(L, horizons, window, type = "har", log = FALSE) {
  L <- as.matrix(L)
  V <- if (log) L else exp(L)
  X <- cbind(1, har_design(V, type))
  y <- V[, 1]
  t <- nrow(V)
  first <- t - window + 1

  sapply(horizons, function(h) {
    # training pairs (row s, target s + h) with s + h <= t. so the newest target is today
    s <- first:(t - h)
    s <- s[stats::complete.cases(X[s, , drop = FALSE])]
    fit <- stats::lm.fit(X[s, , drop = FALSE], y[s + h])
    b <- fit$coefficients
    b[is.na(b)] <- 0 # in case of a collinear column
    pred <- sum(X[t, ] * b)
    # a linear HAR can forecast RV <= 0 (rare, VHAR at long h). then QLIKE is not
    # defined so we use the smallest RV of the window instead ("insanity filter")
    if (!log) return(if (pred > 0) pred else min(y[first:t]))
    s2 <- sum(fit$residuals^2) / (length(s) - fit$rank)
    exp(pred + s2 / 2)
  })
}
