# Risk Lab data: download, parsing and preparation of log-volatility panels.
#
# Source: Dacheng Xiu's Risk Lab (https://dachxiu.chicagobooth.edu/#risklab),
# endpoint data.php?ticker=<PERMNO>. Each file has a 6-line header
# (ticker, permno, company, number of days, first day, last day) followed by
# one line per day with 12 space-separated fields:
#   permno date qmle_trade ma_trade ci_trade rv5_trade rv15_trade
#               qmle_quote ma_quote ci_quote rv5_quote rv15_quote
# Volatilities are annualized (square-root) daily estimates.

risklab_columns <- c(
  "permno", "date", "qmle_trade", "ma_trade", "ci_trade", "rv5_trade", "rv15_trade",
  "qmle_quote", "ma_quote", "ci_quote", "rv5_quote", "rv15_quote"
)

risklab_file <- function(ticker, dir = file.path("data", "raw")) {
  file.path(dir, sprintf("risklab_%s.txt", ticker))
}

# Download one PERMNO (skips existing files unless overwrite = TRUE).
download_risklab <- function(ticker, permno, base_url, dir = file.path("data", "raw"),
                             overwrite = FALSE, pause = 1) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dest <- risklab_file(ticker, dir)
  if (file.exists(dest) && !overwrite) return(invisible(dest))
  url <- paste0(base_url, permno)
  tmp <- tempfile(fileext = ".txt")
  utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
  lines <- readLines(tmp, warn = FALSE)
  if (length(lines) < 10) stop(sprintf("Risk Lab returned no data for %s (%s)", ticker, permno))
  if (trimws(lines[2]) != as.character(permno)) {
    stop(sprintf("PERMNO mismatch for %s: expected %s, got %s", ticker, permno, lines[2]))
  }
  file.copy(tmp, dest, overwrite = TRUE)
  Sys.sleep(pause)
  invisible(dest)
}

# Parse a Risk Lab file into a data.frame (one row per day).
read_risklab <- function(path, ticker = NULL) {
  lines <- readLines(path, warn = FALSE)
  header <- trimws(lines[1:6])
  body <- lines[-(1:6)]
  body <- body[nzchar(trimws(body))]
  df <- utils::read.table(text = body, col.names = risklab_columns,
                          colClasses = c("character", "character", rep("numeric", 10)))
  df$date <- as.Date(df$date, format = "%Y%m%d")
  df$ticker <- if (is.null(ticker)) header[1] else ticker
  df$risklab_ticker <- header[1]
  df$company <- header[3]
  df
}

# Keep valid observations of the chosen measure inside [start, end].
# Mirrors the Risk Lab viewer, which drops values below 1e-6.
clean_risklab <- function(df, measure, start, end) {
  v <- df[[measure]]
  keep <- !is.na(v) & is.finite(v) & v > 1e-6 & df$date >= start & df$date <= end
  out <- df[keep, c("ticker", "date")]
  out$vol <- v[keep]
  out$logvol <- log(v[keep])
  out <- out[order(out$date), ]
  out[!duplicated(out$date), ]
}

# Wide panel of log volatility: one row per date, one column per ticker (NA if missing).
build_logvol_panel <- function(tickers, measure, start, end, dir = file.path("data", "raw")) {
  series <- lapply(tickers, function(tk) {
    clean_risklab(read_risklab(risklab_file(tk, dir), tk), measure, start, end)
  })
  dates <- sort(unique(do.call(c, lapply(series, `[[`, "date"))))
  panel <- data.frame(date = dates)
  for (k in seq_along(tickers)) {
    panel[[tickers[k]]] <- series[[k]]$logvol[match(dates, series[[k]]$date)]
  }
  panel
}

# Increments of log volatility for a set of columns, using only dates where
# all of them are observed (consecutive available observations).
panel_increments <- function(panel, cols) {
  sub <- panel[stats::complete.cases(panel[, cols, drop = FALSE]), c("date", cols), drop = FALSE]
  X <- apply(as.matrix(sub[, cols, drop = FALSE]), 2, diff)
  if (!is.matrix(X)) X <- matrix(X, ncol = length(cols))
  colnames(X) <- cols
  list(X = X, dates = sub$date[-1])
}

# Summary of each series: coverage, sample length and basic moments.
summarise_panel <- function(panel, sample_name) {
  tickers <- setdiff(names(panel), "date")
  do.call(rbind, lapply(tickers, function(tk) {
    v <- panel[[tk]]
    ok <- !is.na(v)
    data.frame(
      sample = sample_name, ticker = tk,
      first_date = min(panel$date[ok]), last_date = max(panel$date[ok]),
      n_obs = sum(ok), share_of_panel_dates = mean(ok),
      mean_vol = mean(exp(v[ok])), mean_logvol = mean(v[ok]), sd_logvol = stats::sd(v[ok]),
      sd_dlogvol = stats::sd(diff(v[ok]))
    )
  }))
}

# Load a processed panel written by scripts/02_prepare_data.R.
load_logvol_panel <- function(name, dir = file.path("data", "processed")) {
  p <- utils::read.csv(file.path(dir, sprintf("%s_logvol.csv", name)), check.names = FALSE)
  p$date <- as.Date(p$date)
  p
}

# Estimates on the full panel:
#   H, sigma2: each series on its own available dates;
#   rho, eta : each pair on the dates where both series are observed.
estimate_panel <- function(panel, delta, R = 20000L) {
  tickers <- setdiff(names(panel), "date")
  d <- length(tickers)
  H <- s2 <- seH <- seS2 <- nobs <- numeric(d)
  for (k in seq_len(d)) {
    x <- panel_increments(panel, tickers[k])$X[, 1]
    H[k] <- hurst_mm(x)
    s2[k] <- sigma2_mm(x, delta, H[k])
    nobs[k] <- length(x)
    seH[k] <- sqrt(avar_hurst(H[k], R) / nobs[k])
    seS2[k] <- se_sigma2(H[k], s2[k], nobs[k], delta, R)
  }
  rho <- eta <- se_rho <- se_eta <- npair <- matrix(NA_real_, d, d, dimnames = list(tickers, tickers))
  diag(rho) <- 1
  diag(eta) <- 0
  for (i in seq_len(d)) for (j in seq_len(d)) if (i < j) {
    X <- panel_increments(panel, tickers[c(i, j)])$X
    hi <- hurst_mm(X[, 1]); hj <- hurst_mm(X[, 2])
    r <- rho_mm(X[, 1], X[, 2])
    e <- eta_mm(X[, 1], X[, 2])
    nn <- nrow(X)
    rho[i, j] <- rho[j, i] <- r
    eta[i, j] <- e
    eta[j, i] <- -e
    se_rho[i, j] <- se_rho[j, i] <- sqrt(avar_rho(hi, hj, r, R) / nn)
    se_eta[i, j] <- se_eta[j, i] <- sqrt(avar_eta(hi, hj, r, R) / nn)
    npair[i, j] <- npair[j, i] <- nn
  }
  names(H) <- names(s2) <- names(seH) <- names(seS2) <- names(nobs) <- tickers
  list(H = H, sigma2 = s2, se_H = seH, se_sigma2 = seS2, n = nobs,
       rho = rho, eta = eta, se_rho = se_rho, se_eta = se_eta, n_pair = npair)
}

# Rolling-window Hurst estimates (Figures 7-8): window of `window` observations,
# evaluated every `step` days; the estimate is dated at the window's last day.
rolling_hurst <- function(panel, tickers, window, step) {
  out <- list()
  for (tk in tickers) {
    v <- panel[[tk]]
    ok <- which(!is.na(v))
    lv <- v[ok]
    dts <- panel$date[ok]
    ends <- seq(window + 1, length(lv), by = step)
    h <- vapply(ends, function(e) hurst_mm(diff(lv[(e - window):e])), 0)
    out[[tk]] <- data.frame(ticker = tk, date = dts[ends], H = h)
  }
  do.call(rbind, out)
}
