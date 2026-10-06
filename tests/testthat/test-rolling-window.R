# tests for the rolling-window pipeline and the HAR models

# small synthetic panel of log vol: 3 correlated random walks around log(0.2)
fake_panel <- function(n = 160, seed = 4) {
  set.seed(seed)
  e <- matrix(rnorm(n * 3, 0, 0.1), n, 3) %*% chol(matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3))
  L <- log(0.2) + apply(e, 2, cumsum) * 0.3
  data.frame(date = as.Date("2010-01-01") + seq_len(n), A = L[, 1], B = L[, 2], C = L[, 3])
}

test_that("HAR features are the daily value and the 5 and 22 day means of the past", {
  v <- (1:30)^1.5
  f <- har_features(v)
  expect_equal(unname(f[10, "day"]), v[10])
  expect_equal(unname(f[10, "week"]), mean(v[6:10]))
  expect_equal(unname(f[25, "month"]), mean(v[4:25]))
  expect_true(all(is.na(f[1:21, "month"])))
})

test_that("VHARF with two assets is the same model as VHAR with two assets", {
  # the paper reports identical numbers for VHAR2 and VHARF2 (tables 3 and 15)
  L <- as.matrix(fake_panel()[, c("A", "B")])
  expect_equal(forecast_har_window(L, c(1, 5), 100, "vhar"),
               forecast_har_window(L, c(1, 5), 100, "vharf"))
})

test_that("HAR forecast of a constant series is that constant", {
  L <- matrix(log(0.3), 80, 1)
  expect_equal(forecast_har_window(L, c(1, 3), 50, "har"), c(0.3, 0.3))
})

test_that("forecasts do not change when the future is changed (no look-ahead)", {
  p <- fake_panel()
  L <- as.matrix(p[, c("A", "B", "C")])
  L2 <- L
  t <- 120
  L2[(t + 1):nrow(L), ] <- L2[(t + 1):nrow(L), ] + 5 # change everything after the origin
  for (m in code_b_models(c("A", "B", "C"), window = 100, delta = 1 / 252)) {
    f1 <- rolling_forecasts(L, p$date, m, t, c(1, 5))
    f2 <- rolling_forecasts(L2, p$date, m, t, c(1, 5))
    expect_equal(f1$forecast, f2$forecast, info = m$name)
  }
})

test_that("all models share the same test dates and origins come before targets", {
  p <- fake_panel()
  cfg <- list(window = 100L, origin_step = 7L, horizons = c(1L, 5L), delta = 1 / 252)
  # run_all_models writes a log under results/, keep it out of the project
  old <- setwd(tempdir())
  on.exit(setwd(old))
  fc <- run_all_models(p, c("A", "B", "C"), cfg)
  # k = 1: fBm, HAR, LHAR. k = 2, 3: one model of each of the 4 classes
  expect_equal(length(unique(fc$model)), 3 + 2 * 4)
  chk <- check_forecasts(fc)
  expect_true(chk$same_test_dates)
  expect_true(chk$no_lookahead)
  expect_true(chk$no_na_inf)
  expect_true(chk$positive_forecasts)
  # the realised value is exp(log vol) of the first asset at the target date
  x <- fc[fc$model == "fBm" & fc$h == 5, ][1, ]
  expect_equal(x$actual, exp(p$A[p$date == x$target_date]))
})

test_that("forecast origins start after a full window and leave room for h = 1", {
  o <- forecast_origins(1000, 500)
  expect_equal(min(o), 500)
  expect_equal(max(o), 999)
  expect_equal(length(forecast_origins(1000, 500, 10)), 50)
})

test_that("check_forecasts stops on look-ahead", {
  bad <- data.frame(model = "x", h = 1, origin_date = as.Date("2020-01-02"),
                    target_date = as.Date("2020-01-02"), forecast = 1, actual = 1)
  expect_error(check_forecasts(bad), "failed")
})
