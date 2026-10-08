# review B1: origins, training windows and h-day targets on the trading calendar.
# every test goes through run_all_models(), not through a pre-sliced matrix.

# synthetic 5-dimensional mfBm log-vol panel (as helper-code-c.R::c_panel(), kept here
# so that the Code B tests do not need Code C) on an explicit Monday-Friday calendar
cal_panel <- function(n = 105L) {
  set.seed(802)
  H <- c(.12, .2, .28, .35, .42, .17, .31)
  r <- matrix(.25, 7, 7); diag(r) <- 1
  x <- simulate_mfbm(n - 1L, H, r, delta = 1 / 252)[[1]]
  L <- log(.2) + increments_to_levels(x) * .18
  colnames(L) <- LETTERS[1:7]
  d <- seq(as.Date("2010-01-04"), by = "day", length.out = 2L * n)
  data.frame(date = utils::head(d[!format(d, "%u") %in% c("6", "7")], n), L[, 1:5], check.names = FALSE)
}
cal_cfg <- function(calendar, t, p) {
  list(window = 60L, delta = 1 / 252, horizons = c(1L, 5L), origin_step = 1L,
       origin_dates = p$date[t], calendar = calendar)
}
run_quiet <- function(...) {
  old <- setwd(tempdir()); on.exit(setwd(old)) # run_all_models logs under results/
  run_all_models(...)
}
key_cols <- c("model", "h", "origin_date", "target_date")

test_that("trading calendar: future values and NA patterns of other assets change nothing", {
  p <- cal_panel(); t <- 70L
  cfg <- cal_cfg("trading", t, p)
  base <- run_quiet(p, LETTERS[1:5], cfg)
  expect_equal(length(unique(base$model)), 1 + 2 + 4 * 4) # fBm, HAR, LHAR + 4 classes x 4
  mutations <- list(
    within(p, B[t + 2L] <- NA),
    within(p, { C[(t + 1L):(t + 5L)] <- NA; E[t + 1L] <- NA }),
    within(p, { B[(t + 1L):nrow(p)] <- 3; D[(t + 1L):nrow(p)] <- -9 }))
  for (q in mutations) {
    fc <- run_quiet(q, LETTERS[1:5], cfg)
    expect_equal(fc[, key_cols], base[, key_cols])
    expect_equal(fc$forecast, base$forecast)
    expect_equal(fc$actual, base$actual)
  }
  # the targets are fixed in advance: h trading days after the origin
  expect_true(all(match(base$target_date, p$date) - t == base$h))
})

test_that("trading calendar: a missing future target is dropped for all models, never moved", {
  p <- cal_panel(); t <- 70L
  cfg <- cal_cfg("trading", t, p)
  base <- run_quiet(p, LETTERS[1:5], cfg)
  q <- p; q$A[t + 5L] <- NA
  fc <- run_quiet(q, LETTERS[1:5], cfg)
  expect_equal(unique(fc$h), 1L)
  b1 <- base[base$h == 1L, c(key_cols, "forecast", "actual")]
  rownames(b1) <- NULL
  expect_equal(fc[, c(key_cols, "forecast", "actual")], b1)
  expect_false(any(fc$target_date == p$date[t + 6L]))
})

test_that("common_obs mode keeps the historical convention, where the target moves", {
  # labelled replication mode: rows are common observation days, so a future NA of
  # another asset turns the 5th common day after the origin into the 6th trading day
  p <- cal_panel(); t <- 70L
  cfg <- cal_cfg("common_obs", t, p)
  q <- p; q$B[t + 2L] <- NA
  a <- run_quiet(p, LETTERS[1:5], cfg, model_names = "fBm")
  b <- run_quiet(q, LETTERS[1:5], cfg, model_names = "fBm")
  expect_equal(a$forecast, b$forecast)
  expect_equal(a$target_date[a$h == 5L], p$date[t + 5L])
  expect_equal(b$target_date[b$h == 5L], p$date[t + 6L])
})

test_that("trading calendar: an earlier missing training day keeps the origin eligible", {
  # eligibility rule: origin day observed for all assets and >= 90% of the window
  # observed. the missing day is dropped from the models that use that asset only;
  # the remaining days keep their true calendar times (no compression, no filling)
  p <- cal_panel(); t <- 70L
  cfg <- cal_cfg("trading", t, p)
  q <- p; q$B[t - 10L] <- NA
  base <- run_quiet(p, LETTERS[1:5], cfg)
  fc <- run_quiet(q, LETTERS[1:5], cfg)
  expect_equal(fc[, key_cols], base[, key_cols])
  same <- fc$model %in% c("fBm", "HAR", "LHAR")
  expect_equal(fc$forecast[same], base$forecast[same]) # B is not used by these models
  expect_true(all(fc$n_train_days[same] == 60L))
  expect_true(all(fc$n_train_days[!same] == 59L))
  expect_false(isTRUE(all.equal(fc$forecast[fc$model == "bfBm"], base$forecast[base$model == "bfBm"])))
  # bfBm is the calendar forecaster on the last 60 trading days, B[t - 10] missing
  Y <- as.matrix(q[(t - 59L):t, c("A", "B")])
  expect_equal(fc$forecast[fc$model == "bfBm"], forecast_mfbm_calendar(Y, c(1, 5), 1 / 252))
  # which differs from compressing the calendar (the common_obs convention)
  expect_false(isTRUE(all.equal(fc$forecast[fc$model == "bfBm"],
    forecast_mfbm_window(as.matrix(q[stats::complete.cases(q[, c("A", "B")]), c("A", "B")][(t - 60L):(t - 1L), ]),
                         c(1, 5), 1 / 252))))
})

test_that("trading calendar: a missing asset at the origin skips the origin for every model", {
  p <- cal_panel(); t <- 70L
  q <- p; q$C[t] <- NA
  expect_error(run_quiet(q, LETTERS[1:5], cal_cfg("trading", t, q)), "not eligible")
  cfg <- list(window = 60L, delta = 1 / 252, horizons = c(1L, 5L), origin_step = 5L, calendar = "trading")
  fc <- run_quiet(q, LETTERS[1:5], cfg)
  expect_false(any(fc$origin_date == q$date[t]))
  expect_true(all(table(fc$model) == table(fc$model)[1])) # same keys for every model
  expect_true(check_forecasts(fc)$same_test_dates)
  # a window with less than 90% observed days is skipped too
  q2 <- p; q2$D[(t - 20L):(t - 10L)] <- NA
  el <- code_b_eligible_origins(as.matrix(q2[, LETTERS[1:5]]), t, 60L)
  expect_false(el$eligible)
  expect_equal(el$reason, "low_window_coverage")
})

test_that("without missing days both conventions give the same forecasts", {
  p <- cal_panel()
  cfg <- list(window = 60L, delta = 1 / 252, horizons = c(1L, 5L), origin_step = 9L)
  a <- run_quiet(p, LETTERS[1:5], modifyList(cfg, list(calendar = "common_obs")))
  b <- run_quiet(p, LETTERS[1:5], modifyList(cfg, list(calendar = "trading")))
  cols <- c(key_cols, "forecast", "actual", "train_start", "n_train_days")
  expect_equal(a[, cols], b[, cols])
  expect_equal(unique(a$calendar), "common_obs")
  expect_equal(unique(b$calendar), "trading")
})

test_that("gap-aware estimators equal the code A estimators when nothing is missing", {
  p <- cal_panel()
  Y <- as.matrix(p[1:80, LETTERS[1:4]])
  a <- window_mfbm_params(apply(Y, 2, diff), 1 / 252)
  b <- window_mfbm_params_calendar(Y, 1 / 252)
  expect_equal(unname(b$H), unname(a$H))
  expect_equal(unname(b$sigma), unname(a$sigma))
  expect_equal(unname(b$rho), unname(a$rho))
  # an increment that would span a missing day is not used
  Y2 <- Y; Y2[40, 2] <- NA
  d1 <- diff(Y2); d1 <- d1[stats::complete.cases(d1), ]
  expect_equal(nrow(d1), 77L)
  expect_equal(unname(window_mfbm_params_calendar(Y2, 1 / 252)$rho[1, 2]),
               sum(d1[, 1] * d1[, 2]) / sqrt(sum(d1[, 1]^2) * sum(d1[, 2]^2)))
})

test_that("observed-day HAR averages equal the usual ones without NA", {
  v <- (1:40)^1.2
  expect_equal(har_features(v, na_avg = TRUE), har_features(v))
  expect_equal(rolling_mean_observed(v, 5)[10], mean(v[6:10]))
  w <- v; w[8] <- NA
  f <- har_features(w, na_avg = TRUE)
  expect_true(is.na(f[8, "day"]))
  expect_equal(unname(f[10, "week"]), mean(v[c(6, 7, 9, 10)]))
  expect_equal(unname(f[30, "month"]), mean(v[c(9:30)]))
  expect_equal(unname(f[25, "month"]), mean(w[4:25], na.rm = TRUE))
})

test_that("checkpoint fingerprints distinguish the calendar and the coverage rule", {
  p <- cal_panel(); cfg <- list(window = 60L, delta = 1 / 252, horizons = 1L, origin_step = 9L)
  k0 <- forecast_fingerprint(p, LETTERS[1:5], cfg)
  k1 <- forecast_fingerprint(p, LETTERS[1:5], modifyList(cfg, list(calendar = "common_obs")))
  k2 <- forecast_fingerprint(p, LETTERS[1:5], modifyList(cfg, list(calendar = "trading")))
  k3 <- forecast_fingerprint(p, LETTERS[1:5], modifyList(cfg, list(calendar = "trading", min_window_coverage = .8)))
  expect_equal(length(unique(c(k1, k2, k3))), 3L)
  expect_false(identical(k0, k2))
  expect_error(code_b_calendar(list(calendar = "weekly")))
})

# --- second review: origins without an evaluable horizon ---------------------------

no_eval <- "code_b_no_evaluable_forecasts"

test_that("one origin whose targets are all missing gives a typed empty table and a clear error", {
  p <- cal_panel(); t <- 70L
  q <- p; q$A[c(t + 1L, t + 5L)] <- NA
  L <- as.matrix(q[, LETTERS[1:5]])
  m <- code_b_models(LETTERS[1:5], 60L, 1 / 252, "trading")
  e <- rolling_forecasts(L, q$date, m$mfBm5, t, c(1L, 5L))
  ok <- rolling_forecasts(as.matrix(p[, LETTERS[1:5]]), p$date, m$mfBm5, t, c(1L, 5L))
  expect_equal(nrow(e), 0L)
  expect_identical(lapply(e, class), lapply(ok, class)) # same schema, Date columns kept
  expect_s3_class(e$target_date, "Date")
  expect_error(run_quiet(q, LETTERS[1:5], cal_cfg("trading", t, q)), class = no_eval)
  expect_error(run_quiet(q, LETTERS[1:5], cal_cfg("trading", t, q)), "no evaluable forecasts")
})

test_that("one origin with a single horizon whose target is missing", {
  p <- cal_panel(); t <- 70L
  q <- p; q$A[t + 1L] <- NA
  cfg <- modifyList(cal_cfg("trading", t, q), list(horizons = 1L))
  expect_error(run_quiet(q, LETTERS[1:5], cfg), class = no_eval)
  # the same origin still forecasts h = 5 when it is configured; the missing h = 1 is not moved
  fc <- run_quiet(q, LETTERS[1:5], modifyList(cfg, list(horizons = c(1L, 5L))), model_names = c("fBm", "VHAR5"))
  expect_equal(unique(fc$h), 5L)
  expect_true(all(fc$target_date == q$date[t + 5L]))
})

test_that("automatic origins with only long horizons skip trailing origins in both calendars", {
  p <- cal_panel()
  for (cal in code_b_calendars()) {
    cfg <- list(window = 60L, delta = 1 / 252, horizons = c(15L, 20L), origin_step = 1L, calendar = cal)
    fc <- run_quiet(p, LETTERS[1:5], cfg, model_names = c("fBm", "bfBm", "HAR", "VHAR3"))
    expect_true(check_forecasts(fc)$same_test_dates)
    expect_equal(max(fc$origin_date), p$date[nrow(p) - 15L]) # later origins have no target in the panel
    expect_true(all(match(fc$target_date, p$date) - match(fc$origin_date, p$date) == fc$h))
    # identical to running only the origins that have an evaluable horizon
    valid <- p$date[60:(nrow(p) - 15L)]
    ref <- run_quiet(p, LETTERS[1:5], modifyList(cfg, list(origin_dates = valid)),
                     model_names = c("fBm", "bfBm", "HAR", "VHAR3"))
    expect_equal(fc, ref)
    # the checkpoint path stores and reloads the same table
    old <- setwd(tempdir()); on.exit(setwd(old), add = TRUE)
    unlink(code_b_results("forecasts", "checkpoints"), recursive = TRUE)
    a <- run_all_models(p, LETTERS[1:5], cfg, model_names = c("fBm", "VHAR3"), checkpoint = "rt")
    b <- run_all_models(p, LETTERS[1:5], cfg, model_names = c("fBm", "VHAR3"), checkpoint = "rt")
    setwd(old)
    expect_equal(a, b)
    sub <- fc[fc$model %in% c("fBm", "VHAR3"), ]; rownames(sub) <- NULL
    expect_equal(a, sub)
  }
})

test_that("valid origins are unchanged and empty origins are omitted for every model", {
  p <- cal_panel(); t <- 70L; o1 <- t - 10L
  q <- p; q$A[c(t + 1L, t + 5L)] <- NA # targets of origin t only; origin o1 trains on rows <= o1
  cfg <- modifyList(cal_cfg("trading", t, q), list(origin_dates = q$date[c(o1, t)]))
  fc <- run_quiet(q, LETTERS[1:5], cfg)
  base <- run_quiet(p, LETTERS[1:5], modifyList(cfg, list(origin_dates = p$date[o1])))
  expect_equal(unique(fc$origin_date), p$date[o1])
  expect_equal(fc, base)
  expect_true(all(table(fc$model) == 2L)) # h = 1 and 5 for each of the 19 models
})

test_that("a run with no evaluable rows raises the domain error in both calendars", {
  p <- cal_panel(); n <- nrow(p)
  for (cal in code_b_calendars()) {
    cfg <- list(window = 60L, delta = 1 / 252, horizons = c(5L, 10L), origin_step = 1L,
                origin_dates = p$date[c(n - 4L, n - 1L)], calendar = cal)
    expect_error(run_quiet(p, LETTERS[1:5], cfg), class = no_eval)
    err <- tryCatch(run_quiet(p, LETTERS[1:5], cfg), error = identity)
    expect_match(conditionMessage(err), paste0("none of the 2 origins.*horizons 5, 10 .*calendar = ", cal))
  }
})

test_that("rolling_mean_observed returns NA when the series is shorter than k", {
  expect_equal(rolling_mean_observed(c(1, NA, 3), 5), rep(NA_real_, 3))
  expect_equal(rolling_mean_observed(numeric(), 5), numeric())
  expect_equal(rolling_mean_observed(c(1, NA, 3, 5, 7, 9), 5), c(NA, NA, NA, NA, 4, 6))
})
