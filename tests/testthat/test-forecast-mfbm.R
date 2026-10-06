# tests for the optimal mfBm forecasts (section 4, appendix B and E.3)

n <- 500
times <- (1:n) / 250
targets <- (n + 1:5) / 250

test_that("theoretical RMSFEs equal the bracketed values of table 8", {
  uni <- mfbm_forecast_weights(mfbm_params(0.1, matrix(1)), times, targets)
  expect_equal(round(sqrt(uni$msfe), 4), c(.4802, .5077, .5254, .5387, .5495))
  biv <- mfbm_forecast_weights(mfbm_params(c(0.1, 0.4), 0.8), times, targets, target = 1)
  expect_equal(round(sqrt(biv$msfe), 4), c(.4246, .4526, .4700, .4827, .4927))
  biv2 <- mfbm_forecast_weights(mfbm_params(c(0.1, 0.4), 0.8), times, targets, target = 2)
  expect_equal(round(sqrt(biv2$msfe), 4), c(.0953, .1242, .1443, .1602, .1734))
})

test_that("general one-observation forecast equals the closed forms of prop. 4.1 and 4.3", {
  for (H2 in c(0.15, 0.4, 0.7)) {
    f <- one_obs_forecast(mfbm_params(c(0.4, H2), 0.5, 0, c(1, 2)), t = 1.5, h = 0.7)
    expect_equal(f$w, prop41_weights(1.5, 0.7, 0.4, H2, 0.5, 1, 2))
  }
  for (d in c(1, 2, 5, 20)) {
    f <- one_obs_forecast(equicorr_params(d, 0.4, 0.1, 0.8), t = 10, h = 1)
    expect_equal(f$msfe, prop43_msfe(d, 10, 1, 0.4, 0.1, 0.8))
  }
})

test_that("no gain from the second component when H are equal or rho = 0 (prop. 4.1, B.1)", {
  # one observation
  expect_equal(prop41_weights(1, 1, 0.3, 0.3, 0.6)[2], 0)
  # full history, prop. B.1: the weights on component 2 are zero
  short <- (1:60) / 250
  W <- mfbm_forecast_weights(mfbm_params(c(0.2, 0.2), 0.6), short, 61 / 250)$W
  on_comp2 <- W[seq(2, length(W), by = 2)]
  expect_lt(max(abs(on_comp2)), 1e-8)
  W0 <- mfbm_forecast_weights(mfbm_params(c(0.1, 0.4), 0), short, 61 / 250)$W
  expect_lt(max(abs(W0[seq(2, length(W0), by = 2)])), 1e-10)
})

test_that("admissibility check agrees with rho_max of code A", {
  for (H2 in c(0.1, 0.3, 0.8)) {
    rm <- rho_max(0.1, H2)
    if (rm < 1) expect_false(is_admissible_mfbm(mfbm_params(c(0.1, H2), min(rm + 0.01, 1))))
    expect_true(is_admissible_mfbm(mfbm_params(c(0.1, H2), rm - 0.01)))
  }
})

test_that("window parameters are shrunk to an admissible mfBm when needed", {
  set.seed(5)
  X <- simulate_mfbm(400, c(0.05, 0.6), 0.3, delta = 1 / 252)[[1]]
  X[, 2] <- X[, 2] + 3 * X[, 1] # makes the increments almost perfectly correlated
  par <- window_mfbm_params(X, 1 / 252)
  expect_true(is_admissible_mfbm(par))
})

test_that("empirical window forecast is positive and starts near the last value", {
  set.seed(6)
  X <- simulate_mfbm(300, c(0.15, 0.3), 0.5, sigma = c(1.5, 1.5), delta = 1 / 252)[[1]]
  Y <- log(0.25) + increments_to_levels(X)
  f <- forecast_mfbm_window(Y, c(1, 5, 20), 1 / 252)
  expect_true(all(is.finite(f) & f > 0))
  # one day ahead the forecast of log RV stays close to today's log RV
  expect_lt(abs(log(f[1]) - Y[nrow(Y), 1]), 0.5)
})

test_that("the window forecast is the conditional mean of eq. (15) plus var / 2", {
  set.seed(8)
  X <- simulate_mfbm(120, c(0.2, 0.35), 0.4, delta = 1 / 252)[[1]]
  Y <- -1.5 + increments_to_levels(X)
  # same thing written out by hand with the stacked vector of section 4
  par <- window_mfbm_params(apply(Y, 2, diff), 1 / 252)
  n <- nrow(Y) - 1
  S <- mfbm_level_cov_matrix(par, (1:(n + 3)) / 252)
  obs <- 1:(2 * n)
  tgt <- 2 * (n + 3) - 1 # component 1 at time n + 3
  z <- as.vector(t(sweep(Y[-1, ], 2, Y[1, ])))
  m <- Y[1, 1] + sum(S[tgt, obs] * solve(S[obs, obs], z))
  v <- S[tgt, tgt] - sum(S[tgt, obs] * solve(S[obs, obs], S[obs, tgt]))
  expect_equal(forecast_mfbm_window(Y, 3, 1 / 252), exp(m + v / 2))
})

test_that("simulated RMSFE is close to the theory with known parameters", {
  set.seed(7)
  paths <- simulate_mfbm(505, c(0.1, 0.4), 0.8, delta = 1 / 250, nsim = 400)
  r <- sim_rmsfe(paths, mfbm_params(c(0.1, 0.4), 0.8), c(1, 2), 1, 500, 1:5, 1 / 250)
  # monte carlo error of an RMSFE with 400 reps is about rmsfe / sqrt(800) ~ 3.5%
  expect_lt(max(abs(r$rmsfe / r$theory - 1)), 0.12)
})
