# optimal forecasts of a time-reversible mfBm, section 4 and appendix B.
# the optimal forecast is the conditional mean of a gaussian vector, eq. (15):
#   E[B_target | X] = gamma' Sigma^-1 X
# all covariances come from code A (mfbm_level_cov_matrix, eq. (5)).

# weights of the optimal forecast of one component at several future times.
#   par          : mfbm parameters from mfbm_params() (code A)
#   obs_times    : times where all d components are observed
#   target_times : times we forecast (one column of W per time)
#   target       : which component we forecast
# returns W (weights for the stacked observations, see below) and the msfe
# of every target. observations are stacked like in section 4:
# (B1_t1, ..., Bd_t1, B1_t2, ..., Bd_t2, ...)
mfbm_forecast_weights <- function(par, obs_times, target_times, target = 1) {
  d <- par$d
  S <- mfbm_level_cov_matrix(par, c(obs_times, target_times))
  n_obs <- length(obs_times) * d
  obs <- seq_len(n_obs)
  tgt <- n_obs + (seq_along(target_times) - 1) * d + target

  S_xx <- S[obs, obs, drop = FALSE]
  S_xy <- S[obs, tgt, drop = FALSE] # this is gamma in (15)
  S_yy <- S[tgt, tgt, drop = FALSE]

  # W = Sigma^-1 gamma via cholesky
  R <- chol(S_xx)
  W <- backsolve(R, forwardsolve(t(R), S_xy))
  msfe <- diag(S_yy) - colSums(S_xy * W)
  list(W = W, msfe = msfe)
}

# --- single observation, section 4.1 ------------------------------------------

# weights and msfe of the forecast of B1_{t+h} from (B1_t, ..., Bd_t), prop. 4.2
one_obs_forecast <- function(par, t, h) {
  fw <- mfbm_forecast_weights(par, obs_times = t, target_times = t + h, target = 1)
  list(w = drop(fw$W), msfe = fw$msfe)
}

# closed form of prop. 4.1 (bivariate). used in tests to check the general code
prop41_weights <- function(t, h, H1, H2, rho, sigma1 = 1, sigma2 = 1) {
  H <- (H1 + H2) / 2
  w11 <- (w_fbm(t, h, H1) / t^(2 * H1) - rho^2 * w_fbm(t, h, H) / t^(2 * H)) / (1 - rho^2)
  w12 <- rho / (1 - rho^2) * sigma1 / sigma2 *
    (w_fbm(t, h, H) / t^(2 * H2) - w_fbm(t, h, H1) / t^(2 * H))
  c(w11, w12)
}

# closed form msfe of prop. 4.3: sigma = 1, all correlations rho,
# H1 for the first component and H for the other d - 1
prop43_msfe <- function(d, t, h, H1, H, rho) {
  Hm <- (H1 + H) / 2 # average exponent for the cross terms
  a <- w_fbm(t, h, H1)
  b <- w_fbm(t, h, Hm)
  c0 <- w_fbm(t, h, H)
  if (d == 1) return((t + h)^(2 * H1) - a^2 / t^(2 * H1))
  k <- 1 + (d - 2) * rho - (d - 1) * rho^2
  (t + h)^(2 * H1) - ((1 + (d - 2) * rho) * a^2 / t^(2 * H1) -
    2 * (d - 1) * rho^2 * a * b / t^(2 * Hm) +
    (d - 1) * rho^2 * b^2 / t^(2 * H)) / k
}

# parameters of the prop. 4.3 setting: equal correlations rho between all pairs
equicorr_params <- function(d, H1, H, rho) {
  R <- matrix(rho, d, d)
  diag(R) <- 1
  mfbm_params(c(H1, rep(H, d - 1)), R, 0, rep(1, d))
}

# --- monte carlo of appendix E.3 (tables 8-11) --------------------------------

# parameters of a sub-vector of an mfbm (the components a forecaster uses)
sub_params <- function(par, use) {
  mfbm_params(par$H[use], par$rho[use, use, drop = FALSE], 0, par$sigma[use])
}

# RMSFE of forecasts of component `target` at times (n + h) delta on simulated paths.
#   paths    : list of (n + max h) x d increment matrices from simulate_mfbm()
#   par      : true parameters of the paths
#   use      : components the forecaster sees, e.g. 1 (fBm) or c(1, 2) (bfBm)
#   estimate : FALSE uses the true parameters, TRUE estimates them on the first
#              n increments of every path (table 11 "unknown")
# returns the monte carlo RMSFE and the theoretical one (true parameters) per h
sim_rmsfe <- function(paths, par, use, target, n, horizons, delta,
                      estimate = FALSE, cores = 1L) {
  par_use <- sub_params(par, use)
  k <- which(use == target) # position of the target among the used components
  obs_times <- (1:n) * delta
  target_times <- (n + horizons) * delta
  known <- mfbm_forecast_weights(par_use, obs_times, target_times, target = k)

  one_path <- function(X) {
    B <- apply(as.matrix(X), 2, cumsum) # levels at delta, 2 delta, ... (B_0 = 0)
    z <- as.vector(t(B[1:n, use, drop = FALSE]))
    W <- known$W
    if (estimate) {
      p_hat <- window_mfbm_params(X[1:n, use, drop = FALSE], delta)
      W <- mfbm_forecast_weights(p_hat, obs_times, target_times, target = k)$W
    }
    drop(z %*% W) - B[n + horizons, target]
  }
  err <- parallel::mclapply(paths, one_path, mc.cores = if (estimate) cores else 1L)
  err <- do.call(rbind, err)
  data.frame(h = horizons, rmsfe = sqrt(colMeans(err^2)), theory = sqrt(known$msfe))
}

# --- empirical forecast on one rolling window ---------------------------------

# forecast of RV of the first column h days ahead from a window of log RV.
#   Y        : window of log vol, rows = days, columns = assets, column 1 is forecast
#   horizons : vector of h
#   delta    : sampling step (1/252)
# steps:
#   1. estimate H, sigma, rho on the increments of the window (code A, eta = 0)
#   2. treat the window as an mfBm started at its first day: Z = Y - Y[1, ]
#   3. conditional mean and variance of log RV at t + h, eq. (15)
#   4. RV is lognormal so E[RV] = exp(mean + var / 2) (footnote 8)
forecast_mfbm_window <- function(Y, horizons, delta) {
  Y <- as.matrix(Y)
  n <- nrow(Y) - 1
  par <- window_mfbm_params(apply(Y, 2, diff), delta)

  fw <- mfbm_forecast_weights(par, (1:n) * delta, (n + horizons) * delta, target = 1)

  Z <- sweep(Y[-1, , drop = FALSE], 2, Y[1, ])
  z <- as.vector(t(Z)) # stacked by time, same order as the covariance matrix
  log_mean <- Y[1, 1] + drop(z %*% fw$W)
  exp(log_mean + fw$msfe / 2)
}

# mfbm parameters estimated on a window of increments.
# H is kept inside (0.01, 0.99) because the covariance needs 0 < H < 1.
# if the estimated correlations are not admissible together with the H
# (no mfBm with these parameters exists) we shrink them by 5% steps until it is.
# par$shrink says how many steps were needed (0 almost always)
window_mfbm_params <- function(X, delta) {
  X <- as.matrix(X)
  est <- estimate_mfbm(X, delta)
  H <- pmin(pmax(est$H, 0.01), 0.99)
  sigma <- sqrt(est$sigma2)
  for (k in 0:100) {
    R <- est$rho * 0.95^k
    diag(R) <- 1
    par <- mfbm_params(H, R, 0, sigma)
    if (is_admissible_mfbm(par)) break
  }
  par$shrink <- k
  par
}

# a time-reversible mfBm exists iff the matrix
#   rho_ij Gamma(Hi + Hj + 1) sin(pi (Hi + Hj) / 2)
# is positive semidefinite (amblard et al. 2013). for d = 2 this is the
# rho_max bound of section 2. we ask for strictly positive to keep chol stable
is_admissible_mfbm <- function(par) {
  if (par$d == 1) return(TRUE)
  Hs <- outer(par$H, par$H, "+")
  M <- par$rho * gamma(Hs + 1) * sin(pi * Hs / 2)
  min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) > 1e-8
}
