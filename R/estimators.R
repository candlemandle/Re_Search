# Method-of-moments (MM) estimators of the mfBm parameters, their asymptotic
# variances, the optimal GMM correlation estimator (Appendix G) and the
# Amblard-Coeurjolly (2011) estimator used for comparison (Appendix A, Table 6).
#
# Inputs are increments: X is an n x d matrix with X[k, j] = Delta_k B^(j).

# ---------------------------------------------------------------------------
# Point estimators, eqs. (6), (7), (8a), (8b)
# ---------------------------------------------------------------------------

# lag-2 increments Delta_{k,2} B = Delta_k B + Delta_{k+1} B, k = 1..n-1
.lag2 <- function(x) x[-length(x)] + x[-1]

# eq. (6)
hurst_mm <- function(x) {
  x <- as.numeric(x)
  log(sum(.lag2(x)^2) / sum(x^2)) / (2 * log(2))
}

# eq. (7)
sigma2_mm <- function(x, delta, H = hurst_mm(x)) {
  x <- as.numeric(x)
  sum(x^2) / (length(x) * delta^(2 * H))
}

# eq. (8a)
rho_mm <- function(x1, x2) {
  sum(x1 * x2) / sqrt(sum(x1^2) * sum(x2^2))
}

# eq. (8b): asymmetry eta_{1,2}.
# Sign conventions differ inside the paper:
#   * "printed": (8b) verbatim, numerator sum(D_{k+1}B2 D_kB1 - D_{k+1}B1 D_kB2).
#     This is what the authors' BYZestimator.R computes; it is consistent with
#     their simulator simMFBM (Amblard-Coeurjolly), whose increment covariance
#     uses (rho - eta sign(h)), i.e. eta of the opposite sign to eq. (4).
#   * "proof" (default): numerator with the order of Appendix C.2, consistent for
#     eta_{1,2} as defined in Section 2 / eq. (4) / eq. (26) and used by our
#     simulator.
# eta_mm(x1, x2, "printed") == eta_mm(x2, x1, "proof"); Table 14 entry (row i,
# column j) equals eta_mm(col j, row i, "proof"). The two-sided test (11) is unaffected.
eta_mm <- function(x1, x2, convention = c("proof", "printed")) {
  convention <- match.arg(convention)
  n <- length(x1)
  num <- sum(x1[-1] * x2[-n] - x2[-1] * x1[-n])
  den <- sqrt(sum(.lag2(x1)^2) * sum(.lag2(x2)^2)) - 2 * sqrt(sum(x1^2) * sum(x2^2))
  if (convention == "printed") -num / den else num / den
}

# All MM estimates for a d-variate increment matrix.
estimate_mfbm <- function(X, delta, eta_convention = "proof") {
  X <- as.matrix(X)
  d <- ncol(X)
  H <- apply(X, 2, hurst_mm)
  s2 <- vapply(seq_len(d), function(j) sigma2_mm(X[, j], delta, H[j]), 0)
  rho <- diag(d)
  eta <- matrix(0, d, d)
  for (i in seq_len(d)) for (j in seq_len(d)) if (i < j) {
    rho[i, j] <- rho[j, i] <- rho_mm(X[, i], X[, j])
    eta[i, j] <- eta_mm(X[, i], X[, j], eta_convention)
    eta[j, i] <- -eta[i, j]
  }
  nm <- colnames(X)
  if (!is.null(nm)) {
    names(H) <- names(s2) <- nm
    dimnames(rho) <- dimnames(eta) <- list(nm, nm)
  }
  list(H = H, sigma2 = s2, rho = rho, eta = eta, n = nrow(X), delta = delta)
}

# ---------------------------------------------------------------------------
# Infinite series with power-law tail correction
# ---------------------------------------------------------------------------

# sum_{r=1}^{inf} f(r), evaluated exactly up to R and with an integral tail
# approximation sum_{r>R} t_R (r/R)^{-p} where p is the local decay exponent.
series_sum <- function(f, R = 20000L) {
  r <- seq_len(R)
  v <- f(r)
  s <- sum(v)
  tR <- v[R]
  tH <- v[R %/% 2]
  if (is.finite(tR) && tR != 0 && tH != 0 && sign(tR) == sign(tH)) {
    p <- log(tH / tR) / log(R / (R %/% 2))
    if (p > 1.0001) s <- s + tR * (R / (p - 1) - 0.5)
  }
  s
}

.g2 <- function(r, a) .apow(r + 1, a) + .apow(r - 1, a) - 2 * .apow(r, a) # = 2 gamma(r)

# ---------------------------------------------------------------------------
# Asymptotic variances (Theorems 3.1-3.3)
# ---------------------------------------------------------------------------

# eq. (9)/(25): AVAR of sqrt(n)(H_hat - H); requires H < 3/4.
avar_hurst <- function(H, R = 20000L) {
  a <- 2 * H
  s1 <- series_sum(function(r) .g2(r, a)^2, R)
  s2 <- series_sum(function(r) (.apow(r + 2, a) + .apow(r - 2, a) - 2 * .apow(r, a))^2, R)
  s3 <- series_sum(function(r) (.apow(r + 1, a) + .apow(r - 2, a) - .apow(r, a) - .apow(r - 1, a))^2, R)
  (4 + s1 + 2^(-4 * H) * s2 - 2^(1 - 2 * H) * s3) / (4 * log(2)^2)
}

# eq. (10): standard error of sigma2_hat.
# The CLT is stated for Delta = 1/n; the proof uses log(1/Delta), which is the
# correct factor when n and Delta are chosen separately (Tables 4-5 use log(1/Delta)).
se_sigma2 <- function(H, sigma2, n, delta, R = 20000L) {
  2 * sigma2 * log(1 / delta) * sqrt(avar_hurst(H, R) / n)
}

# eq. (14): AVAR of sqrt(n)(rho_hat - rho) for a time-reversible bivariate fBm.
avar_rho <- function(H1, H2, rho, R = 20000L) {
  ups1 <- function(a, b) 0.5 * series_sum(function(r) .g2(r, a + b)^2, R)
  ups2 <- 0.5 * series_sum(function(r) .g2(r, 2 * H1) * .g2(r, 2 * H2), R)
  ups3 <- function(a, b) series_sum(function(r) .g2(r, 2 * a) * .g2(r, a + b), R)
  (1 - rho^2)^2 +
    rho^2 * ((1 + rho^2) * ups1(H1, H2) + ups1(H1, H1) / 2 + ups1(H2, H2) / 2 -
               ups3(H1, H2) - ups3(H2, H1)) +
    ups2
}

# eq. (27): AVAR of sqrt(n) eta_hat under H0: eta = 0.
avar_eta <- function(H1, H2, rho, R = 20000L) {
  a <- H1 + H2
  g1 <- (2^(2 * H1) - 2) / 2
  g2 <- (2^(2 * H2) - 2) / 2
  g12 <- (2^a - 2) / 2
  P <- function(r, e) .apow(r, e)
  s1 <- series_sum(function(r) {
    (P(r + 1, a) + P(r - 1, a) - 2 * P(r, a))^2 -
      (P(r + 2, a) + P(r, a) - 2 * P(r + 1, a)) * (P(r, a) + P(r - 2, a) - 2 * P(r - 1, a))
  }, R)
  A <- function(r, h) P(r + 1, 2 * h) + P(r - 1, 2 * h) - 2 * P(r, 2 * h)
  D <- function(r, h) P(r + 2, 2 * h) + P(r, 2 * h) - 2 * P(r + 1, 2 * h)
  E <- function(r, h) P(r, 2 * h) + P(r - 2, 2 * h) - 2 * P(r - 1, 2 * h)
  s2 <- series_sum(function(r) {
    2 * A(r, H1) * A(r, H2) - D(r, H1) * E(r, H2) - D(r, H2) * E(r, H1)
  }, R)
  (2 * (1 - g1 * g2) + 2 * rho^2 * (g12^2 - 1) - rho^2 * s1 + 0.5 * s2) / (2^a - 2)^2
}

# Standard errors of all MM estimates (bivariate pairs for rho and eta).
mfbm_standard_errors <- function(est, R = 20000L) {
  d <- length(est$H)
  n <- est$n
  se_H <- vapply(est$H, function(h) sqrt(avar_hurst(h, R) / n), 0)
  se_s2 <- vapply(seq_len(d), function(j) se_sigma2(est$H[j], est$sigma2[j], n, est$delta, R), 0)
  se_rho <- se_eta <- matrix(NA_real_, d, d, dimnames = dimnames(est$rho))
  for (i in seq_len(d)) for (j in seq_len(d)) if (i != j) {
    se_rho[i, j] <- sqrt(avar_rho(est$H[i], est$H[j], est$rho[i, j], R) / n)
    se_eta[i, j] <- sqrt(avar_eta(est$H[i], est$H[j], est$rho[i, j], R) / n)
  }
  names(se_s2) <- names(se_H)
  list(H = se_H, sigma2 = se_s2, rho = se_rho, eta = se_eta)
}

# ---------------------------------------------------------------------------
# Generic Isserlis-based asymptotic covariances (cross-check of the closed
# forms and the covariance matrix needed by the optimal GMM, eqs. (42)-(45))
# ---------------------------------------------------------------------------

# kappa(j) = Cov(Delta_{k,l1} B^p, Delta_{k-j,l2} B^q) at Delta = 1, sigma = 1,
# for the time-reversible mfBm with correlation r and exponents Hp, Hq.
.lagged_inc_cov <- function(j, Hp, Hq, r, l1, l2) {
  f <- function(s, t) mfbm_cross_cov_unit(s, t, Hp, Hq, r, 0)
  f(j, 0) - f(j, -l2) - f(j - l1, 0) + f(j - l1, -l2)
}

# Two-sided sum over j in Z with tail corrections on both sides.
.sum_z <- function(v0, fpos, fneg, R) {
  v0 + series_sum(fpos, R) + series_sum(fneg, R)
}

# Asymptotic covariance matrix C (L x L) of sqrt(n)(rho_hat_1, ..., rho_hat_L),
# eq. (46), via Isserlis' theorem and the delta method at the true moments.
gmm_rho_cov <- function(H1, H2, rho, L, R = 5000L) {
  Hs <- c(H1, H2)
  rr <- function(p, q) if (p == q) 1 else rho
  # statistic index: 1 = QV1, 2 = QV2, 3 = QC, each a pair of components
  pairs <- list(c(1, 1), c(2, 2), c(1, 2))
  kap <- function(j, p, q, l1, l2) .lagged_inc_cov(j, Hs[p], Hs[q], rr(p, q), l1, l2)
  av_entry <- function(a, b, l1, l2) {
    p <- pairs[[a]][1]; q <- pairs[[a]][2]
    r <- pairs[[b]][1]; s <- pairs[[b]][2]
    term <- function(j) {
      kap(j, p, r, l1, l2) * kap(j, q, s, l1, l2) + kap(j, p, s, l1, l2) * kap(j, q, r, l1, l2)
    }
    .sum_z(term(0), function(j) term(j), function(j) term(-j), R)
  }
  psi <- function(l) {
    y1 <- l^(2 * H1); y2 <- l^(2 * H2); x <- rho * l^(H1 + H2)
    c(-x / (2 * y1^1.5 * sqrt(y2)), -x / (2 * y2^1.5 * sqrt(y1)), 1 / sqrt(y1 * y2))
  }
  C <- matrix(0, L, L)
  for (l1 in seq_len(L)) for (l2 in l1:L) {
    AV <- matrix(0, 3, 3)
    for (a in 1:3) for (b in 1:3) AV[a, b] <- av_entry(a, b, l1, l2)
    C[l1, l2] <- C[l2, l1] <- drop(t(psi(l1)) %*% AV %*% psi(l2))
  }
  C
}

# Asymptotic variance of sqrt(n) eta_hat under H0 by direct Isserlis summation
# (independent check of eq. (27)).
avar_eta_isserlis <- function(H1, H2, rho, R = 5000L) {
  par <- mfbm_params(c(H1, H2), rho, 0)
  # c_pq(l) = Cov(X^p_{t+l}, X^q_t); time-reversible => c_pq(l) = c_pq(-l)
  cf <- function(l, p, q) {
    r <- if (p == q) 1 else rho
    mfbm_increment_cov_unit(l, par$H[p], par$H[q], r, 0)
  }
  # T_k = X2_{k+1} X1_k - X1_{k+1} X2_k, Cov(T_0, T_j) by Isserlis
  cov4 <- function(a, i, b, jj, c, k, dd, l) {
    cf(i - k, a, c) * cf(jj - l, b, dd) + cf(i - l, a, dd) * cf(jj - k, b, c)
  }
  term <- function(j) {
    cov4(2, 1, 1, 0, 2, j + 1, 1, j) - cov4(2, 1, 1, 0, 1, j + 1, 2, j) -
      cov4(1, 1, 2, 0, 2, j + 1, 1, j) + cov4(1, 1, 2, 0, 1, j + 1, 2, j)
  }
  v <- .sum_z(term(0), function(j) term(j), function(j) term(-j), R)
  v / (2^(H1 + H2) - 2)^2
}

# ---------------------------------------------------------------------------
# Appendix G: lagged and second-order correlation estimators, optimal GMM
# ---------------------------------------------------------------------------

# eq. (36): rho_hat_lag with lag l (l = 1 gives eq. (8a)); x are increments.
rho_lag <- function(x1, x2, l) {
  if (l == 1) return(rho_mm(x1, x2))
  y1 <- stats::filter(x1, rep(1, l), sides = 1)[l:length(x1)]
  y2 <- stats::filter(x2, rep(1, l), sides = 1)[l:length(x2)]
  sum(y1 * y2) / sqrt(sum(y1^2) * sum(y2^2))
}

rho_lags <- function(x1, x2, L) {
  c1 <- c(0, cumsum(x1)); c2 <- c(0, cumsum(x2))
  n <- length(x1)
  vapply(seq_len(L), function(l) {
    y1 <- c1[(l + 1):(n + 1)] - c1[1:(n + 1 - l)]
    y2 <- c2[(l + 1):(n + 1)] - c2[1:(n + 1 - l)]
    sum(y1 * y2) / sqrt(sum(y1^2) * sum(y2^2))
  }, 0)
}

# eq. (41): correlation of second-order increments.
rho_second_order <- function(x1, x2) {
  rho_mm(diff(x1), diff(x2))
}

# eq. (38): optimal linear combination given a covariance matrix C.
gmm_weights <- function(C) {
  w <- solve(C, rep(1, nrow(C)))
  w / sum(w)
}

# Pilot Hurst estimates for plug-in steps, clipped to [0.001, 0.999] as in the
# authors' estMFBM (with H1 = 0.1 the MM estimate is occasionally negative).
.pilot_hurst <- function(x) min(max(hurst_mm(x), 0.001), 0.999)

# Feasible two-step rho_opt: pilot H1, H2, rho_1 plugged into C.
rho_opt <- function(x1, x2, L, R = 5000L) {
  rl <- rho_lags(x1, x2, L)
  H1 <- .pilot_hurst(x1); H2 <- .pilot_hurst(x2)
  C <- gmm_rho_cov(H1, H2, rl[1], L, R)
  w <- gmm_weights(C)
  list(estimate = sum(w * rl), weights = w, rho_lags = rl, avar = 1 / sum(solve(C, rep(1, L))))
}

# Bootstrap version rho_opt_boot (Appendix G.3): C estimated from B parametric
# bootstrap samples drawn from the fitted model.
rho_opt_boot <- function(x1, x2, L, B, delta) {
  n <- length(x1)
  rl <- rho_lags(x1, x2, L)
  H <- c(.pilot_hurst(x1), .pilot_hurst(x2))
  s2 <- c(sigma2_mm(x1, delta, H[1]), sigma2_mm(x2, delta, H[2]))
  r0 <- max(min(rl[1], 0.99 * rho_max(H[1], H[2])), -0.99 * rho_max(H[1], H[2]))
  emb <- mfbm_embedding(mfbm_params(H, r0, 0, sqrt(s2)), n, delta)
  Xb <- mfbm_simulate_increments(emb, B)
  R <- t(vapply(seq_len(B), function(b) rho_lags(Xb[, 1, b], Xb[, 2, b], L), numeric(L)))
  C <- n * stats::cov(R)
  w <- gmm_weights(C)
  list(estimate = sum(w * rl), weights = w)
}

# ---------------------------------------------------------------------------
# Amblard and Coeurjolly (2011) estimator (Appendix A)
# ---------------------------------------------------------------------------
# Implementation choices (the authors used SimEstFBM.R with default settings,
# which is not distributed with the paper):
#   * filter a = (1, -2, 1) (second-order difference), dilations m in M = 1..5;
#   * weights w_c = w_dd = 0, so (H_i, sigma_i) solve the log-variance
#     regression over m, and (rho, eta) are then fitted to the log-absolute
#     moment equations eps_c and eps_d given (H, sigma);
#   * the log-absolute criterion identifies |eta| only, so eta_hat >= 0
#     (this is the source of the large positive bias at eta = 0 in Table 6);
#     the sign of rho is taken from the empirical lag-0 cross-covariance
#     (estMFBM returns |rho|; identical for the positive rho of Table 6).
# These choices were later confirmed against the authors' BYZestimator.R, which
# calls estMFBM(nma = "i2", M1 = 1, M2 = 5, w = c(1, 0, 0)); see
# scripts/00_check_author_code.R.

.ac_filter <- c(1, -2, 1)

.ac_filtered <- function(B, m, a = .ac_filter, author_quirk = TRUE) {
  # B: levels (n + 1 values incl. B_0); B^m_k = sum_t a_t B_{k - m t}
  # author_quirk: Coeurjolly's dilatation() returns a filter of length m * (l + 1) - 1
  # with trailing zeros, so stats::filter() drops m * (l + 1) - 2 instead of m * l
  # leading values; reproduced by default to match the authors' estMFBM exactly.
  l <- length(a) - 1
  N <- length(B)
  skip <- if (author_quirk && m > 1) m * (l + 1) - 2 else m * l
  idx <- (skip + 1):N
  out <- 0
  for (t in 0:l) out <- out + a[t + 1] * B[idx - m * t]
  out
}

# -1/2 sum_{t,l} a_t a_l w(h + m(t - l)) |h + m(t - l)|^H  (unit Delta)
.ac_kernel <- function(h, m, H, weight = function(u) 1, a = .ac_filter) {
  l <- length(a) - 1
  s <- 0
  for (t in 0:l) for (u in 0:l) {
    z <- h + m * (t - u)
    s <- s + a[t + 1] * a[u + 1] * weight(z) * .apow(z, H)
  }
  -0.5 * s
}

estimate_ac <- function(X, delta, dilations = 1:5, a = .ac_filter) {
  X <- as.matrix(X)
  d <- ncol(X)
  B <- rbind(0, apply(X, 2, cumsum))
  l <- length(a) - 1
  Fm <- lapply(dilations, function(m) {
    apply(B, 2, .ac_filtered, m = m, a = a)
  })
  emp_cov <- function(Y, i, j, h) {
    nn <- nrow(Y)
    if (h >= nn) return(NA_real_)
    sum(Y[seq_len(nn - h), i] * Y[seq_len(nn - h) + h, j]) / (nn - h)
  }
  H <- s2 <- numeric(d)
  lm_ <- log(dilations)
  for (i in seq_len(d)) {
    lc <- vapply(Fm, function(Y) log(emp_cov(Y, i, i, 0)), 0)
    fit <- stats::lm.fit(cbind(1, lm_), lc)$coefficients
    H[i] <- min(max(fit[2] / 2, 0.001), 0.999) # clipped as in estMFBM
    K <- vapply(dilations, function(m) .ac_kernel(0, m, 2 * H[i], a = a), 0)
    s2[i] <- exp(mean(lc - 2 * H[i] * log(dilations * delta) - log(K / dilations^(2 * H[i]))))
  }
  rho <- diag(d)
  eta <- matrix(0, d, d)
  for (i in seq_len(d)) for (j in seq_len(d)) if (i < j) {
    Hs <- H[i] + H[j]
    ss <- sqrt(s2[i] * s2[j])
    c0 <- vapply(Fm, function(Y) emp_cov(Y, i, j, 0), 0)
    k0 <- ss * delta^Hs * vapply(dilations, function(m) .ac_kernel(0, m, Hs, a = a), 0)
    ratio <- c0 / k0
    r_abs <- exp(mean(log(abs(ratio))))
    rho[i, j] <- rho[j, i] <- sign(sum(ratio)) * min(r_abs, 1)
    # antisymmetric part at lag h = m * l
    dd <- vapply(seq_along(dilations), function(k) {
      Y <- Fm[[k]]; h <- dilations[k] * l
      0.5 * abs(emp_cov(Y, i, j, h) - emp_cov(Y, j, i, h))
    }, 0)
    kd <- vapply(dilations, function(m) {
      # COV_ij(h) - COV_ji(h) = eta * sigma_i sigma_j * 2 * (-1/2) sum a a (-sign) |.|^H
      0.5 * abs(2 * .ac_kernel(m * l, m, Hs, weight = function(u) -sign(u), a = a))
    }, 0) * ss * delta^Hs
    eta[i, j] <- exp(mean(log(dd / kd)))
    eta[j, i] <- -eta[i, j]
  }
  list(H = H, sigma2 = s2, rho = rho, eta = eta, n = nrow(X), delta = delta)
}
