# Covariance structure of the multivariate fractional Brownian motion (mfBm)
# and exact simulation of its increments.
#
# Paper references: eq. (1) w(t, h, H); eq. (4) general bivariate cross-covariance
# (Amblard et al. 2013, Prop. 3); eq. (5) time-reversible case; eq. (18)/(26)
# covariances of increments; Section 2 rho_max(H1, H2).

# eq. (1): Cov(B^H_{t+h}, B^H_t) = w(t, h, H)
w_fbm <- function(t, h, H) {
  0.5 * (abs(t + h)^(2 * H) + abs(t)^(2 * H) - abs(h)^(2 * H))
}

# Absolute power with the convention 0^a = 0 for a > 0.
.apow <- function(x, a) abs(x)^a

# eq. (4) and the H1 + H2 = 1 case: Cov(B^(1)_s, B^(2)_t) / (sigma1 sigma2).
# Components of the same process use rho = 1, eta = 0, H1 = H2.
mfbm_cross_cov_unit <- function(s, t, H1, H2, rho, eta = 0) {
  Hs <- H1 + H2
  if (abs(Hs - 1) > 1e-12) {
    0.5 * ((rho + eta * sign(s)) * .apow(s, Hs) +
             (rho - eta * sign(t)) * .apow(t, Hs) -
             (rho - eta * sign(t - s)) * .apow(t - s, Hs))
  } else {
    xlogx <- function(x) ifelse(x == 0, 0, x * log(abs(x)))
    0.5 * (rho * (abs(s) + abs(t) - abs(s - t)) +
             eta * (xlogx(t) - xlogx(s) - xlogx(t - s)))
  }
}

# Cross-covariance of unit-step increments:
# Cov(Delta_{j+k} B^(p), Delta_j B^(q)) / (sigma_p sigma_q) at Delta = 1.
# For p = q (rho = 1, eta = 0) this is the fGn autocovariance gamma_p(k), eq. (19).
mfbm_increment_cov_unit <- function(k, Hp, Hq, rho, eta = 0) {
  f <- function(s, t) mfbm_cross_cov_unit(s, t, Hp, Hq, rho, eta)
  f(k, 0) - f(k, -1) - f(k - 1, 0) + f(k - 1, -1)
}

# eq. (19): gamma_{p,q}(l) for the time-reversible case (exponent a = Hp + Hq).
gamma_pq <- function(l, a) {
  0.5 * (.apow(l + 1, a) + .apow(l - 1, a) - 2 * .apow(l, a))
}

# Validate and complete a parameter set for a d-dimensional mfBm.
#   H      : length-d Hurst exponents
#   rho    : d x d correlation matrix (or scalar for d = 2)
#   eta    : d x d antisymmetric asymmetry matrix (or scalar eta_{1,2} for d = 2)
#   sigma  : length-d scale parameters (standard deviations of B_1)
mfbm_params <- function(H, rho, eta = 0, sigma = rep(1, length(H))) {
  d <- length(H)
  if (any(H <= 0 | H >= 1)) stop("Hurst exponents must lie in (0, 1)")
  if (length(sigma) == 1) sigma <- rep(sigma, d)
  if (!is.matrix(rho)) {
    if (d != 2) stop("scalar rho only allowed for d = 2")
    rho <- matrix(c(1, rho, rho, 1), 2, 2)
  }
  if (!is.matrix(eta)) {
    if (d == 2) {
      eta <- matrix(c(0, -eta, eta, 0), 2, 2) # eta[1, 2] = eta_{1,2}
    } else if (all(eta == 0)) {
      eta <- matrix(0, d, d)
    } else {
      stop("scalar eta only allowed for d = 2")
    }
  }
  if (max(abs(rho - t(rho))) > 1e-12) stop("rho must be symmetric")
  if (max(abs(eta + t(eta))) > 1e-12) stop("eta must be antisymmetric")
  if (any(abs(diag(rho) - 1) > 1e-12)) stop("diag(rho) must be 1")
  if (any(abs(rho) > 1 + 1e-12)) stop("|rho| must be <= 1")
  list(H = H, rho = rho, eta = eta, sigma = sigma, d = d)
}

# Autocovariance sequence of the d-variate increment process on grid Delta:
# G[p, q, k + 1] = Cov(Delta_{j+k} B^(p), Delta_j B^(q)), k = 0..max_lag.
mfbm_increment_acf <- function(par, max_lag, delta = 1) {
  d <- par$d
  k <- 0:max_lag
  G <- array(0, c(d, d, max_lag + 1))
  for (p in seq_len(d)) for (q in seq_len(d)) {
    if (p == q) {
      g <- mfbm_increment_cov_unit(k, par$H[p], par$H[p], 1, 0)
    } else {
      g <- mfbm_increment_cov_unit(k, par$H[p], par$H[q], par$rho[p, q], par$eta[p, q])
    }
    G[p, q, ] <- g * par$sigma[p] * par$sigma[q] * delta^(par$H[p] + par$H[q])
  }
  G
}

# Section 2: maximal admissible |rho| of a bivariate time-reversible fBm.
rho_max <- function(H1, H2) {
  H <- (H1 + H2) / 2
  sqrt(sin(pi * H1) * sin(pi * H2) * gamma(2 * H1 + 1) * gamma(2 * H2 + 1)) /
    (sin(pi * H) * gamma(2 * H + 1))
}

# Full covariance matrix of the stacked increments (X_1, ..., X_n), X_k in R^d,
# ordered as (X_1^(1), ..., X_1^(d), X_2^(1), ...). Used for Cholesky checks.
mfbm_increment_cov_matrix <- function(par, n, delta = 1) {
  d <- par$d
  G <- mfbm_increment_acf(par, n - 1, delta)
  S <- matrix(0, n * d, n * d)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    k <- i - j
    blk <- if (k >= 0) G[, , k + 1] else t(G[, , -k + 1])
    S[(i - 1) * d + seq_len(d), (j - 1) * d + seq_len(d)] <- blk
  }
  S
}

# Covariance matrix of levels (B_{t_1}, ..., B_{t_m}) of a time-reversible mfBm,
# stacked as in Section 4 (all components at t_1, then at t_2, ...). eq. (5).
mfbm_level_cov_matrix <- function(par, times) {
  d <- par$d
  m <- length(times)
  S <- matrix(0, m * d, m * d)
  st <- expand.grid(s = times, t = times)
  for (p in seq_len(d)) for (q in seq_len(d)) {
    r <- if (p == q) 1 else par$rho[p, q]
    e <- if (p == q) 0 else par$eta[p, q]
    blk <- par$sigma[p] * par$sigma[q] *
      mfbm_cross_cov_unit(st$s, st$t, par$H[p], par$H[q], r, e)
    S[(seq_len(m) - 1) * d + p, (seq_len(m) - 1) * d + q] <- matrix(blk, m, m)
  }
  S
}

# ---------------------------------------------------------------------------
# Exact simulation by multivariate circulant embedding (Wood and Chan 1994;
# Helgason, Pipiras and Abry 2011 for mfGn). The factorisation depends only on
# (par, n, delta), so it is computed once and reused across replications.
# ---------------------------------------------------------------------------

# The embedding length 2m (m >= n) is doubled until the circulant is PSD
# (at most `max_doublings` times); near the admissibility boundary a longer
# embedding is often needed.
# If no length gives an exactly PSD circulant, the approximate Wood-Chan method
# is used: eigenvalues with |lambda| < approx_tol * max(lambda) are set to zero
# (flagged by emb$approximate = TRUE). This happens only very close to the
# admissibility boundary, e.g. eta = 0.65 in Table 7 (boundary 0.693).
mfbm_embedding <- function(par, n, delta = 1, tol = 1e-10, max_doublings = 2L,
                           approx_tol = 1e-4) {
  m <- n
  for (k in 0:max_doublings) {
    emb <- try(.mfbm_embedding_m(par, n, m, delta, tol), silent = TRUE)
    if (!inherits(emb, "try-error")) return(emb)
    m <- 2L * m
  }
  emb <- .mfbm_embedding_m(par, n, n, delta, approx_tol)
  emb$approximate <- TRUE
  emb
}

# Admissibility of a bivariate mfBm (Amblard et al. 2013, Prop. 9), H1 + H2 != 1:
# rho^2 sin^2(pi/2 (H1+H2)) + eta^2 cos^2(pi/2 (H1+H2))
#   <= Gamma(2H1+1) Gamma(2H2+1) sin(pi H1) sin(pi H2) / Gamma(H1+H2+1)^2.
mfbm_admissible <- function(H1, H2, rho, eta = 0) {
  Hs <- H1 + H2
  lhs <- rho^2 * sin(pi * Hs / 2)^2 + eta^2 * cos(pi * Hs / 2)^2
  rhs <- gamma(2 * H1 + 1) * gamma(2 * H2 + 1) * sin(pi * H1) * sin(pi * H2) / gamma(Hs + 1)^2
  lhs <= rhs
}

.mfbm_embedding_m <- function(par, n, m, delta, tol) {
  d <- par$d
  M <- 2L * m
  G <- mfbm_increment_acf(par, m, delta)
  # circulant first row c(j), j = 0..M-1 (d x d blocks)
  C <- array(0, c(d, d, M))
  for (j in 0:(m - 1)) C[, , j + 1] <- G[, , j + 1]
  C[, , m + 1] <- 0.5 * (G[, , m + 1] + t(G[, , m + 1]))
  for (j in (m + 1):(M - 1)) C[, , j + 1] <- t(G[, , M - j + 1])
  Fq <- array(0 + 0i, c(d, d, M))
  for (p in seq_len(d)) for (q in seq_len(d)) Fq[p, q, ] <- fft(C[p, q, ])
  A <- array(0 + 0i, c(d, d, M))
  min_eig <- Inf
  max_eig <- 0
  for (l in seq_len(M)) {
    Fl <- Fq[, , l, drop = TRUE]
    if (d == 1) Fl <- matrix(Fl, 1, 1)
    Fl <- 0.5 * (Fl + Conj(t(Fl)))
    if (max(abs(Im(Fl))) < 1e-14 * max(1, max(abs(Fl)))) Fl <- Re(Fl)
    e <- eigen(Fl, symmetric = TRUE)
    lam <- Re(e$values)
    min_eig <- min(min_eig, lam)
    max_eig <- max(max_eig, lam)
    A[, , l] <- e$vectors %*% diag(sqrt(pmax(lam, 0)), d, d)
  }
  if (min_eig < -tol * max_eig) {
    stop(sprintf("circulant embedding not PSD (min eigenvalue %.3g): parameters may be inadmissible",
                 min_eig))
  }
  list(A = A, n = n, M = M, d = d, delta = delta, par = par, min_eig = min_eig,
       rel_min_eig = min_eig / max_eig, approximate = FALSE)
}

# Simulate nsim independent increment paths: array n x d x nsim.
# Each complex draw yields two independent real paths (real and imaginary parts).
mfbm_simulate_increments <- function(emb, nsim = 1L) {
  d <- emb$d
  M <- emb$M
  n <- emb$n
  ncplx <- ceiling(nsim / 2)
  out <- array(0, c(n, d, 2 * ncplx))
  # Z[q] is an M x ncplx complex matrix with independent N(0,1) real/imag parts
  Z <- lapply(seq_len(d), function(q) {
    matrix(complex(real = rnorm(M * ncplx), imaginary = rnorm(M * ncplx)), M, ncplx)
  })
  for (p in seq_len(d)) {
    V <- matrix(0 + 0i, M, ncplx)
    for (q in seq_len(d)) V <- V + emb$A[p, q, ] * Z[[q]]
    Y <- mvfft(V, inverse = TRUE)[seq_len(n), , drop = FALSE] / sqrt(M)
    out[, p, seq_len(ncplx)] <- Re(Y)
    out[, p, ncplx + seq_len(ncplx)] <- Im(Y)
  }
  out[, , seq_len(nsim), drop = FALSE]
}

# Convenience wrapper: n increments of a d-variate mfBm on grid Delta.
# Returns a list of n x d increment matrices (length nsim).
simulate_mfbm <- function(n, H, rho, eta = 0, sigma = rep(1, length(H)),
                          delta = 1 / n, nsim = 1L, emb = NULL) {
  if (is.null(emb)) emb <- mfbm_embedding(mfbm_params(H, rho, eta, sigma), n, delta)
  X <- mfbm_simulate_increments(emb, nsim)
  lapply(seq_len(nsim), function(s) X[, , s, drop = TRUE])
}

# Levels B_{0}, B_{Delta}, ..., B_{n Delta} from increments (B_0 = 0).
increments_to_levels <- function(X) {
  X <- as.matrix(X)
  rbind(0, apply(X, 2, cumsum))
}

# Reference simulator by Cholesky factorisation (small n; used in tests).
simulate_mfbm_cholesky <- function(n, par, delta = 1, nsim = 1L) {
  S <- mfbm_increment_cov_matrix(par, n, delta)
  R <- chol(S)
  lapply(seq_len(nsim), function(s) {
    x <- drop(crossprod(R, rnorm(n * par$d)))
    matrix(x, n, par$d, byrow = TRUE)
  })
}
