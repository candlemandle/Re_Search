# Test of time-reversibility H0: eta_{1,2} = 0 (Corollary 3.1, eq. (11)).

# x1, x2: increments of the two components.
# Returns the statistic sqrt(n)|eta_hat| / sqrt(AVAR_eta(H1_hat, H2_hat, rho_hat)),
# its two-sided p-value and rejection decisions at the requested levels.
time_reversibility_test <- function(x1, x2, alpha = c(0.01, 0.05), R = 20000L) {
  n <- length(x1)
  H1 <- hurst_mm(x1)
  H2 <- hurst_mm(x2)
  rho <- rho_mm(x1, x2)
  eta <- eta_mm(x1, x2)
  # The CLT needs 0 < H < 3/4 and H1 + H2 != 1. In small samples H_hat can be
  # negative (e.g. H1 = 0.1, n = 500); then AVAR is undefined and the test is
  # reported as invalid (statistic NA) instead of silently returning NaN.
  valid <- min(H1, H2) > 0 && max(H1, H2) < 0.75 && abs(H1 + H2 - 1) > 1e-8
  av <- if (valid) avar_eta(H1, H2, rho, R) else NA_real_
  stat <- sqrt(n) * abs(eta) / sqrt(av)
  out <- list(
    n = n, H1 = H1, H2 = H2, rho = rho, eta = eta,
    se_eta = sqrt(av / n), statistic = stat,
    p_value = 2 * stats::pnorm(-stat),
    valid = valid
  )
  for (a in alpha) out[[sprintf("reject_%g", a)]] <- stat > stats::qnorm(1 - a / 2)
  out
}

# Pairwise tests for all columns of a d-variate increment matrix.
time_reversibility_table <- function(X, alpha = c(0.01, 0.05), R = 20000L) {
  X <- as.matrix(X)
  nm <- colnames(X)
  if (is.null(nm)) nm <- paste0("B", seq_len(ncol(X)))
  rows <- list()
  for (i in seq_len(ncol(X))) for (j in seq_len(ncol(X))) if (i < j) {
    ok <- stats::complete.cases(X[, i], X[, j])
    t <- time_reversibility_test(X[ok, i], X[ok, j], alpha, R)
    rows[[length(rows) + 1]] <- data.frame(
      asset1 = nm[i], asset2 = nm[j], n = t$n, eta = t$eta, se_eta = t$se_eta,
      statistic = t$statistic, p_value = t$p_value,
      t[grep("^reject_", names(t))], valid = t$valid, check.names = FALSE
    )
  }
  do.call(rbind, rows)
}
