# Tests for Code A: covariance structure, simulation, estimators, asymptotic
# variances, time-reversibility test, data processing and reproducibility.

# --- covariance structure ------------------------------------------------------

test_that("fBm covariance (1) and fGn autocovariance (19) are correct", {
  expect_equal(w_fbm(1, 0, 0.3), 1)
  expect_equal(w_fbm(2, 1, 0.5), 2) # Brownian motion: Cov(B_3, B_2) = 2
  expect_equal(gamma_pq(0, 0.6), 1)
  expect_equal(mfbm_increment_cov_unit(0:5, 0.3, 0.3, 1, 0), gamma_pq(0:5, 0.6))
  # Brownian increments are uncorrelated
  expect_equal(mfbm_increment_cov_unit(1:5, 0.5, 0.5, 1, 0), rep(0, 5))
})

test_that("increment covariance matrices are symmetric and positive definite", {
  for (eta in c(0, 0.3)) {
    S <- mfbm_increment_cov_matrix(mfbm_params(c(0.1, 0.4), 0.4, eta), 60, 1 / 250)
    expect_equal(S, t(S), tolerance = 1e-14)
    expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  S3 <- mfbm_increment_cov_matrix(
    mfbm_params(c(0.1, 0.25, 0.4), matrix(c(1, .3, .2, .3, 1, .4, .2, .4, 1), 3)), 30)
  expect_equal(S3, t(S3), tolerance = 1e-14)
  expect_gt(min(eigen(S3, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("level covariance (5) is consistent with the increment covariance", {
  par <- mfbm_params(c(0.1, 0.4), 0.4, 0, c(1, 2))
  n <- 40
  S <- mfbm_level_cov_matrix(par, (1:n) / 250)
  expect_equal(S, t(S), tolerance = 1e-14)
  expect_equal(diag(S)[c(1, 2)], c(1, 4) * (1 / 250)^(2 * c(0.1, 0.4)))
  # increments D = A B with B_0 = 0 must reproduce the increment covariance matrix
  A <- diag(n * 2)
  for (k in 2:n) A[(k - 1) * 2 + 1:2, (k - 2) * 2 + 1:2] <- -diag(2)
  expect_equal(A %*% S %*% t(A), mfbm_increment_cov_matrix(par, n, 1 / 250), tolerance = 1e-12)
})

test_that("rho_max reproduces the values quoted in Section 2", {
  expect_equal(round(rho_max(0.2, 0.8), 3), 0.662)
  expect_equal(round(rho_max(0.1, 0.9), 3), 0.383)
  expect_equal(rho_max(0.3, 0.3), 1)
  expect_equal(rho_max(0.1, 0.4), rho_max(0.4, 0.1))
})

test_that("admissibility condition agrees with rho_max and with the exact covariance", {
  expect_true(mfbm_admissible(0.1, 0.4, 0.99 * rho_max(0.1, 0.4)))
  expect_false(mfbm_admissible(0.1, 0.4, 1.01 * rho_max(0.1, 0.4)))
  # boundary for eta at rho = 0.4 is 0.693 (Table 7 uses eta up to 0.65)
  expect_true(mfbm_admissible(0.1, 0.4, 0.4, 0.69))
  expect_false(mfbm_admissible(0.1, 0.4, 0.4, 0.71))
  bad <- mfbm_increment_cov_matrix(mfbm_params(c(0.1, 0.4), 0.4, 0.75), 300, 1 / 250)
  expect_lt(min(eigen(bad, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("invalid parameters are rejected", {
  expect_error(mfbm_params(c(0.1, 1.2), 0.3), "Hurst")
  expect_error(mfbm_params(c(0.1, 0.4, 0.3), matrix(c(1, .2, .1, .3, 1, .2, .1, .2, 1), 3)), "symmetric")
  expect_error(mfbm_params(c(0.1, 0.4), 1.5), "rho")
})

# --- simulation ------------------------------------------------------------------

test_that("circulant embedding reproduces the target increment covariances", {
  set.seed(11)
  par <- mfbm_params(c(0.1, 0.4), 0.4, 0.5)
  emb <- mfbm_embedding(par, 400, 1 / 250)
  expect_false(emb$approximate)
  X <- mfbm_simulate_increments(emb, 4000)
  G <- mfbm_increment_acf(par, 2, 1 / 250)
  for (k in 0:2) for (p in 1:2) for (q in 1:2) {
    prod <- X[(1 + k):400, p, ] * X[1:(400 - k), q, ]
    m <- mean(prod)
    se <- stats::sd(colMeans(prod)) / sqrt(dim(X)[3])
    expect_lt(abs(m - G[p, q, k + 1]), 5 * se + 1e-6)
  }
})

test_that("embedding and Cholesky simulators agree in distribution", {
  par <- mfbm_params(c(0.2, 0.35), -0.5, 0)
  set.seed(3)
  a <- simulate_mfbm(50, par$H, par$rho, 0, delta = 1 / 50, nsim = 3000)
  set.seed(3)
  b <- simulate_mfbm_cholesky(50, par, 1 / 50, nsim = 3000)
  ra <- vapply(a, function(x) rho_mm(x[, 1], x[, 2]), 0)
  rb <- vapply(b, function(x) rho_mm(x[, 1], x[, 2]), 0)
  expect_lt(abs(mean(ra) - mean(rb)), 4 * sqrt(stats::var(ra) / 3000 + stats::var(rb) / 3000))
})

test_that("near-boundary parameters fall back to the flagged approximate embedding", {
  emb <- mfbm_embedding(mfbm_params(c(0.1, 0.4), 0.4, 0.65), 1000, 1 / 250)
  expect_true(emb$approximate)
  expect_gt(emb$rel_min_eig, -1e-4)
})

test_that("same seed gives identical paths; levels start at zero", {
  set.seed(5); a <- simulate_mfbm(100, c(0.1, 0.4), 0.3)[[1]]
  set.seed(5); b <- simulate_mfbm(100, c(0.1, 0.4), 0.3)[[1]]
  expect_identical(a, b)
  lv <- increments_to_levels(a)
  expect_equal(lv[1, ], c(0, 0))
  expect_equal(nrow(lv), 101)
})

# --- point estimators ------------------------------------------------------------

test_that("estimators recover parameters and lie in admissible ranges", {
  set.seed(21)
  emb <- mfbm_embedding(mfbm_params(c(0.1, 0.4), 0.4, 0.3, c(1, 2)), 5000, 1 / 250)
  est <- lapply(seq_len(40), function(s) estimate_mfbm(mfbm_simulate_increments(emb, 1)[, , 1], 1 / 250))
  H <- rowMeans(sapply(est, `[[`, "H"))
  expect_equal(H, c(0.1, 0.4), tolerance = 0.01)
  expect_equal(mean(sapply(est, function(e) e$rho[1, 2])), 0.4, tolerance = 0.01)
  expect_equal(mean(sapply(est, function(e) e$eta[1, 2])), 0.3, tolerance = 0.03)
  expect_equal(mean(sapply(est, function(e) e$sigma2[2])), 4, tolerance = 0.3)
  e <- est[[1]]
  expect_true(all(e$H > 0 & e$H < 1))
  expect_true(all(abs(e$rho) <= 1))
  expect_equal(e$rho, t(e$rho))
  expect_equal(e$eta, -t(e$eta))
})

test_that("eta convention: (8b) as printed has the opposite sign of Appendix C.2", {
  set.seed(8)
  x <- simulate_mfbm(20000, c(0.1, 0.4), 0.4, 0.5, delta = 1 / 250)[[1]]
  expect_gt(eta_mm(x[, 1], x[, 2]), 0.4)
  expect_equal(eta_mm(x[, 1], x[, 2], "printed"), -eta_mm(x[, 1], x[, 2]))
  expect_equal(eta_mm(x[, 2], x[, 1]), -eta_mm(x[, 1], x[, 2]))
})

test_that("H, rho, eta are scale invariant; sigma2 scales quadratically", {
  set.seed(2)
  x <- simulate_mfbm(500, c(0.15, 0.3), 0.5, 0.1)[[1]]
  expect_equal(hurst_mm(3 * x[, 1]), hurst_mm(x[, 1]))
  expect_equal(rho_mm(2 * x[, 1], 5 * x[, 2]), rho_mm(x[, 1], x[, 2]))
  expect_equal(eta_mm(2 * x[, 1], 5 * x[, 2]), eta_mm(x[, 1], x[, 2]))
  expect_equal(sigma2_mm(3 * x[, 1], 0.01), 9 * sigma2_mm(x[, 1], 0.01))
})

test_that("Brownian motion gives H close to 1/2", {
  set.seed(9)
  expect_equal(hurst_mm(rnorm(2e5)), 0.5, tolerance = 0.01)
})

# --- asymptotic variances -----------------------------------------------------------

test_that("asymptotic standard errors match Tables 4-5 (parentheses)", {
  expect_equal(round(sqrt(avar_hurst(0.1) / 500), 4), 0.0431)
  expect_equal(round(sqrt(avar_hurst(0.4) / 1000), 4), 0.0248)
  expect_equal(round(se_sigma2(0.1, 1, 500, 1 / 52), 4), 0.3404)
  # the paper prints this quantity as 0.3876 (Table 4) and 0.3877 (Table 5): 4th digit is rounding
  expect_equal(se_sigma2(0.4, 1, 1000, 1 / 250), 0.2741, tolerance = 2e-4 / 0.2741)
  expect_equal(round(sqrt(avar_rho(0.1, 0.4, 0) / 500), 4), 0.0472)
  expect_equal(round(sqrt(avar_rho(0.1, 0.4, 0.4) / 1000), 4), 0.0279)
  expect_equal(round(sqrt(avar_eta(0.1, 0.4, 0) / 500), 4), 0.1137)
  expect_equal(round(sqrt(avar_eta(0.1, 0.4, 0.4) / 1000), 4), 0.0733)
})

test_that("closed forms agree with direct Isserlis computations", {
  for (p in list(c(0.1, 0.4, 0.4), c(0.25, 0.3, -0.6), c(0.6, 0.2, 0.2))) {
    expect_equal(avar_eta(p[1], p[2], p[3]), avar_eta_isserlis(p[1], p[2], p[3]), tolerance = 1e-6)
    expect_equal(gmm_rho_cov(p[1], p[2], p[3], 1)[1, 1], avar_rho(p[1], p[2], p[3]), tolerance = 1e-6)
  }
  # Remark 3.4: Brownian case reduces to (1 - rho^2)^2
  expect_equal(avar_rho(0.5, 0.5, 0.3), (1 - 0.09)^2)
  # tail correction makes the truncation point (almost) irrelevant; convergence
  # slows down as H -> 3/4 (terms ~ r^(4H - 4))
  expect_equal(avar_hurst(0.4, 2000), avar_hurst(0.4, 20000), tolerance = 1e-6)
  expect_equal(avar_hurst(0.7, 2000), avar_hurst(0.7, 20000), tolerance = 3e-3)
})

test_that("optimal GMM weights sum to one and do not increase the variance", {
  C <- gmm_rho_cov(0.1, 0.4, 0.4, 5)
  expect_equal(C, t(C))
  w <- gmm_weights(C)
  expect_equal(sum(w), 1)
  expect_lt(1 / sum(solve(C, rep(1, 5))), C[1, 1])
  expect_true(all(diff(diag(C)) > 0)) # Proposition G.1: AVAR increasing in the lag
})

# --- time-reversibility test ------------------------------------------------------

test_that("test statistic equals sqrt(n)|eta| / sqrt(AVAR_eta)", {
  set.seed(4)
  x <- simulate_mfbm(800, c(0.1, 0.4), 0.4, 0, delta = 1 / 250)[[1]]
  t <- time_reversibility_test(x[, 1], x[, 2])
  manual <- sqrt(800) * abs(eta_mm(x[, 1], x[, 2])) /
    sqrt(avar_eta(hurst_mm(x[, 1]), hurst_mm(x[, 2]), rho_mm(x[, 1], x[, 2])))
  expect_equal(t$statistic, manual)
  expect_equal(t$p_value, 2 * pnorm(-manual))
  expect_true(t$valid)
})

test_that("the test is flagged invalid (not NaN) when H_hat is outside (0, 3/4)", {
  set.seed(1)
  x1 <- c(1, -1)[rep(1:2, 200)] + rnorm(400, sd = 0.01) # strongly anti-persistent: H_hat < 0
  x2 <- rnorm(400)
  expect_lt(hurst_mm(x1), 0)
  t <- time_reversibility_test(x1, x2)
  expect_false(t$valid)
  expect_true(is.na(t$statistic))
  expect_true(is.na(t$reject_0.05))
})

test_that("the test rejects a strongly asymmetric process", {
  set.seed(6)
  x <- simulate_mfbm(2000, c(0.1, 0.4), 0.4, 0.6, delta = 1 / 250)[[1]]
  expect_lt(time_reversibility_test(x[, 1], x[, 2])$p_value, 0.01)
})

# --- data processing ---------------------------------------------------------------

test_that("Risk Lab files are parsed and cleaned correctly", {
  f <- testthat::test_path("fixtures", "risklab_sample.txt")
  df <- read_risklab(f, "AAPL")
  expect_equal(nrow(df), 10)
  expect_equal(df$date[1], as.Date("1996-01-02"))
  expect_equal(df$qmle_trade[1], 0.170735)
  cl <- clean_risklab(df, "qmle_trade", as.Date("1996-01-03"), as.Date("1996-01-10"))
  expect_true(all(cl$date >= as.Date("1996-01-03") & cl$date <= as.Date("1996-01-10")))
  expect_equal(cl$logvol, log(cl$vol))
  df$qmle_trade[3] <- 0 # invalid values are dropped
  expect_false(df$date[3] %in% clean_risklab(df, "qmle_trade", as.Date("1990-01-01"), Sys.Date())$date)
})

test_that("panel increments use only jointly observed dates", {
  p <- data.frame(date = as.Date("2020-01-01") + 0:5,
                  A = c(0, 1, NA, 3, 4, 5), B = c(0, 2, 4, 6, NA, 10))
  inc <- panel_increments(p, c("A", "B"))
  expect_equal(unname(inc$X[, "A"]), c(1, 2, 2))
  expect_equal(unname(inc$X[, "B"]), c(2, 4, 4))
  expect_false(anyNA(inc$X))
})

# --- configuration and reproducibility --------------------------------------------------

test_that("smoke and full configurations differ only in sizes", {
  s <- code_a_config("smoke"); f <- code_a_config("full")
  expect_identical(names(s), names(f))
  expect_identical(names(s$mc), names(f$mc))
  expect_identical(s$dj30, f$dj30)
  expect_identical(s$seed, f$seed)
  expect_lte(s$mc$e1_reps, f$mc$e1_reps)
})

test_that("smoke runs write to results/smoke and never overwrite full results", {
  old <- options(code_a.results = NULL)
  on.exit(options(old))
  code_a_set_mode("smoke")
  expect_equal(code_a_results("tables"), file.path("results", "smoke", "tables"))
  code_a_set_mode("full")
  expect_equal(code_a_results("tables"), file.path("results", "tables"))
})

test_that("Monte Carlo results do not depend on the number of cores", {
  skip_on_os("windows")
  fun <- function(nsim) matrix(rnorm(nsim), ncol = 1)
  a <- mc_run(120, 77, fun, cores = 1)
  b <- mc_run(120, 77, fun, cores = 2)
  expect_identical(a, b)
  expect_equal(nrow(a), 120)
})

test_that("Code A sources contain no absolute paths", {
  files <- c(list.files(file.path(code_a_root, "R"), "\\.R$", full.names = TRUE),
             list.files(file.path(code_a_root, "scripts"), "^0[1-4].*\\.R$", full.names = TRUE))
  hits <- unlist(lapply(files, function(f) grep("(/Users/|/home/|C:\\\\)", readLines(f), value = TRUE)))
  expect_length(hits, 0)
})

# --- empirical regression test (skipped without the downloaded data) -----------------------

test_that("Table 1 point estimates are reproduced on Risk Lab data", {
  p <- file.path(code_a_root, "data", "processed", "dj30_logvol.csv")
  skip_if_not(file.exists(p), "processed Risk Lab data not available")
  panel <- load_logvol_panel("dj30", file.path(code_a_root, "data", "processed"))
  est <- estimate_mfbm(panel_increments(panel, c("AAPL", "ALD", "AMGN", "AXP", "BA"))$X, 1 / 252)
  expect_equal(round(unname(est$H), 4), c(0.2609, 0.2050, 0.1645, 0.2058, 0.2153), tolerance = 1e-4)
  expect_equal(round(est$rho["AAPL", "ALD"], 4), 0.3902)
  expect_equal(round(est$eta["AAPL", "ALD"], 4), 0.0950)
  expect_true(all(est$H > 0 & est$H < 0.75))
})
