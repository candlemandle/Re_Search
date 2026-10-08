# Monte Carlo experiments of Code A (Appendix E.1, E.2, E.4 and G).
# Every experiment takes its sizes from code_a_config(); smoke and full runs
# execute identical code. Replications are simulated in fixed-size chunks with
# a seed per chunk, so results do not depend on the number of cores.

MC_CHUNK <- 50L

# Run `fun(nsim)` (returning a matrix with nsim rows) over `reps` replications.
mc_run <- function(reps, seed, fun, cores = code_a_cores()) {
  n_chunks <- ceiling(reps / MC_CHUNK)
  sizes <- rep(MC_CHUNK, n_chunks)
  sizes[n_chunks] <- reps - MC_CHUNK * (n_chunks - 1)
  job <- function(k) {
    set.seed(seed + 7919L * k)
    fun(sizes[k])
  }
  res <- if (cores > 1) {
    parallel::mclapply(seq_len(n_chunks), job, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(n_chunks), job)
  }
  err <- vapply(res, inherits, TRUE, what = "try-error")
  if (any(err)) stop(res[[which(err)[1]]])
  do.call(rbind, res)
}

# Scenario-specific seed so that adding scenarios does not change others.
mc_seed <- function(base, ...) {
  key <- paste(..., sep = "|")
  base + sum(utf8ToInt(key) * seq_along(utf8ToInt(key))) %% 100000L
}

summarise_mc <- function(est, truth) {
  data.frame(
    parameter = colnames(est),
    true = truth[colnames(est)],
    mean = colMeans(est),
    bias = colMeans(est) - truth[colnames(est)],
    sd = apply(est, 2, stats::sd),
    rmse = sqrt(colMeans(sweep(est, 2, truth[colnames(est)])^2)),
    row.names = NULL
  )
}

# --- Tables 4-5: finite-sample performance of the MM estimators ---------------
mc_estimators <- function(cfg, cores = code_a_cores()) {
  mc <- cfg$mc
  H <- c(0.1, 0.4)
  out <- list()
  for (rho in mc$e1_rho) for (n in mc$e1_n) for (delta in mc$e1_delta) {
    emb <- mfbm_embedding(mfbm_params(H, rho, 0, c(1, 1)), n, delta)
    seed <- mc_seed(cfg$seed, "E1", rho, n, delta)
    est <- mc_run(mc$e1_reps, seed, function(nsim) {
      X <- mfbm_simulate_increments(emb, nsim)
      t(vapply(seq_len(nsim), function(s) {
        e <- estimate_mfbm(X[, , s], delta)
        c(H1 = e$H[1], H2 = e$H[2], sigma2_1 = e$sigma2[1], sigma2_2 = e$sigma2[2],
          rho = e$rho[1, 2], eta = e$eta[1, 2])
      }, numeric(6)))
    }, cores)
    truth <- c(H1 = H[1], H2 = H[2], sigma2_1 = 1, sigma2_2 = 1, rho = rho, eta = 0)
    s <- summarise_mc(est, truth)
    s$asym_se <- c(sqrt(avar_hurst(H[1], cfg$series_terms) / n),
                   sqrt(avar_hurst(H[2], cfg$series_terms) / n),
                   se_sigma2(H[1], 1, n, delta, cfg$series_terms),
                   se_sigma2(H[2], 1, n, delta, cfg$series_terms),
                   sqrt(avar_rho(H[1], H[2], rho, cfg$series_terms) / n),
                   sqrt(avar_eta(H[1], H[2], rho, cfg$series_terms) / n))
    out[[length(out) + 1]] <- cbind(
      table = if (rho == 0) "Table 4" else "Table 5", method = "BYZ",
      n = n, delta = sprintf("1/%d", round(1 / delta)), rho = rho, eta = 0,
      reps = mc$e1_reps, s)
  }
  do.call(rbind, out)
}

# --- Table 6: BYZ vs Amblard-Coeurjolly ----------------------------------------
mc_ac_comparison <- function(cfg, cores = code_a_cores()) {
  mc <- cfg$mc
  H <- c(0.1, 0.4)
  delta <- 1 / 250
  out <- list()
  for (g in seq_len(nrow(mc$t6_grid))) for (n in mc$t6_n) {
    rho <- mc$t6_grid$rho[g]
    eta <- mc$t6_grid$eta[g]
    emb <- mfbm_embedding(mfbm_params(H, rho, eta, c(1, 1)), n, delta)
    seed <- mc_seed(cfg$seed, "T6", rho, eta, n)
    est <- mc_run(mc$t6_reps, seed, function(nsim) {
      X <- mfbm_simulate_increments(emb, nsim)
      t(vapply(seq_len(nsim), function(s) {
        x <- X[, , s]
        ac <- estimate_ac(x, delta, mc$ac_dilations)
        c(BYZ_rho = rho_mm(x[, 1], x[, 2]), BYZ_eta = eta_mm(x[, 1], x[, 2]),
          AC_rho = ac$rho[1, 2], AC_eta = ac$eta[1, 2])
      }, numeric(4)))
    }, cores)
    s <- summarise_mc(est, c(BYZ_rho = rho, BYZ_eta = eta, AC_rho = rho, AC_eta = eta))
    s$method <- sub("_.*", "", s$parameter)
    s$parameter <- sub(".*_", "", s$parameter)
    out[[length(out) + 1]] <- cbind(table = "Table 6", n = n, delta = "1/250", rho = rho, eta = eta,
                                    reps = mc$t6_reps, s)
  }
  do.call(rbind, out)
}

# --- Table 7: size and power of the time-reversibility test ---------------------
mc_test_power <- function(cfg, cores = code_a_cores()) {
  mc <- cfg$mc
  H <- c(0.1, 0.4)
  rho <- 0.4
  delta <- 1 / 250
  R <- min(cfg$series_terms, 5000L)
  out <- list()
  for (eta in mc$t7_eta) for (n in mc$t7_n) {
    emb <- mfbm_embedding(mfbm_params(H, rho, eta, c(1, 1)), n, delta)
    seed <- mc_seed(cfg$seed, "T7", eta, n)
    res <- mc_run(mc$t7_reps, seed, function(nsim) {
      X <- mfbm_simulate_increments(emb, nsim)
      t(vapply(seq_len(nsim), function(s) {
        tt <- time_reversibility_test(X[, 1, s], X[, 2, s], R = R)
        c(stat = tt$statistic, p = tt$p_value)
      }, numeric(2)))
    }, cores)
    # rates over replications where the test is defined (H_hat in (0, 3/4));
    # the share of undefined cases is reported separately
    ok <- !is.na(res[, "p"])
    r05 <- mean(res[ok, "p"] < 0.05)
    out[[length(out) + 1]] <- data.frame(
      table = "Table 7", n = n, delta = "1/250", rho = rho, eta = eta, reps = mc$t7_reps,
      rejection_rate_0.01 = mean(res[ok, "p"] < 0.01),
      rejection_rate_0.05 = r05,
      mc_se_0.05 = sqrt(r05 * (1 - r05) / sum(ok)),
      share_invalid = mean(!ok),
      # alternative convention: undefined cases counted as non-rejections
      rejection_rate_0.01_all = sum(res[ok, "p"] < 0.01) / nrow(res),
      rejection_rate_0.05_all = sum(res[ok, "p"] < 0.05) / nrow(res),
      approximate_simulation = isTRUE(emb$approximate))
  }
  do.call(rbind, out)
}

# --- Table 12: measurement error in log RV (Appendix E.4) -----------------------
# Observed log RV = latent log IV (bivariate fBm) + i.i.d. N(0, 2 / m) noise,
# m = number of intraday returns (Fukasawa et al. 2022, Thm 2.1).
mc_measurement_error <- function(cfg, cores = code_a_cores()) {
  mc <- cfg$mc
  H <- c(0.1, 0.4)
  delta <- 1 / 250
  noise_sd <- sqrt(2 / mc$t12_intraday)
  grid <- expand.grid(rho = c(0, 0.4), eta = c(0, 0.5))
  out <- list()
  for (g in seq_len(nrow(grid))) for (n in mc$t12_n) {
    rho <- grid$rho[g]; eta <- grid$eta[g]
    emb <- mfbm_embedding(mfbm_params(H, rho, eta, c(1, 1)), n, delta)
    seed <- mc_seed(cfg$seed, "T12", rho, eta, n)
    est <- mc_run(mc$t12_reps, seed, function(nsim) {
      X <- mfbm_simulate_increments(emb, nsim)
      t(vapply(seq_len(nsim), function(s) {
        x <- X[, , s]
        e <- matrix(stats::rnorm(2 * (n + 1), sd = noise_sd), n + 1, 2)
        y <- x + apply(e, 2, diff)
        c(bfBm_rho = rho_mm(x[, 1], x[, 2]), bfBm_eta = eta_mm(x[, 1], x[, 2]),
          noise_rho = rho_mm(y[, 1], y[, 2]), noise_eta = eta_mm(y[, 1], y[, 2]))
      }, numeric(4)))
    }, cores)
    s <- summarise_mc(est, c(bfBm_rho = rho, bfBm_eta = eta, noise_rho = rho, noise_eta = eta))
    s$method <- ifelse(grepl("^noise", s$parameter), "bfBm+noise", "bfBm")
    s$parameter <- sub(".*_", "", s$parameter)
    out[[length(out) + 1]] <- cbind(table = "Table 12", n = n, delta = "1/250", rho = rho, eta = eta,
                                    reps = mc$t12_reps, s)
  }
  do.call(rbind, out)
}

# --- Appendix G: lagged / second-order / optimal GMM correlation estimators ----
mc_gmm <- function(cfg, cores = code_a_cores()) {
  mc <- cfg$mc
  H <- c(0.1, 0.4)
  rho <- 0.4
  n <- mc$g_n
  L <- mc$g_lags
  delta <- 1 / n
  R <- min(cfg$series_terms, 2000L)
  emb <- mfbm_embedding(mfbm_params(H, rho, 0, c(1, 1)), n, delta)
  seed <- mc_seed(cfg$seed, "G", n, L)
  est <- mc_run(mc$g_reps, seed, function(nsim) {
    X <- mfbm_simulate_increments(emb, nsim)
    t(vapply(seq_len(nsim), function(s) {
      x1 <- X[, 1, s]; x2 <- X[, 2, s]
      opt <- rho_opt(x1, x2, L, R)
      boot <- if (mc$g_boot > 0) rho_opt_boot(x1, x2, L, mc$g_boot, delta)$estimate else NA_real_
      c(opt$rho_lags, rho_second_order(x1, x2), opt$estimate, boot)
    }, numeric(L + 3)))
  }, cores)
  colnames(est) <- c(sprintf("rho_lag%d", seq_len(L)), "rho_2nd_order", "rho_opt", "rho_opt_boot")
  # The 1000-replication variance estimates have ~5% relative MC error; the cheap
  # estimators are re-estimated on an independent, larger set of paths.
  cheap <- mc_run(mc$g_reps_cheap, mc_seed(cfg$seed, "G-cheap", n, L), function(nsim) {
    X <- mfbm_simulate_increments(emb, nsim)
    t(vapply(seq_len(nsim), function(s) {
      c(rho_lags(X[, 1, s], X[, 2, s], L), rho_second_order(X[, 1, s], X[, 2, s]))
    }, numeric(L + 1)))
  }, cores)
  C_true <- gmm_rho_cov(H[1], H[2], rho, L, max(R, 5000L))
  list(
    estimates = est,
    summary = data.frame(
      method = colnames(est), mean = colMeans(est),
      n_variance = n * apply(est, 2, stats::var),
      n_rmse = sqrt(n * colMeans((est - rho)^2)),
      n_variance_large_mc = c(n * apply(cheap, 2, stats::var), NA, NA),
      reps = mc$g_reps,
      reps_large_mc = c(rep(mc$g_reps_cheap, L + 1), NA, NA),
      asymptotic_n_variance = c(diag(C_true), NA, 1 / sum(solve(C_true, rep(1, L))), NA),
      row.names = NULL),
    weights_true = gmm_weights(C_true)
  )
}

# ---------------------------------------------------------------------------
# Cell-level agreement between our Monte Carlo tables and the paper
# ---------------------------------------------------------------------------
# For each cell: z = (ours - paper) / sqrt(se_ours^2 + se_paper^2), where both
# Monte Carlo standard errors use the respective number of replications and
# plug-in moments (paper sd when the paper reports it, otherwise ours):
#   bias / mean     sd / sqrt(R)
#   sd              sd / sqrt(2 (R - 1))                (normal approximation)
#   rmse            sqrt(2 sd^4 + 4 bias^2 sd^2) / (2 rmse sqrt(R))  (delta method)
#   n * variance    n var sqrt(2 / (R - 1));  sqrt(n * MSE) analogous to rmse
#   rejection rate  sqrt(p (1 - p) / R)
# Asymptotic standard errors are deterministic: agreement means equal after
# rounding to the paper's 4 decimals. Cells share simulated paths (correlated)
# and there are many of them, so about 5% of |z| > 2 is expected by chance;
# the criterion is a screening device, not a joint test.
PAPER_REPS <- c("Table 4" = 1000, "Table 5" = 1000, "Table 6" = 1000, "Table 7" = 5000,
                "Table 12" = 1000, "Table 20" = 1000)

mc_agreement <- function(tables_dir, ref_path = file.path("docs", "paper_values", "monte_carlo_reference.csv")) {
  ref <- utils::read.csv(ref_path, stringsAsFactors = FALSE)
  rd <- function(f) {
    p <- file.path(tables_dir, f)
    if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
  }
  key <- function(d) paste(d$table, d$method, d$n, d$delta, d$rho, d$eta, d$parameter)
  out <- list()
  se_rmse <- function(sd, bias, rmse, R) sqrt(2 * sd^4 + 4 * bias^2 * sd^2) / (2 * rmse * sqrt(R))

  wide <- rbind(rd("code_a_mc_table04_05_estimators.csv")[, c("table", "method", "n", "delta", "rho", "eta",
                                                              "reps", "parameter", "bias", "sd", "rmse", "asym_se")],
                {
                  w <- rbind(rd("code_a_mc_table06_byz_vs_ac.csv"), rd("code_a_mc_table12_measurement_error.csv"))
                  if (!is.null(w)) {
                    w$asym_se <- NA
                    w[, c("table", "method", "n", "delta", "rho", "eta", "reps", "parameter", "bias", "sd", "rmse", "asym_se")]
                  }
                })
  if (!is.null(wide)) {
    for (m in c("bias", "sd", "rmse", "asym_se")) {
      r <- ref[ref$metric == m, ]
      w <- wide[!is.na(wide[[m]]), ]
      idx <- match(key(r), key(w))
      r <- r[!is.na(idx), ]; w <- w[idx[!is.na(idx)], ]
      if (nrow(r) == 0) next
      Rp <- PAPER_REPS[r$table]
      # paper sd for the same cell when published (Tables 4-6), else ours
      psd <- ref$paper_value[match(paste(key(r), "sd"), paste(key(ref), ref$metric))]
      psd[is.na(psd)] <- w$sd[is.na(psd)]
      se_o <- switch(m,
        bias = w$sd / sqrt(w$reps),
        sd = w$sd / sqrt(2 * (w$reps - 1)),
        rmse = se_rmse(w$sd, w$bias, w$rmse, w$reps),
        asym_se = rep(0, nrow(w)))
      se_p <- switch(m,
        bias = psd / sqrt(Rp),
        sd = psd / sqrt(2 * (Rp - 1)),
        rmse = se_rmse(psd, w$bias, w$rmse, Rp),
        asym_se = rep(0, nrow(w)))
      out[[length(out) + 1]] <- data.frame(r[, c("table", "method", "n", "delta", "rho", "eta", "parameter", "metric")],
                                           ours = w[[m]], paper = r$paper_value, reps_ours = w$reps, reps_paper = Rp,
                                           se_ours = se_o, se_paper = se_p)
    }
  }

  t7 <- rd("code_a_mc_table07_test_size_power.csv")
  if (!is.null(t7)) for (a in c("0.01", "0.05")) {
    r <- ref[ref$metric == paste0("rejection_rate_", a), ]
    k7 <- paste(t7$n, t7$eta); idx <- match(paste(r$n, r$eta), k7)
    r <- r[!is.na(idx), ]; w <- t7[idx[!is.na(idx)], ]
    if (nrow(r) == 0) next
    valid_reps <- round(w$reps * (1 - w$share_invalid))
    p_o <- w[[paste0("rejection_rate_", a)]]
    out[[length(out) + 1]] <- data.frame(r[, c("table", "method", "n", "delta", "rho", "eta", "parameter", "metric")],
                                         ours = p_o, paper = r$paper_value, reps_ours = valid_reps, reps_paper = 5000,
                                         se_ours = sqrt(p_o * (1 - p_o) / valid_reps),
                                         se_paper = sqrt(r$paper_value * (1 - r$paper_value) / 5000))
  }

  t20 <- rd("code_a_mc_table20_gmm.csv")
  if (!is.null(t20) && any(ref$table == "Table 20")) {
    t20$method[t20$method == "rho_lag1"] <- "rho_1"
    r <- ref[ref$table == "Table 20" & ref$method %in% t20$method, ]
    w <- t20[match(r$method, t20$method), ]
    n <- r$n; R <- w$reps; v <- w$n_variance / n; bias <- w$mean - r$rho
    nm <- r$metric
    ours <- ifelse(nm == "mean", w$mean, ifelse(nm == "n_variance", w$n_variance, w$n_rmse))
    se_fun <- function(RR) ifelse(nm == "mean", sqrt(v / RR),
                           ifelse(nm == "n_variance", w$n_variance * sqrt(2 / (RR - 1)),
                                  sqrt(n) * se_rmse(sqrt(v), bias, w$n_rmse / sqrt(n), RR)))
    out[[length(out) + 1]] <- data.frame(r[, c("table", "method", "n", "delta", "rho", "eta", "parameter", "metric")],
                                         ours = ours, paper = r$paper_value, reps_ours = R, reps_paper = 1000,
                                         se_ours = se_fun(R), se_paper = se_fun(1000))
  }

  res <- do.call(rbind, out)
  res$diff <- res$ours - res$paper
  res$se_diff <- sqrt(res$se_ours^2 + res$se_paper^2)
  det <- res$metric == "asym_se"
  res$z <- ifelse(det, NA_real_, ifelse(res$se_diff > 0, res$diff / res$se_diff, ifelse(res$diff == 0, 0, Inf)))
  # paper values are rounded to 4 decimals (Table 20: 2-3); rounding tolerance
  tol <- ifelse(res$table == "Table 20", 0.005, 5e-5)
  res$agree <- ifelse(det, abs(round(res$ours, 4) - res$paper) < 1e-9,
                      abs(res$z) <= 2 | abs(res$diff) <= tol)
  res$criterion <- ifelse(det, "equal after rounding to 4 decimals",
                          "|z| <= 2 (or |diff| within paper rounding)")
  res
}

mc_agreement_summary <- function(a) {
  do.call(rbind, lapply(split(a, paste(a$table, a$metric == "asym_se")), function(d) {
    data.frame(table = d$table[1],
               cells = if (d$metric[1] == "asym_se") "asymptotic SE (deterministic)" else "Monte Carlo",
               n_cells = nrow(d), n_agree = sum(d$agree), share_agree = mean(d$agree),
               n_abs_z_gt_3 = sum(abs(d$z) > 3, na.rm = TRUE),
               max_abs_z = if (all(is.na(d$z))) NA_real_ else max(abs(d$z), na.rm = TRUE))
  }))
}
