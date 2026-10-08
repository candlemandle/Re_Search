# Monte Carlo experiments of Code A.
# Usage: Rscript scripts/04_run_monte_carlo.R [smoke|full] [experiments...]
#   experiments: e1 (Tables 4-5), t6 (Table 6), t7 (Table 7), t12 (Table 12), g (Table 20,
#   Figures 9-10); default: all. Cores: env MFBM_CORES (default detectCores() - 1).
#
# Outputs (smoke mode: under results/smoke/): results/tables/code_a_mc_*.csv (each with paper values and differences),
#          results/figures/code_a_fig09_*.pdf, code_a_fig10_*.pdf, code_a_mc_power_curve.pdf

args <- commandArgs(trailingOnly = TRUE)

# Forked workers + multithreaded BLAS (e.g. Homebrew OpenBLAS) oversubscribe the CPU
# badly (load ~70 on 8 cores). BLAS reads its thread count at load time, so relaunch
# once with single-threaded BLAS.
if (Sys.getenv("OPENBLAS_NUM_THREADS") == "") {
  status <- system2(file.path(R.home("bin"), "Rscript"), c("scripts/04_run_monte_carlo.R", args),
                    env = c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "VECLIB_MAXIMUM_THREADS=1"))
  quit(save = "no", status = status)
}

source("R/code_a_config.R")
source_code_a()
mode <- code_a_mode(args)
cfg <- code_a_config(mode)
code_a_set_mode(mode)
all_exp <- c("e1", "t6", "t7", "t12", "g")
todo <- intersect(args[-1], all_exp)
if (length(todo) == 0) todo <- all_exp
code_a_log(sprintf("04_run_monte_carlo: mode = %s, experiments = %s, cores = %d",
                   mode, paste(todo, collapse = ","), code_a_cores()))

ref <- utils::read.csv("docs/paper_values/monte_carlo_reference.csv", stringsAsFactors = FALSE)

# Join paper values on (table, method, n, delta, rho, eta, parameter, metric).
attach_paper <- function(long) {
  k <- function(d) paste(d$table, d$method, d$n, d$delta, d$rho, d$eta, d$parameter, d$metric)
  long$paper_value <- ref$paper_value[match(k(long), k(ref))]
  long$diff <- long$value - long$paper_value
  long
}
to_long <- function(df, metrics) {
  do.call(rbind, lapply(metrics, function(m) {
    data.frame(df[, setdiff(names(df), c("true", "mean", "bias", "sd", "rmse", "asym_se"))],
               metric = m, value = df[[m]], stringsAsFactors = FALSE)
  }))
}
timed <- function(label, expr) {
  t0 <- Sys.time()
  val <- force(expr)
  code_a_log(sprintf("  %s finished in %.1f s", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  val
}

# --- Tables 4-5 --------------------------------------------------------------
if ("e1" %in% todo) {
  e1 <- timed("Tables 4-5", mc_estimators(cfg))
  code_a_write(e1, "code_a_mc_table04_05_estimators.csv")
  code_a_write(attach_paper(to_long(e1, c("bias", "sd", "asym_se", "rmse"))),
               "code_a_mc_table04_05_vs_paper.csv")
}

# --- Table 6 -------------------------------------------------------------------
if ("t6" %in% todo) {
  t6 <- timed("Table 6", mc_ac_comparison(cfg))
  code_a_write(t6, "code_a_mc_table06_byz_vs_ac.csv")
  code_a_write(attach_paper(to_long(t6, c("bias", "sd", "rmse"))), "code_a_mc_table06_vs_paper.csv")
}

# --- Table 7 -------------------------------------------------------------------
if ("t7" %in% todo) {
  t7 <- timed("Table 7", mc_test_power(cfg))
  code_a_write(t7, "code_a_mc_table07_test_size_power.csv")
  long7 <- rbind(
    data.frame(t7[, c("table", "n", "delta", "rho", "eta")], method = "BYZ", parameter = "eta",
               metric = "rejection_rate_0.01", value = t7$rejection_rate_0.01),
    data.frame(t7[, c("table", "n", "delta", "rho", "eta")], method = "BYZ", parameter = "eta",
               metric = "rejection_rate_0.05", value = t7$rejection_rate_0.05))
  code_a_write(attach_paper(long7), "code_a_mc_table07_vs_paper.csv")

  grDevices::pdf(code_a_fig("code_a_mc_power_curve.pdf"), width = 7, height = 4.5)
  graphics::plot(range(t7$eta), c(0, 1), type = "n", xlab = expression(eta[12]),
                 ylab = "Rejection rate", main = "Size and power of the time-reversibility test")
  cols <- c("#2c7bb6", "#d7191c")
  for (k in seq_along(cfg$mc$t7_n)) {
    s <- t7[t7$n == cfg$mc$t7_n[k], ]
    graphics::lines(s$eta, s$rejection_rate_0.05, col = cols[k], lwd = 2, type = "b", pch = 16)
    graphics::lines(s$eta, s$rejection_rate_0.01, col = cols[k], lwd = 2, lty = 2, type = "b", pch = 1)
  }
  graphics::abline(h = c(0.01, 0.05), col = "grey60", lty = 3)
  graphics::legend("bottomright", bty = "n", cex = 0.85,
                   legend = c(sprintf("n = %d, 5%%", cfg$mc$t7_n), sprintf("n = %d, 1%%", cfg$mc$t7_n)),
                   col = rep(cols, 2), lty = rep(1:2, each = 2), pch = rep(c(16, 1), each = 2))
  grDevices::dev.off()
  code_a_log("wrote ", code_a_fig("code_a_mc_power_curve.pdf"))
}

# --- Table 12 ------------------------------------------------------------------
if ("t12" %in% todo) {
  t12 <- timed("Table 12", mc_measurement_error(cfg))
  code_a_write(t12, "code_a_mc_table12_measurement_error.csv")
  code_a_write(attach_paper(to_long(t12, "bias")), "code_a_mc_table12_vs_paper.csv")
}

# --- Appendix G: Table 20, Figures 9-10 ----------------------------------------------
if ("g" %in% todo) {
  g <- timed("Table 20 / Figures 9-10", mc_gmm(cfg))
  code_a_write(g$summary, "code_a_mc_table20_gmm.csv")
  code_a_write(data.frame(lag = seq_along(g$weights_true), weight = g$weights_true),
               "code_a_mc_gmm_optimal_weights_true.csv")
  s20 <- g$summary[g$summary$method %in% c("rho_lag1", "rho_2nd_order", "rho_opt", "rho_opt_boot"), ]
  s20$method[s20$method == "rho_lag1"] <- "rho_1"
  long20 <- do.call(rbind, lapply(c("mean", "n_variance", "n_rmse"), function(m) {
    data.frame(table = "Table 20", method = s20$method, n = cfg$mc$g_n, delta = sprintf("1/%d", cfg$mc$g_n),
               rho = 0.4, eta = 0, parameter = "rho", metric = m, value = s20[[m]])
  }))
  code_a_write(attach_paper(long20), "code_a_mc_table20_vs_paper.csv")

  est <- g$estimates
  grDevices::pdf(code_a_fig("code_a_fig09_gmm_estimates.pdf"), width = 8, height = 5)
  it <- seq_len(nrow(est))
  graphics::plot(it, est[, "rho_lag1"], pch = 16, cex = 0.4, col = "black",
                 ylim = range(est[, c("rho_lag1", "rho_2nd_order", "rho_opt")]),
                 xlab = "Iteration", ylab = "Estimates", main = "Realized correlation estimates")
  graphics::points(it, est[, "rho_2nd_order"], pch = 16, cex = 0.4, col = "#4daf4a")
  graphics::points(it, est[, "rho_opt"], pch = 16, cex = 0.4, col = "#377eb8")
  graphics::abline(h = 0.4)
  graphics::legend("topright", c(expression(hat(rho)[1]), expression(hat(rho)[1]^(2)), expression(hat(rho)^opt)),
                   col = c("black", "#4daf4a", "#377eb8"), pch = 16, bty = "n")
  grDevices::dev.off()

  L <- cfg$mc$g_lags
  sv <- g$summary
  grDevices::pdf(code_a_fig("code_a_fig10_gmm_variances.pdf"), width = 7, height = 4.5)
  # lag estimators and rho^(2): variances from the larger independent MC run
  v_lag <- sv$n_variance_large_mc[seq_len(L)]
  graphics::plot(seq_len(L), v_lag, type = "b", pch = 16, xlab = "Lag", ylab = "n x variance",
                 ylim = range(c(v_lag, sv$n_variance, sv$asymptotic_n_variance), na.rm = TRUE),
                 main = "Normalized variances of correlation estimators")
  graphics::lines(seq_len(L), sv$asymptotic_n_variance[seq_len(L)], lty = 3, col = "grey40")
  graphics::abline(h = sv$n_variance[sv$method == "rho_opt"], lty = 2)
  graphics::points(1, sv$n_variance_large_mc[sv$method == "rho_2nd_order"], pch = 18, cex = 1.5)
  graphics::legend("topleft", bty = "n", cex = 0.85,
                   legend = c(expression(hat(rho)[lag] ~ "(MC)"), expression(hat(rho)[lag] ~ "(asymptotic)"),
                              expression(hat(rho)^opt), expression(hat(rho)[1]^(2))),
                   lty = c(1, 3, 2, NA), pch = c(16, NA, NA, 18), col = c("black", "grey40", "black", "black"))
  grDevices::dev.off()
  code_a_log("wrote ", code_a_fig("code_a_fig09_gmm_estimates.pdf"), ", ", code_a_fig("code_a_fig10_gmm_variances.pdf"))
}
# cell-level comparison with the paper from the tables just written
a <- mc_agreement(code_a_results("tables"))
a$source_mode <- mode
code_a_write(a, "code_a_mc_agreement.csv")
code_a_write(cbind(mc_agreement_summary(a), source_mode = mode), "code_a_mc_agreement_summary.csv")
code_a_log_run(paste("04_run_monte_carlo", paste(todo, collapse = ",")), mode, cfg)
code_a_log("04_run_monte_carlo: done")
