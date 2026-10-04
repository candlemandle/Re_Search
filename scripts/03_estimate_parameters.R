# Parameter estimation on Risk Lab data and model-structure figures.
# Usage: Rscript scripts/03_estimate_parameters.R [smoke|full]
#
# Outputs (results/tables/code_a_*.csv, results/figures/code_a_*.pdf; smoke mode writes to
# results/smoke/ instead so it never overwrites full results):
#   Table 1      code_a_table01_estimates.csv     (5 stocks, common sample, with SEs)
#   Tables 13-14 code_a_table13_14_dj30.csv       (30 stocks: H, rho, eta, SEs)
#   Test (11)    code_a_time_reversibility_{dj30,mag7}.csv, code_a_time_reversibility_summary.csv
#   Mag7         code_a_mag7_estimates.csv
#   Figures 1-3  code_a_fig01_sample_paths.pdf, code_a_fig02_rho_max.pdf, code_a_fig03_admissible.pdf
#   Figures 7-8  code_a_fig07_rolling_H_dj30.pdf, code_a_fig08_rolling_H_mag7.pdf (+ csv)
#   Comparison   code_a_original_vs_ours_empirical.csv

source("R/code_a_config.R")
source_code_a()
mode <- code_a_mode()
cfg <- code_a_config(mode)
code_a_set_mode(mode)
R <- cfg$series_terms
code_a_log("03_estimate_parameters: mode = ", mode)

dj30 <- load_logvol_panel("dj30")
mag7 <- load_logvol_panel("mag7")
delta <- cfg$delta_daily

long_estimates <- function(est, se, sample) {
  tk <- names(est$H)
  rows <- list(data.frame(sample = sample, asset1 = tk, asset2 = "", statistic = "H",
                          estimate = est$H, se = se$H, n = est$n),
               data.frame(sample = sample, asset1 = tk, asset2 = "", statistic = "sigma2",
                          estimate = est$sigma2, se = se$sigma2, n = est$n))
  for (i in seq_along(tk)) for (j in seq_along(tk)) if (i < j) {
    rows[[length(rows) + 1]] <- data.frame(
      sample = sample, asset1 = tk[i], asset2 = tk[j], statistic = c("rho", "eta"),
      estimate = c(est$rho[i, j], est$eta[i, j]), se = c(se$rho[i, j], se$eta[i, j]), n = est$n)
  }
  do.call(rbind, rows)
}

# --- Table 1: five-dimensional mfBm on the dates common to all five stocks ---
inc5 <- panel_increments(dj30, cfg$table1_tickers)
est5 <- estimate_mfbm(inc5$X, delta)
se5 <- mfbm_standard_errors(est5, R)
t1 <- long_estimates(est5, se5, "DJ30 Table 1 (common sample of 5)")
code_a_log(sprintf("Table 1: n = %d increments, %s to %s", est5$n, min(inc5$dates), max(inc5$dates)))

# --- Tables 13-14: H per series on own dates, (rho, eta) pairwise -------------
estp <- estimate_panel(dj30, delta, R)
tk <- names(estp$H)
t13 <- list(data.frame(sample = "DJ30 Tables 13-14", asset1 = tk, asset2 = "", statistic = "H",
                       estimate = estp$H, se = estp$se_H, n = estp$n))
for (i in seq_along(tk)) for (j in seq_along(tk)) if (i < j) {
  t13[[length(t13) + 1]] <- data.frame(
    sample = "DJ30 Tables 13-14", asset1 = tk[i], asset2 = tk[j], statistic = c("rho", "eta"),
    estimate = c(estp$rho[i, j], estp$eta[i, j]), se = c(estp$se_rho[i, j], estp$se_eta[i, j]),
    n = estp$n_pair[i, j])
}
t13 <- do.call(rbind, t13)

# --- attach paper values -------------------------------------------------------
ref13 <- utils::read.csv("docs/paper_values/table13_H_rho.csv", stringsAsFactors = FALSE)
ref14 <- utils::read.csv("docs/paper_values/table14_eta.csv", stringsAsFactors = FALSE)
ref13$asset2[is.na(ref13$asset2)] <- ""
key_sym <- function(a, b, s) ifelse(s == "rho", paste(pmin(a, b), pmax(a, b), s), paste(a, b, s))
paper_lookup <- function(df) {
  ref <- rbind(ref13, ref14)
  k_ref <- key_sym(ref$asset1, ref$asset2, ref$statistic)
  k <- key_sym(df$asset1, df$asset2, df$statistic)
  df$paper_value <- ref$paper_value[match(k, k_ref)]
  df$abs_diff <- abs(df$estimate - df$paper_value)
  df
}
t1 <- paper_lookup(t1)
ref_se <- utils::read.csv("docs/paper_values/table01_se.csv", stringsAsFactors = FALSE)
ref_se$asset2[is.na(ref_se$asset2)] <- ""
t1$paper_se <- ref_se$paper_se[match(paste(t1$asset1, t1$asset2, t1$statistic),
                                     paste(ref_se$asset1, ref_se$asset2, ref_se$statistic))]
t13 <- paper_lookup(t13)
code_a_write(t1, "code_a_table01_estimates.csv")
code_a_write(t13, "code_a_table13_14_dj30.csv")

# --- time-reversibility tests, eq. (11) ---------------------------------------
tr_dj30 <- do.call(rbind, lapply(seq_along(tk), function(i) {
  do.call(rbind, lapply(seq_along(tk), function(j) {
    if (i >= j) return(NULL)
    X <- panel_increments(dj30, tk[c(i, j)])$X
    t <- time_reversibility_test(X[, 1], X[, 2], R = R)
    data.frame(asset1 = tk[i], asset2 = tk[j], n = t$n, eta = t$eta, se_eta = t$se_eta,
               statistic = t$statistic, p_value = t$p_value,
               reject_1pct = t$reject_0.01, reject_5pct = t$reject_0.05)
  }))
}))
tk7 <- names(cfg$mag7)
tr_mag7 <- time_reversibility_table(panel_increments(mag7, tk7)$X, R = R)
names(tr_mag7)[names(tr_mag7) == "reject_0.01"] <- "reject_1pct"
names(tr_mag7)[names(tr_mag7) == "reject_0.05"] <- "reject_5pct"
tr_dj30$valid <- NULL; tr_mag7$valid <- NULL
code_a_write(tr_dj30, "code_a_time_reversibility_dj30.csv")
code_a_write(tr_mag7, "code_a_time_reversibility_mag7.csv")
tr_sum <- data.frame(
  sample = c("DJ30 (all 435 pairs)", "DJ30 Table 1 (10 pairs, common sample)", "Mag7 (10 pairs)"),
  pairs = c(nrow(tr_dj30), NA, nrow(tr_mag7)),
  rejected_1pct = c(sum(tr_dj30$reject_1pct), NA, sum(tr_mag7$reject_1pct)),
  rejected_5pct = c(sum(tr_dj30$reject_5pct), NA, sum(tr_mag7$reject_5pct))
)
tr5 <- time_reversibility_table(inc5$X, R = R)
tr_sum$pairs[2] <- nrow(tr5)
tr_sum$rejected_1pct[2] <- sum(tr5$reject_0.01)
tr_sum$rejected_5pct[2] <- sum(tr5$reject_0.05)
tr_sum$share_rejected_1pct <- tr_sum$rejected_1pct / tr_sum$pairs
tr_sum$share_rejected_5pct <- tr_sum$rejected_5pct / tr_sum$pairs
code_a_write(tr_sum, "code_a_time_reversibility_summary.csv")
code_a_log(sprintf("Reversibility: DJ30 rejected at 1%%: %d / %d", sum(tr_dj30$reject_1pct), nrow(tr_dj30)))

# --- Magnificent 7 subset (Appendix F.3), common sample ------------------------
inc7 <- panel_increments(mag7, tk7)
est7 <- estimate_mfbm(inc7$X, delta)
code_a_write(long_estimates(est7, mfbm_standard_errors(est7, R), "Mag7 (common sample)"),
             "code_a_mag7_estimates.csv")

# --- original vs ours (headline numbers) ----------------------------------------
ovo <- rbind(
  data.frame(item = "Table 1", quantity = paste(t1$statistic, t1$asset1, t1$asset2),
             paper = t1$paper_value, ours = t1$estimate),
  data.frame(item = "Table 13-14 (max |diff| H)", quantity = "max over 30 stocks",
             paper = NA, ours = max(t13$abs_diff[t13$statistic == "H"], na.rm = TRUE)),
  data.frame(item = "Table 13-14 (max |diff| rho)", quantity = "max over 435 pairs",
             paper = NA, ours = max(t13$abs_diff[t13$statistic == "rho"], na.rm = TRUE)),
  data.frame(item = "Table 13-14 (max |diff| eta)", quantity = "max over 435 pairs",
             paper = NA, ours = max(t13$abs_diff[t13$statistic == "eta"], na.rm = TRUE))
)
ovo <- ovo[!is.na(ovo$paper) | grepl("max", ovo$item), ]
ovo$abs_diff <- abs(ovo$ours - ovo$paper)
code_a_write(ovo, "code_a_original_vs_ours_empirical.csv")

# --- Figures 7-8: rolling two-year Hurst estimates ------------------------------
plot_rolling <- function(rh, file, main, vlines = NULL) {
  grDevices::pdf(file, width = 8, height = 4.5)
  tks <- unique(rh$ticker)
  cols <- c("#d7191c", "#1a9641", "#2c7bb6", "#fdae61", "#5e3c99")[seq_along(tks)]
  graphics::plot(range(rh$date), range(rh$H), type = "n", xlab = "Year", ylab = "H", main = main)
  for (k in seq_along(tks)) {
    s <- rh[rh$ticker == tks[k], ]
    graphics::lines(s$date, s$H, col = cols[k], lwd = 1, lty = k)
  }
  if (!is.null(vlines)) graphics::abline(v = vlines, lty = 2, col = "grey40")
  graphics::legend("topright", legend = tks, col = cols, lty = seq_along(tks), bty = "n", cex = 0.8)
  grDevices::dev.off()
  code_a_log("wrote ", file)
}
rh7 <- rolling_hurst(dj30, cfg$table1_tickers, cfg$rolling_window, cfg$rolling_step)
rh8 <- rolling_hurst(mag7, tk7, cfg$rolling_window, cfg$rolling_step)
code_a_write(rh7, "code_a_fig07_rolling_H_dj30.csv")
code_a_write(rh8, "code_a_fig08_rolling_H_mag7.csv")
plot_rolling(rh7, code_a_fig("code_a_fig07_rolling_H_dj30.pdf"),
             "Rolling two-year estimates of H, log RV (DJ30 subset)",
             as.Date(c("2017-04-11", "2021-07-30")))
plot_rolling(rh8, code_a_fig("code_a_fig08_rolling_H_mag7.pdf"),
             "Rolling two-year estimates of H, log RV (Magnificent 7 subset)")

# --- Figure 1: sample paths of a bivariate fBm ---------------------------------
set.seed(cfg$seed)
paths <- increments_to_levels(simulate_mfbm(1000, c(0.1, 0.4), 0.8, 0, delta = 1 / 250)[[1]])
tt <- (0:1000) / 250
grDevices::pdf(code_a_fig("code_a_fig01_sample_paths.pdf"), width = 8, height = 5)
graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
graphics::plot(tt, paths[, 1], type = "l", col = "#2c7bb6", xlab = "Time", ylab = "", main = expression(H[1] == 0.1))
graphics::plot(tt, paths[, 2], type = "l", col = "#2c7bb6", xlab = "Time", ylab = "", main = expression(H[2] == 0.4))
grDevices::dev.off()
code_a_log("wrote ", code_a_fig("code_a_fig01_sample_paths.pdf"))

# --- Figure 2: rho_max(H1, H2) ----------------------------------------------------
hg <- seq(0.005, 0.995, length.out = 199)
rm_grid <- outer(hg, hg, rho_max)
grDevices::pdf(code_a_fig("code_a_fig02_rho_max.pdf"), width = 6, height = 5)
graphics::image(hg, hg, rm_grid, col = grDevices::gray.colors(100, start = 0, end = 1),
                xlab = expression(H[1]), ylab = expression(H[2]), zlim = c(0, 1),
                main = expression(rho[max](H[1], H[2])))
graphics::contour(hg, hg, rm_grid, levels = c(0.25, 0.5, 0.75, 0.9), add = TRUE, col = "#d7191c")
grDevices::dev.off()
code_a_log("wrote ", code_a_fig("code_a_fig02_rho_max.pdf"))

# --- Figure 3: admissible (H1, H2) for rho = 0.5, 0.75, 0.9 ----------------------
grDevices::pdf(code_a_fig("code_a_fig03_admissible.pdf"), width = 10, height = 3.6)
graphics::par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (r in c(0.5, 0.75, 0.9)) {
  graphics::image(hg, hg, (rm_grid >= r) * 1, col = c("white", "grey60"), zlim = c(0, 1),
                  xlab = expression(H[1]), ylab = expression(H[2]), main = bquote(rho == .(r)))
}
grDevices::dev.off()
code_a_log("wrote ", code_a_fig("code_a_fig03_admissible.pdf"))
code_a_log("03_estimate_parameters: done")
