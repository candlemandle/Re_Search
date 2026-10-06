# code B step 6: empirical forecasts of AAPL RV with DJ30 stocks (sections 5.2, 5.3, appendix F).
#   table 2  : fBm, bfBm, mfBm3-5, full period and two subperiods (MSFE + MCS)
#   table 3  : HAR and VHAR2-5
#   table 15 : HAR and VHARF2-5 (one common factor)
#   table 16 : LHAR and VLHAR2-5 (log HAR)
#   table 18 : QLIKE of all four classes
# assets are added in the paper order AAPL, ALD, AMGN, AXP, BA.
#
# run from the project root (needs data/processed/dj30_logvol.csv from code A):
#   Rscript scripts/06_run_empirical_forecasts.R smoke   (every 40th day, results/smoke/)
#   Rscript scripts/06_run_empirical_forecasts.R full    (every day as in the paper)

source("R/code_b_config.R")
code_b_relaunch_single_thread("scripts/06_run_empirical_forecasts.R")
source_code_b()
mode <- code_b_mode()
cfg <- code_b_config(mode)
code_b_set_mode(mode)
cores <- code_b_cores()
code_b_log("06 empirical forecasts DJ30, mode = ", mode, ", cores = ", cores)

# --- forecasts -----------------------------------------------------------------

panel <- load_logvol_panel("dj30")
t0 <- Sys.time()
fc <- run_all_models(panel, cfg$dj30_assets, cfg, cores, checkpoint = "dj30")
code_b_log(sprintf("all models: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

dir.create(code_b_results("forecasts"), recursive = TRUE, showWarnings = FALSE)
saveRDS(fc, code_b_results("forecasts", "code_b_forecasts_dj30.rds"))
code_b_write(check_forecasts(fc), "code_b_checks_dj30.csv")

# --- tables --------------------------------------------------------------------

specs <- rbind(
  data.frame(table = "Table 2", panel = c("full", "period1", "period2"), class = "mfbm",
             period = c("full", "period1", "period2"), metric = "MSFE"),
  data.frame(table = c("Table 3", "Table 15", "Table 16"), panel = "full",
             class = c("vhar", "vharf", "vlhar"), period = "full", metric = "MSFE"),
  data.frame(table = "Table 18", panel = c("mfbm", "vhar", "vharf", "vlhar"),
             class = c("mfbm", "vhar", "vharf", "vlhar"), period = "full", metric = "QLIKE")
)
tabs <- paper_tables(fc, specs, cfg$dj30_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)

paper <- utils::read.csv("docs/paper_values/code_b_forecast_tables.csv")
cmp <- compare_with_paper(tabs, paper)
code_b_write(cmp, "code_b_dj30_tables_vs_paper.csv")
ovo <- original_vs_ours(cmp)
code_b_write(ovo, "code_b_dj30_original_vs_ours.csv")

# one wide file per paper table, easy to read and to put on slides
for (tb in unique(tabs$table)) {
  x <- tabs[tabs$table == tb, ]
  wide <- do.call(rbind, lapply(split(x, x$panel), function(p) {
    cbind(panel = p$panel[1], merge(to_wide(p, "value"), to_wide(p, "mcs_p"), sort = FALSE))
  }))
  code_b_write(wide, sprintf("code_b_table%02d_dj30.csv", as.integer(sub("Table ", "", tb))))
}

code_b_log("mfBm5 better than fBm (MSFE, full period): ",
           sum(ovo$ours < 1 & ovo$table == "Table 2" & ovo$panel == "full"), " of 8 horizons")
code_b_log("same direction as the paper: ", sum(ovo$same_direction), " of ", nrow(ovo))

# --- figures -------------------------------------------------------------------

# msfe of every mfBm model relative to fBm, by horizon, for the three periods
pdf(code_b_fig("code_b_fig_dj30_msfe_ratio.pdf"), width = 10, height = 3.5)
par(mfrow = c(1, 3))
cols <- c(bfBm = "#1b9e77", mfBm3 = "#d95f02", mfBm4 = "#7570b3", mfBm5 = "#e7298a")
for (p in c("full", "period1", "period2")) {
  x <- tabs[tabs$table == "Table 2" & tabs$panel == p, ]
  base <- x$value[x$model == "fBm"]
  ratios <- sapply(names(cols), function(m) x$value[x$model == m] / base)
  matplot(cfg$horizons, ratios, type = "b", pch = 20, lty = 1, col = cols,
          xlab = "horizon h (days)", ylab = "MSFE / MSFE fBm", main = p,
          ylim = range(c(ratios, 0.96, 1.01)))
  abline(h = 1, col = "grey")
}
legend("bottomleft", legend = names(cols), col = cols, lty = 1, pch = 20, bty = "n")
dev.off()

# the largest model of every class relative to its univariate model (full period)
x <- ovo[ovo$table %in% c("Table 2", "Table 3", "Table 15", "Table 16") & ovo$panel == "full", ]
plot_class_ratio(x, code_b_fig("code_b_fig_dj30_class_ratio.pdf"), "DJ30: AAPL + ALD, AMGN, AXP, BA")

# forecasts and realised RV of AAPL, h = 1 and h = 20
pdf(code_b_fig("code_b_fig_dj30_aapl_forecasts.pdf"), width = 10, height = 6)
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
for (h in c(1, 20)) {
  a <- fc[fc$h == h & fc$model == "fBm", ]
  b <- fc[fc$h == h & fc$model == "mfBm5", ]
  r <- fc[fc$h == h & fc$model == "HAR", ]
  plot(a$target_date, a$actual, type = "l", col = "grey70", log = "y",
       xlab = "", ylab = "RV (annualised vol)", main = paste0("AAPL, h = ", h))
  lines(a$target_date, a$forecast, col = "black")
  lines(b$target_date, b$forecast, col = "red")
  lines(r$target_date, r$forecast, col = "blue", lty = 3)
  abline(v = as.Date(c("2017-04-11", "2021-07-30")), lty = 2)
  legend("topright", c("realised", "fBm", "mfBm5", "HAR"), col = c("grey70", "black", "red", "blue"),
         lty = c(1, 1, 1, 3), bty = "n", cex = 0.8)
}
dev.off()

code_b_log("06 done")
