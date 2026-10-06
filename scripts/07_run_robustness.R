# code B step 7: robustness.
#   1. magnificent 7 (appendix F.3, F.5): AAPL + AMZN, FB, GOOG, MSFT, tables 17 and 19
#   2. DJ30 HAR classes in the two subperiods of table 2 (not in the paper, uses the
#      forecasts saved by script 06)
#   3. DJ30 with a one-year window (250 days) for the mfBm and VHAR classes
#
# run from the project root, after script 06:
#   Rscript scripts/07_run_robustness.R smoke
#   Rscript scripts/07_run_robustness.R full

source("R/code_b_config.R")
code_b_relaunch_single_thread("scripts/07_run_robustness.R")
source_code_b()
mode <- code_b_mode()
cfg <- code_b_config(mode)
code_b_set_mode(mode)
cores <- code_b_cores()
code_b_log("07 robustness, mode = ", mode, ", cores = ", cores)
paper <- utils::read.csv("docs/paper_values/code_b_forecast_tables.csv")
classes <- code_b_classes()

# --- 1. magnificent 7 ----------------------------------------------------------

mag7 <- load_logvol_panel("mag7")
fc_mag7 <- run_all_models(mag7, cfg$mag7_assets, cfg, cores, checkpoint = "mag7")
dir.create(code_b_results("forecasts"), recursive = TRUE, showWarnings = FALSE)
saveRDS(fc_mag7, code_b_results("forecasts", "code_b_forecasts_mag7.rds"))
code_b_write(check_forecasts(fc_mag7), "code_b_checks_mag7.csv")

specs <- rbind(
  data.frame(table = "Table 17", panel = names(classes), class = names(classes),
             period = "full", metric = "MSFE"),
  data.frame(table = "Table 19", panel = names(classes), class = names(classes),
             period = "full", metric = "QLIKE")
)
tabs <- paper_tables(fc_mag7, specs, cfg$mag7_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
cmp <- compare_with_paper(tabs, paper)
code_b_write(cmp, "code_b_mag7_tables_vs_paper.csv")
ovo <- original_vs_ours(cmp)
code_b_write(ovo, "code_b_mag7_original_vs_ours.csv")
plot_class_ratio(ovo[ovo$table == "Table 17", ], code_b_fig("code_b_fig_mag7_class_ratio.pdf"),
                 "Mag7: AAPL + AMZN, FB, GOOG, MSFT")
code_b_log("mag7 same direction as the paper: ", sum(ovo$same_direction), " of ", nrow(ovo))

# --- 2. DJ30 HAR classes in the subperiods --------------------------------------

fc_dj30 <- readRDS(code_b_results("forecasts", "code_b_forecasts_dj30.rds"))
specs <- expand.grid(class = names(classes), period = names(cfg$dj30_periods),
                     stringsAsFactors = FALSE)
specs$table <- "robust_subperiods"
specs$panel <- paste(specs$class, specs$period)
specs$metric <- "MSFE"
sub <- paper_tables(fc_dj30, specs, cfg$dj30_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
code_b_write(sub, "code_b_robust_dj30_subperiods.csv")

# --- 3. DJ30 with a 250-day window ----------------------------------------------

cfg250 <- cfg
cfg250$window <- 250L
dj30 <- load_logvol_panel("dj30")
fc250 <- run_all_models(dj30, cfg$dj30_assets, cfg250, cores,
                        model_names = c(classes$mfbm, classes$vhar), checkpoint = "dj30_window250")
code_b_write(check_forecasts(fc250), "code_b_checks_dj30_window250.csv")
# same evaluation period as table 2 so the numbers are comparable with 500 days
specs <- data.frame(table = "robust_window250", panel = c("mfbm", "vhar"),
                    class = c("mfbm", "vhar"), period = "full", metric = "MSFE")
w250 <- paper_tables(fc250, specs, cfg$dj30_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
code_b_write(w250, "code_b_robust_dj30_window250.csv")

code_b_log("07 done")
