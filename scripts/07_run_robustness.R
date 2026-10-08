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
specs <- rbind(
  data.frame(table = "Table 17", panel = names(classes), class = names(classes),
             period = "full", metric = "MSFE"),
  data.frame(table = "Table 19", panel = names(classes), class = names(classes),
             period = "full", metric = "QLIKE")
)
dir.create(code_b_results("forecasts"), recursive = TRUE, showWarnings = FALSE)
fcs <- tabs_all <- list()
for (cal in cfg$calendars) { # both time conventions of review B1, primary first
  cfg_cal <- code_b_calendar_config(cfg, cal)
  sfx <- code_b_calendar_suffix(cfg, cal)
  code_b_log("calendar = ", cal, if (sfx == "") " (primary)" else " (sensitivity / replication mode)")
  fc_mag7 <- run_all_models(mag7, cfg$mag7_assets, cfg_cal, cores, checkpoint = paste0("mag7_", cal))
  saveRDS(fc_mag7, code_b_results("forecasts", paste0("code_b_forecasts_mag7", sfx, ".rds")))
  code_b_write(cbind(calendar = cal, check_forecasts(fc_mag7)), paste0("code_b_checks_mag7", sfx, ".csv"))
  tabs <- paper_tables(fc_mag7, specs, cfg$mag7_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
  tabs$calendar <- cal
  cmp <- compare_with_paper(tabs, paper)
  code_b_write(cmp, paste0("code_b_mag7_tables_vs_paper", sfx, ".csv"))
  ovo <- original_vs_ours(cmp)
  code_b_write(ovo, paste0("code_b_mag7_original_vs_ours", sfx, ".csv"))
  if (sfx == "") plot_class_ratio(ovo[ovo$table == "Table 17", ], code_b_fig("code_b_fig_mag7_class_ratio.pdf"),
                                  "Mag7: AAPL + AMZN, FB, GOOG, MSFT")
  code_b_log("[", cal, "] mag7 same direction as the paper: ", sum(ovo$same_direction), " of ", nrow(ovo))
  fcs[[cal]] <- fc_mag7
  tabs_all[[cal]] <- tabs
}
if (all(c("trading", "common_obs") %in% names(fcs))) {
  code_b_write(calendar_sensitivity(tabs_all, paper), "code_b_calendar_sensitivity_mag7.csv")
  code_b_write(calendar_key_audit(fcs$trading, fcs$common_obs, mag7$date), "code_b_calendar_keys_mag7.csv")
}
code_b_write(code_b_provenance("07_run_robustness.R",
  c("data/processed/mag7_logvol.csv", "data/processed/dj30_logvol.csv"), cfg), "code_b_provenance_07.csv")

# --- 2. DJ30 HAR classes in the subperiods --------------------------------------

# sections 2 and 3 use the primary convention (forecasts of script 06)
fc_dj30 <- readRDS(code_b_results("forecasts", "code_b_forecasts_dj30.rds"))
stopifnot(identical(unique(fc_dj30$calendar), cfg$calendar))
specs <- expand.grid(class = names(classes), period = names(cfg$dj30_periods),
                     stringsAsFactors = FALSE)
specs$table <- "robust_subperiods"
specs$panel <- paste(specs$class, specs$period)
specs$metric <- "MSFE"
sub <- paper_tables(fc_dj30, specs, cfg$dj30_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
code_b_write(sub, "code_b_robust_dj30_subperiods.csv")

# --- 3. DJ30 with a 250-day window ----------------------------------------------

cfg250 <- code_b_calendar_config(cfg, cfg$calendar)
cfg250$window <- 250L
cfg250$origin_dates <- sort(unique(fc_dj30$origin_date))
dj30 <- load_logvol_panel("dj30")
fc250 <- run_all_models(dj30, cfg$dj30_assets, cfg250, cores,
                        model_names = c(classes$mfbm, classes$vhar), checkpoint = "dj30_window250")
code_b_write(check_forecasts(fc250), "code_b_checks_dj30_window250.csv")
# Match exact keys too: the start-of-period counts can differ at h > 1.
common <- align_forecast_keys(
  fc_dj30[fc_dj30$model %in% c(classes$mfbm, classes$vhar), ], fc250)
fc250 <- common$comparison
specs <- data.frame(table = "robust_window250", panel = c("mfbm", "vhar"),
                    class = c("mfbm", "vhar"), period = "full", metric = "MSFE")
w250 <- paper_tables(fc250, specs, cfg$dj30_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
code_b_write(w250, "code_b_robust_dj30_window250.csv")
w500 <- paper_tables(common$reference, specs, cfg$dj30_periods, cfg$mcs_boot, cfg$mcs_block, cfg$seed)
code_b_write(w500, "code_b_robust_dj30_window500_matched.csv")
matched <- merge(w250[, c("class","model","h","value","n_forecasts")],
  w500[, c("class","model","h","value","n_forecasts")],
  by = c("class","model","h"), suffixes = c("_250","_500"))
stopifnot(all(matched$n_forecasts_250 == matched$n_forecasts_500))
matched$ratio_250_to_500 <- matched$value_250 / matched$value_500
code_b_write(matched, "code_b_robust_window_comparison_matched.csv")

code_b_log("07 done")
