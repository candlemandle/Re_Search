# Build log-volatility panels for the DJ30 and Magnificent 7 samples.
# Usage: Rscript scripts/02_prepare_data.R
# Outputs: data/processed/{dj30,mag7}_logvol.csv, results/tables/code_a_data_summary.csv

source("R/code_a_config.R")
source_code_a()
cfg <- code_a_config("full")
dir.create(file.path("data", "processed"), recursive = TRUE, showWarnings = FALSE)

code_a_log("02_prepare_data: measure = ", cfg$rv_measure)
dj30 <- build_logvol_panel(names(cfg$dj30), cfg$rv_measure, cfg$dj30_start, cfg$dj30_end)
mag7 <- build_logvol_panel(names(cfg$mag7), cfg$rv_measure, cfg$mag7_start, cfg$mag7_end)

utils::write.csv(dj30, file.path("data", "processed", "dj30_logvol.csv"), row.names = FALSE)
utils::write.csv(mag7, file.path("data", "processed", "mag7_logvol.csv"), row.names = FALSE)
code_a_log(sprintf("DJ30 panel: %d dates, %s to %s", nrow(dj30), min(dj30$date), max(dj30$date)))
code_a_log(sprintf("Mag7 panel: %d dates, %s to %s", nrow(mag7), min(mag7$date), max(mag7$date)))

summ <- rbind(summarise_panel(dj30, "DJ30"), summarise_panel(mag7, "Mag7"))
code_a_write(summ, "code_a_data_summary.csv")

# Fail loudly on obviously broken inputs.
stopifnot(all(summ$n_obs > 500), all(is.finite(summ$mean_logvol)))
code_a_log("02_prepare_data: done")
