# Build log-volatility panels for the DJ30 and Magnificent 7 samples.
# Usage: Rscript scripts/02_prepare_data.R [--refresh-snapshot]
# Outputs: data/processed/{dj30,mag7}_logvol.csv, results/tables/code_a_data_summary.csv
#
# Snapshot protection: the processed panels used for all reported results are
# recorded in data/processed/MD5SUMS (raw files: data/raw_snapshot.csv). Risk Lab
# extends and may revise its history, so a fresh download can change the panels.
# If the rebuilt panels differ from the recorded snapshot the script stops and
# leaves the snapshot untouched; --refresh-snapshot replaces it deliberately.

source("R/code_a_config.R")
source_code_a()
cfg <- code_a_config("full")
refresh <- "--refresh-snapshot" %in% commandArgs(trailingOnly = TRUE)
proc_dir <- file.path("data", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

code_a_log("02_prepare_data: measure = ", cfg$rv_measure)

# raw files vs the recorded download snapshot (information only)
raw_rec_path <- file.path("data", "raw_snapshot.csv")
if (file.exists(raw_rec_path)) {
  raw_rec <- utils::read.csv(raw_rec_path, stringsAsFactors = FALSE)
  now <- code_a_md5(file.path("data", "raw", raw_rec$file))
  changed <- raw_rec$file[is.na(now[file.path("data", "raw", raw_rec$file)]) |
                            now[file.path("data", "raw", raw_rec$file)] != raw_rec$md5]
  if (length(changed)) code_a_log("  raw files differ from data/raw_snapshot.csv: ", paste(changed, collapse = ", "))
}

dj30 <- build_logvol_panel(names(cfg$dj30), cfg$rv_measure, cfg$dj30_start, cfg$dj30_end)
mag7 <- build_logvol_panel(names(cfg$mag7), cfg$rv_measure, cfg$mag7_start, cfg$mag7_end)

tmp <- tempfile("panels")
dir.create(tmp)
utils::write.csv(dj30, file.path(tmp, "dj30_logvol.csv"), row.names = FALSE)
utils::write.csv(mag7, file.path(tmp, "mag7_logvol.csv"), row.names = FALSE)
new_md5 <- unname(tools::md5sum(file.path(tmp, c("dj30_logvol.csv", "mag7_logvol.csv"))))
sums_path <- file.path(proc_dir, "MD5SUMS")
if (file.exists(sums_path) && !refresh) {
  rec <- utils::read.table(sums_path, col.names = c("md5", "file"), stringsAsFactors = FALSE)
  old_md5 <- rec$md5[match(c("dj30_logvol.csv", "mag7_logvol.csv"), rec$file)]
  if (!identical(old_md5, new_md5)) {
    stop("rebuilt panels differ from the recorded snapshot (data/processed/MD5SUMS); ",
         "the snapshot was left unchanged. Rerun with --refresh-snapshot to replace it deliberately.")
  }
  code_a_log("  rebuilt panels identical to the recorded snapshot")
}
file.copy(file.path(tmp, c("dj30_logvol.csv", "mag7_logvol.csv")), proc_dir, overwrite = TRUE)
writeLines(paste(new_md5, c("dj30_logvol.csv", "mag7_logvol.csv")), sums_path)
code_a_log(sprintf("DJ30 panel: %d dates, %s to %s", nrow(dj30), min(dj30$date), max(dj30$date)))
code_a_log(sprintf("Mag7 panel: %d dates, %s to %s", nrow(mag7), min(mag7$date), max(mag7$date)))

summ <- rbind(summarise_panel(dj30, "DJ30"), summarise_panel(mag7, "Mag7"))
code_a_write(summ, "code_a_data_summary.csv")

# Fail loudly on obviously broken inputs.
stopifnot(all(summ$n_obs > 500), all(is.finite(summ$mean_logvol)))
code_a_log_run("02_prepare_data", "full", cfg)
code_a_log("02_prepare_data: done")
