# B1 probe: does future missingness of a non-target asset move an h-step target?
# Usage (project root): Rscript docs/review_B/b1_probe.R [common_obs|trading|both]
# Writes docs/review_B/b1_probe.csv. Uses the full run_all_models() pipeline.
source("R/code_b_config.R"); source_code_b()
args <- commandArgs(trailingOnly = TRUE)
modes <- if (!length(args) || args[1] == "both") c("common_obs", "trading") else args[1]
supports_calendar <- "calendar" %in% all.names(body(run_all_models)) ||
  exists("code_b_calendar_index", mode = "function")
if (!supports_calendar) modes <- "common_obs" # code before the B1 change

# Same synthetic mfBm panel as tests/testthat/helper-code-c.R::c_panel().
probe_panel <- function(dates) {
  set.seed(802)
  H <- c(.12, .2, .28, .35, .42, .17, .31)
  r <- matrix(.25, 7, 7); diag(r) <- 1
  x <- simulate_mfbm(length(dates) - 1L, H, r, delta = 1 / 252)[[1]]
  L <- log(.2) + increments_to_levels(x) * .18
  colnames(L) <- LETTERS[1:7]
  data.frame(date = dates, L, check.names = FALSE)[, c("date", LETTERS[1:5])]
}
calendars <- list(
  daily = as.Date("2010-01-01") + seq_len(105L),
  # explicit Monday-Friday trading calendar: no weekend arithmetic involved
  weekday = { d <- seq(as.Date("2010-01-04"), by = "day", length.out = 160L)
              utils::head(d[!format(d, "%u") %in% c("6", "7")], 105L) })

old <- setwd(tempdir()); on.exit(setwd(old)) # run_all_models logs under results/
rows <- list()
for (cal in names(calendars)) for (mode in modes) {
  p <- probe_panel(calendars[[cal]])
  t <- 70L # origin row, as in the reviewer probe
  cfg <- list(window = 60L, delta = 1 / 252, horizons = c(1L, 5L),
              origin_step = 1L, origin_dates = p$date[t], calendar = mode)
  scenarios <- list(
    baseline = p,
    future_B_na = within(p, B[t + 2L] <- NA),             # non-target, after origin
    future_B_values = within(p, B[(t + 1L):nrow(p)] <- 9), # non-target values after origin
    future_target_na = within(p, A[t + 5L] <- NA))         # the h = 5 target itself
  for (sc in names(scenarios)) {
    fc <- tryCatch(run_all_models(scenarios[[sc]], LETTERS[1:5], cfg, model_names = c("fBm", "mfBm5", "HAR", "VHAR5")),
                   error = function(e) { message(cal, "/", mode, "/", sc, ": ", conditionMessage(e)); NULL })
    if (is.null(fc) || !nrow(fc)) next
    rows[[length(rows) + 1L]] <- data.frame(calendar = cal, mode = mode, scenario = sc,
      fc[, c("model", "h", "origin_date", "target_date", "forecast", "actual")])
  }
}
setwd(old)
res <- do.call(rbind, rows)
utils::write.csv(res, "docs/review_B/b1_probe.csv", row.names = FALSE)
show <- res[res$h == 5L & res$model %in% c("fBm", "mfBm5"), ]
show$origin_date <- format(show$origin_date); show$target_date <- format(show$target_date)
print(show, row.names = FALSE, digits = 7)
