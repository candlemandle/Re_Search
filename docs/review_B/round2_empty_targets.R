# Second review: origins with no evaluable horizon.
# Usage (project root): Rscript --vanilla docs/review_B/round2_empty_targets.R <label>
# Writes docs/review_B/round2_<label>.csv (outcome of every case) and
# docs/review_B/round2_<label>_dj30_vhar4.rds (forecasts of the DJ30 case, if any).
source("R/code_b_config.R"); source_code_b()
label <- commandArgs(trailingOnly = TRUE)[1]
outcome <- function(case, expr) {
  r <- tryCatch(list(value = force(expr), msg = ""), error = function(e) list(value = NULL, msg = conditionMessage(e)))
  data.frame(case = case, status = if (nzchar(r$msg)) "error" else "ok",
             rows = if (is.data.frame(r$value)) nrow(r$value) else NA_integer_,
             message = r$msg)
}
dj30 <- load_logvol_panel("dj30")
old <- setwd(tempdir()) # run_all_models logs under results/
on.exit(setwd(old))

# synthetic 5-asset mfBm panel on a Monday-Friday calendar (as tests/testthat/test-calendar.R)
set.seed(802)
H <- c(.12, .2, .28, .35, .42, .17, .31); r <- matrix(.25, 7, 7); diag(r) <- 1
x <- simulate_mfbm(104L, H, r, delta = 1 / 252)[[1]]
Lv <- log(.2) + increments_to_levels(x) * .18; colnames(Lv) <- LETTERS[1:7]
d <- seq(as.Date("2010-01-04"), by = "day", length.out = 210L)
p <- data.frame(date = head(d[!format(d, "%u") %in% c("6", "7")], 105L), Lv[, 1:5], check.names = FALSE)
t <- 70L
q <- p; q$A[c(t + 1L, t + 5L)] <- NA
cfg <- function(h, calendar) list(window = 60L, delta = 1 / 252, horizons = h, origin_step = 1L,
                                  origin_dates = p$date[t], calendar = calendar)
res <- list(
  outcome("1 synthetic trading: one origin, h = 1,5, both targets NA",
          run_all_models(q, LETTERS[1:5], cfg(c(1L, 5L), "trading"))),
  outcome("1b synthetic trading: one origin, h = 1, target NA",
          run_all_models(q, LETTERS[1:5], cfg(1L, "trading"))))

# actual DJ30 panel, full configuration, automatic origins, only long horizons
full <- code_b_config("full")
fcs <- list()
for (cal in c("trading", "common_obs")) {
  cf <- code_b_calendar_config(full, cal); cf$horizons <- c(15L, 20L)
  o <- outcome(paste0("2 DJ30 ", cal, ": automatic origins, h = 15,20, VHAR4"),
               fcs[[cal]] <- run_all_models(dj30, full$dj30_assets, cf, cores = code_b_cores(), model_names = "VHAR4"))
  res[[length(res) + 1L]] <- o
}
res[[length(res) + 1L]] <- outcome("optional: rolling_mean_observed(c(1, NA, 3), 5)",
  { v <- rolling_mean_observed(c(1, NA, 3), 5); stopifnot(length(v) == 3L, all(is.na(v))); data.frame(v = v) })
setwd(old)
out <- do.call(rbind, res)
utils::write.csv(out, sprintf("docs/review_B/round2_%s.csv", label), row.names = FALSE)
if (length(fcs) && all(vapply(fcs, is.data.frame, TRUE)))
  saveRDS(fcs, sprintf("docs/review_B/round2_%s_dj30_vhar4.rds", label))
print(out[, c("case", "status", "rows", "message")], right = FALSE, row.names = FALSE)
