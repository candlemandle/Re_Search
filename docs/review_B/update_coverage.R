# Updates the Code B rows of docs/paper_coverage_code_b.csv after the review-B full run.
p <- "docs/paper_coverage_code_b.csv"; x <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
set <- function(item, status, match, notes) { i <- x$item == item; stopifnot(sum(i) == 1)
  x$status[i] <<- status; x$match_with_paper[i] <<- match; x$notes[i] <<- notes }
cal <- " Primary = trading calendar; *_common_obs files = historical replication mode (review B1). Full run 2026-10-08."
set("Table 2", "reproduced_with_differences",
 "common_obs: h = 1 equal to 4 decimals in full period and period 1 (10/10 cells; trading 3/10); h > 1 within 0.7% (both conventions); mfBm5 gain over fBm within 0.06 pp in full period and period 1. Period 2 levels -1.4% to +9.1% (median +6.6%) (0.4099 vs 0.3842 at fBm h = 1)",
 paste0("Period 2 gap unresolved: paper period 2 is inconsistent with its own full/period-1 values on our forecasts; April 11 vs 12 changes it by 0.08 pp; no single alternative boundary reproduces all 40 cells (exploratory scan). Subperiods chosen after Figure 7: exploratory.", cal))
set("Table 3", "reproduced_with_differences",
 "h <= 5 within 2%; VHAR4/VHAR5 at h >= 10 up to -12% (VHAR4 h = 20: 1.6981 vs 1.9388)",
 paste0("Long-horizon VHAR cells depend on the positivity filter: without it (exploratory, not adopted) VHAR4 h = 20 is 1.929. Not an exact replication.", cal))
set("Table 15", "reproduced_with_differences", "median |rel diff| 0.8%, max 4.3%; long-h cells filter-sensitive",
 paste0("VHARF2 = VHAR2 exactly as in the paper (unit test).", cal))
set("Table 16", "reproduced_with_differences", "median |rel diff| 0.6%, max 2.2%",
 paste0("Lognormal correction exp(fit + s2/2).", cal))
set("Table 17", "reproduced_with_differences",
 "mfBm class +2.4% to +5.3% (median +4.1%), VHARF/VLHAR +0.4% to +4.4%; linear VHAR down to -49.7% (VHAR5 h = 20: 1.2826 vs 2.5475). The 'about 4% higher' statement applies only to the first two groups",
 paste0("mfBm5 beats fBm for h >= 4 (paper h >= 3); VHAR long-h levels depend on the positivity filter (no-filter variant 2.209, exploratory).", cal))
set("Table 18", "reproduced_with_differences",
 "mfBm within 0.44%, VLHAR within 2.6%; VHAR/VHARF at h >= 10 differ by up to -52% (VHAR4 h = 15: 0.1383 vs 0.2896)",
 paste0("Linear-HAR QLIKE is dominated by a few near-zero positive forecasts that the <= 0 filter does not catch; level not robust.", cal))
set("Table 19", "reproduced_with_differences",
 "mfBm +2.1% to +5.3%; VHAR3 h = 10 +128% (0.3390 vs 0.1487)",
 paste0("One forecast (origin 2020-04-08, 0.000485 vs actual 0.270) gives 74% of the VHAR3 h = 10 loss; the level is not robustly estimable.", cal))
write.csv(x, p, row.names = FALSE)
cat("updated", p, "\n")
