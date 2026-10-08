# Cell-level check of the 25 Table 1 standard errors.
# Usage: Rscript scripts/03b_check_table01_se.R
# Output: results/tables/code_a_table01_se_check.csv
#
# Columns
#   se_limit        our SE: closed-form limits (9), (14), (27), N = 5007 increments
#   se_limit_R5000 / se_limit_R1e6   same with 5000 / 1e6 series terms (convergence check)
#   se_author_n2000_N5007 / _N5008   authors' BYZinference.R functions (finite sums with
#                   (1 - r/n) weights, truncation n = 2000) with N = 5007 increments or
#                   N = 5008 observed levels
#   paper_se        Table 1 value; match_* = equal after rounding to 4 decimals
# The author-function / N = 5008 variant is a post hoc explanation of the two cells
# that differ in the 4th decimal; it is not stated in the paper.

source("R/code_a_config.R")
source_code_a()
cfg <- code_a_config("full")
code_a_set_mode("full")
code_a_log("03b_check_table01_se: start")

auth <- new.env()
sys.source(file.path("author_code", "BYZ", "Setup code.R"), envir = auth)
inf <- new.env(parent = auth)
# function definitions only (the rest of BYZinference.R computes an example)
exprs <- parse(file.path("author_code", "BYZ", "BYZinference.R"))
for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                     is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, inf)
quiet <- function(expr) { utils::capture.output(v <- expr); v } # author functions print()

panel <- load_logvol_panel("dj30")
tk <- cfg$table1_tickers
X <- panel_increments(panel, tk)$X
N <- nrow(X)
est <- estimate_mfbm(X, cfg$delta_daily)
ref <- utils::read.csv(file.path("docs", "paper_values", "table01_se.csv"), stringsAsFactors = FALSE)
ref$asset2[is.na(ref$asset2)] <- ""

se_ours <- function(i, j, stat, R) {
  if (stat == "H") return(sqrt(avar_hurst(est$H[i], R) / N))
  f <- if (stat == "rho") avar_rho else avar_eta
  sqrt(f(est$H[i], est$H[j], est$rho[i, j], R) / N)
}
se_author <- function(i, j, stat, NN) {
  if (stat == "H") return(sqrt(inf$varHest(2000, est$H[i]) / NN))
  f <- if (stat == "rho") inf$varrho2 else inf$vareta
  sqrt(quiet(f(2000, est$rho[i, j], est$H[i], est$H[j])) / NN)
}
out <- do.call(rbind, lapply(seq_len(nrow(ref)), function(k) {
  i <- ref$asset1[k]; j <- ref$asset2[k]; s <- ref$statistic[k]
  data.frame(asset1 = i, asset2 = j, statistic = s, N = N,
             se_limit = se_ours(i, j, s, cfg$series_terms),
             se_limit_R5000 = se_ours(i, j, s, 5000L),
             se_limit_R1e6 = se_ours(i, j, s, 1000000L),
             se_author_n2000_N5007 = se_author(i, j, s, N),
             se_author_n2000_N5008 = se_author(i, j, s, N + 1),
             paper_se = ref$paper_se[k])
}))
for (v in c("se_limit", "se_author_n2000_N5007", "se_author_n2000_N5008")) {
  out[[paste0("match_", v)]] <- abs(round(out[[v]], 4) - out$paper_se) < 1e-9
}
out$convergence_max_rel_diff <- pmax(abs(out$se_limit_R5000 / out$se_limit - 1),
                                     abs(out$se_limit_R1e6 / out$se_limit - 1))
code_a_write(out, "code_a_table01_se_check.csv")
code_a_log(sprintf("  rounded matches: limit %d/25, author n=2000 N=%d %d/25, author n=2000 N=%d %d/25; max rel. truncation effect %.1e",
                   sum(out$match_se_limit), N, sum(out$match_se_author_n2000_N5007), N + 1,
                   sum(out$match_se_author_n2000_N5008), max(out$convergence_max_rel_diff)))
code_a_log_run("03b_check_table01_se", "full", cfg)
