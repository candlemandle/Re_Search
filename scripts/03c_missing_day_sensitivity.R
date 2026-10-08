# Sensitivity of the empirical estimates to the missing-day convention.
# Usage: Rscript scripts/03c_missing_day_sensitivity.R
# Outputs: results/tables/code_a_missing_day_sensitivity.csv (agreement with Tables 1, 13, 14
#          under alternative rules) and code_a_missing_day_gaps.csv (bridged gaps per series).
#
# All rules treat consecutive *available* observations as one sampling step Delta,
# i.e. a missing day is collapsed rather than modelled as a 2*Delta increment. This is a
# replication convention, not an exactly equally spaced model. Which rule the paper used
# is not stated; the rules below are hypotheses compared by their agreement with the
# published tables, not established facts.
#   own_or_pairwise   H on each series' own dates; rho/eta on each pair's common dates (default)
#   common_sample     every estimate on dates common to the whole set of series compared
#   na_grid           increments on the full panel calendar; increments spanning a missing
#                     day are dropped instead of bridged
#   hybrid_table1     default, except the five Table 1 stocks use the five-stock common sample
#                     for H and for every pair involving one of them

source("R/code_a_config.R")
source_code_a()
cfg <- code_a_config("full")
code_a_set_mode("full")
code_a_log("03c_missing_day_sensitivity: start")

panel <- load_logvol_panel("dj30")
tk <- names(cfg$dj30)
f5 <- cfg$table1_tickers
ref13 <- utils::read.csv(file.path("docs", "paper_values", "table13_H_rho.csv"), stringsAsFactors = FALSE)
ref14 <- utils::read.csv(file.path("docs", "paper_values", "table14_eta.csv"), stringsAsFactors = FALSE)
ref13$asset2[is.na(ref13$asset2)] <- ""

grid_inc <- sapply(tk, function(t) diff(panel[[t]])) # NA when either day is missing
na_rho <- function(a, b) { ok <- !is.na(a) & !is.na(b); rho_mm(a[ok], b[ok]) }
na_eta <- function(a, b) {
  # lagged products need both days of both series; drop pairs with any NA
  n <- length(a)
  ok <- !is.na(a) & !is.na(b)
  ok2 <- ok[-1] & ok[-n]
  num <- sum((a[-1] * b[-n] - b[-1] * a[-n])[ok2])
  l2 <- function(x) x[-n] + x[-1]
  okl <- ok2
  den <- sqrt(sum(l2(a)[okl]^2) * sum(l2(b)[okl]^2)) - 2 * sqrt(sum(a[ok]^2) * sum(b[ok]^2))
  num / den
}
na_hurst <- function(a) {
  n <- length(a); ok <- !is.na(a); ok2 <- ok[-1] & ok[-n]
  log(sum(((a[-n] + a[-1])[ok2])^2) / sum(a[ok]^2)) / (2 * log(2))
}

est_rule <- function(rule) {
  common <- if (rule == "common_sample") panel_increments(panel, tk)$X else NULL
  common5 <- panel_increments(panel, f5)$X
  H <- vapply(tk, function(t) switch(rule,
    own_or_pairwise = hurst_mm(panel_increments(panel, t)$X[, 1]),
    common_sample = hurst_mm(common[, t]),
    na_grid = na_hurst(grid_inc[, t]),
    hybrid_table1 = if (t %in% f5) hurst_mm(common5[, t]) else hurst_mm(panel_increments(panel, t)$X[, 1])), 0)
  pairs <- utils::combn(tk, 2)
  pr <- apply(pairs, 2, function(p) {
    X <- switch(rule,
      own_or_pairwise = panel_increments(panel, p)$X,
      common_sample = common[, p],
      na_grid = NULL,
      hybrid_table1 = if (any(p %in% f5)) panel_increments(panel, unique(c(f5, p)))$X[, p]
                      else panel_increments(panel, p)$X)
    if (is.null(X)) c(na_rho(grid_inc[, p[1]], grid_inc[, p[2]]), na_eta(grid_inc[, p[1]], grid_inc[, p[2]]))
    else c(rho_mm(X[, 1], X[, 2]), eta_mm(X[, 1], X[, 2]))
  })
  rbind(data.frame(rule = rule, asset1 = tk, asset2 = "", statistic = "H", estimate = H),
        data.frame(rule = rule, asset1 = pairs[1, ], asset2 = pairs[2, ], statistic = "rho", estimate = pr[1, ]),
        data.frame(rule = rule, asset1 = pairs[1, ], asset2 = pairs[2, ], statistic = "eta", estimate = pr[2, ]))
}

rules <- c("own_or_pairwise", "common_sample", "na_grid", "hybrid_table1")
all_est <- do.call(rbind, lapply(rules, est_rule))
sym <- function(a, b, s) ifelse(s == "rho", paste(pmin(a, b), pmax(a, b), s), paste(a, b, s))
ref <- rbind(ref13, ref14)
all_est$paper_value <- ref$paper_value[match(sym(all_est$asset1, all_est$asset2, all_est$statistic),
                                             sym(ref$asset1, ref$asset2, ref$statistic))]
all_est$abs_diff <- abs(all_est$estimate - all_est$paper_value)
all_est$exact4 <- abs(round(all_est$estimate, 4) - all_est$paper_value) < 1e-9

summ <- do.call(rbind, lapply(split(all_est, paste(all_est$rule, all_est$statistic)), function(d) {
  data.frame(rule = d$rule[1], table = ifelse(d$statistic[1] == "eta", "Table 14", "Table 13"),
             statistic = d$statistic[1], n_cells = nrow(d), share_equal_4dp = mean(d$exact4),
             median_abs_diff = stats::median(d$abs_diff), max_abs_diff = max(d$abs_diff))
}))
# Table 1: H of the five stocks, own dates vs five-stock common sample
X5 <- panel_increments(panel, f5)$X
t1 <- data.frame(rule = c("own_or_pairwise", "common_sample"), table = "Table 1", statistic = "H",
                 n_cells = 5,
                 share_equal_4dp = c(mean(abs(round(vapply(f5, function(t) hurst_mm(panel_increments(panel, t)$X[, 1]), 0), 4) -
                                              ref13$paper_value[match(f5, ref13$asset1)]) < 1e-9),
                                     mean(abs(round(apply(X5, 2, hurst_mm), 4) - ref13$paper_value[match(f5, ref13$asset1)]) < 1e-9)),
                 median_abs_diff = NA, max_abs_diff = c(
                   max(abs(vapply(f5, function(t) hurst_mm(panel_increments(panel, t)$X[, 1]), 0) - ref13$paper_value[match(f5, ref13$asset1)])),
                   max(abs(apply(X5, 2, hurst_mm) - ref13$paper_value[match(f5, ref13$asset1)]))))
summ <- rbind(t1, summ)
summ$status <- "hypothesis comparison (paper does not state its rule)"
code_a_write(summ, "code_a_missing_day_sensitivity.csv")

gaps <- do.call(rbind, lapply(tk, function(t) {
  ok <- which(!is.na(panel[[t]]))
  data.frame(ticker = t, first = panel$date[min(ok)], observations = length(ok),
             missing_panel_dates_inside_span = (max(ok) - min(ok) + 1) - length(ok),
             bridged_increments = sum(diff(ok) > 1))
}))
code_a_write(gaps, "code_a_missing_day_gaps.csv")
print(summ[, c("rule", "table", "statistic", "share_equal_4dp", "max_abs_diff")], digits = 3, row.names = FALSE)
code_a_log_run("03c_missing_day_sensitivity", "full", cfg)
