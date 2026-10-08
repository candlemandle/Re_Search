# code B step 5: optimal forecast theory and simulations.
#   figures 4-5 : bivariate forecast weights and msfe with one observation (prop. 4.1)
#   figure 6    : msfe ratio vs dimension d (prop. 4.3)
#   tables 8-11 : monte carlo RMSFE of fBm / bfBm / mfBm3 / mfBm4 (appendix E.3)
#
# run from the project root:
#   Rscript scripts/05_run_forecast_simulations.R smoke   (50 reps, results/smoke/)
#   Rscript scripts/05_run_forecast_simulations.R full    (10000 reps as in the paper)

source("R/code_b_config.R")
code_b_relaunch_single_thread("scripts/05_run_forecast_simulations.R")
source_code_b()
mode <- code_b_mode()
cfg <- code_b_config(mode)
code_b_set_mode(mode)
cores <- code_b_cores()
code_b_log("05 forecast simulations, mode = ", mode, ", cores = ", cores)

# =============================================================================
# figures 4 and 5: one observation (B1_t, B2_t), t = 1, h = 1, H1 = 0.4
# =============================================================================

H2_grid <- seq(0.01, 0.99, by = 0.01)
H1 <- 0.4

bivariate_curves <- function(rho) {
  do.call(rbind, lapply(H2_grid, function(H2) {
    f2 <- one_obs_forecast(mfbm_params(c(H1, H2), rho), t = 1, h = 1)
    f1 <- one_obs_forecast(mfbm_params(H1, matrix(1)), t = 1, h = 1)
    data.frame(
      rho = rho, H2 = H2,
      w11 = f2$w[1], w12 = f2$w[2],
      rel_w11 = f2$w[1] / sum(abs(f2$w)),
      rel_w12 = f2$w[2] / sum(abs(f2$w)),
      w11_vs_univariate = f2$w[1] / w_fbm(1, 1, H1), # t^(2 H1) = 1 for t = 1
      rel_msfe = f2$msfe / f1$msfe,
      admissible = abs(rho) <= rho_max(H1, H2)
    )
  }))
}
curves <- rbind(bivariate_curves(0.5), bivariate_curves(0.9))
code_b_write(curves, "code_b_fig04_05_weights_msfe.csv")

# solid black where the mfBm exists, dotted red where it does not (as in the paper)
plot_curve <- function(x, y, ok, ylab, main) {
  plot(x, y, type = "n", xlab = "H2", ylab = ylab, main = main)
  y_ok <- ifelse(ok, y, NA)
  y_bad <- ifelse(ok, NA, y)
  lines(x, y_ok, lwd = 2)
  lines(x, y_bad, lwd = 2, lty = 3, col = "red")
  abline(v = H1, col = "grey")
}

c05 <- curves[curves$rho == 0.5, ]
pdf(code_b_fig("code_b_fig04_weights.pdf"), width = 10, height = 3.5)
par(mfrow = c(1, 3))
plot_curve(c05$H2, c05$rel_w11, c05$admissible, "relative weight", "w11 / (|w11| + |w12|)")
plot_curve(c05$H2, c05$rel_w12, c05$admissible, "relative weight", "w12 / (|w11| + |w12|)")
plot_curve(c05$H2, c05$w11_vs_univariate, c05$admissible, "ratio", "w11 / w11 (rho = 0)")
dev.off()

pdf(code_b_fig("code_b_fig05_rel_msfe.pdf"), width = 8, height = 3.5)
par(mfrow = c(1, 2))
for (r in c(0.5, 0.9)) {
  cc <- curves[curves$rho == r, ]
  plot_curve(cc$H2, cc$rel_msfe, cc$admissible, "MSFE bfBm / MSFE fBm", paste("rho =", r))
}
dev.off()

# =============================================================================
# figure 6: prop. 4.3, equal correlations rho = 0.8, h = 1, t = 1 and 10.
# the text says H1 = 0.1 and H = 0.4, the caption says H1 = 0.4 and H = 0.1.
# only the caption gives the numbers quoted in the text (89% and 77.4%) so we use
# it. the limits for the other order are saved too
# =============================================================================

H1_fig6 <- 0.4
H_fig6 <- 0.1
fig6 <- do.call(rbind, lapply(c(1, 10), function(t) {
  do.call(rbind, lapply(cfg$fig6_dims, function(d) {
    general <- one_obs_forecast(equicorr_params(d, H1_fig6, H_fig6, 0.8), t, 1)$msfe
    closed <- prop43_msfe(d, t, 1, H1_fig6, H_fig6, 0.8)
    data.frame(t = t, d = d, msfe = general, msfe_closed_form = closed,
               ratio = general / prop43_msfe(1, t, 1, H1_fig6, H_fig6, 0.8))
  }))
}))
# d -> infinity: the closed form at a huge d
limit_ratio <- function(t, H1, H) prop43_msfe(1e6, t, 1, H1, H, 0.8) / prop43_msfe(1, t, 1, H1, H, 0.8)
fig6_limit <- data.frame(
  t = c(1, 10),
  ratio_limit = sapply(c(1, 10), limit_ratio, H1 = H1_fig6, H = H_fig6),
  ratio_limit_text_order = sapply(c(1, 10), limit_ratio, H1 = H_fig6, H = H1_fig6),
  paper_text = c(0.89, 0.774)
)
code_b_write(fig6, "code_b_fig06_dimension.csv")
code_b_write(fig6_limit, "code_b_fig06_limits.csv")
code_b_log("figure 6 max |general - closed form| = ",
           signif(max(abs(fig6$msfe - fig6$msfe_closed_form)), 3))

pdf(code_b_fig("code_b_fig06_dimension.pdf"), width = 8, height = 3.5)
par(mfrow = c(1, 2))
for (t in c(1, 10)) {
  f <- fig6[fig6$t == t, ]
  plot(f$d, f$ratio, type = "b", pch = 20, xlab = "dimension d", ylab = "MSFE ratio",
       main = paste0("t = ", t, ", h = 1"))
  abline(h = fig6_limit$ratio_limit[fig6_limit$t == t], lty = 2)
}
dev.off()

# =============================================================================
# tables 8-11: monte carlo, n = 500, delta = 1/250, h = 1..5
# =============================================================================

n <- cfg$sim_n
hs <- cfg$sim_horizons
delta <- cfg$sim_delta

# one experiment: simulate paths once, then evaluate every forecaster on them
#   models: named list, each element = components the forecaster uses
# every finished experiment is saved to <results>/simulations/checkpoints and a
# rerun loads it only when settings, simulation count and model code agree.
run_experiment <- function(table, setting, H, rho, models, targets, seed, estimate = FALSE) {
  dir.create(code_b_results("simulations", "checkpoints"), recursive = TRUE, showWarnings = FALSE)
  identity <- code_b_fingerprint(list(table = table, setting = setting, H = H,
    rho = rho, models = models, targets = targets, seed = seed, estimate = estimate,
    n = n, horizons = hs, delta = delta, reps = cfg$sim_reps))
  file <- code_b_results("simulations", "checkpoints",
                         paste0(gsub("[^A-Za-z0-9.]+", "_", paste(table, setting)), "_", identity, ".rds"))
  if (file.exists(file)) {
    code_b_log("  ", table, " ", setting, " loaded from checkpoint")
    return(readRDS(file))
  }
  t0 <- Sys.time()
  set.seed(seed)
  par <- mfbm_params(H, rho)
  paths <- simulate_mfbm(n + max(hs), H, rho, delta = delta, nsim = cfg$sim_reps)
  out <- list()
  for (i in targets) for (m in names(models)) {
    use <- models[[m]]
    if (!(i %in% use)) use <- i # fBm of component 2 only uses component 2
    kinds <- if (estimate) c(known = FALSE, unknown = TRUE) else c(mc = FALSE)
    for (kind in names(kinds)) {
      r <- sim_rmsfe(paths, par, use, i, n, hs, delta, estimate = kinds[[kind]], cores = cores)
      out[[length(out) + 1]] <- data.frame(table = table, setting = setting, model = m,
                                           i = i, kind = kind, r)
    }
  }
  code_b_log(sprintf("  %s %s: %.1f s", table, setting,
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res <- do.call(rbind, out)
  saveRDS(res, file)
  res
}

biv <- list(fBm = 1, bfBm = c(1, 2))
sim <- list()

# table 8: H = (0.1, 0.4), rho = 0, 0.4, 0.8
for (r in c(0, 0.4, 0.8)) {
  sim[[length(sim) + 1]] <- run_experiment("Table 8", paste0("rho=", r), c(0.1, 0.4), r,
                                           biv, 1:2, cfg$seed + 10 * r)
}

# table 9: rho = 0.4, H1 = 0.1, H2 = 0.1, 0.2, 0.4
for (H2 in c(0.1, 0.2, 0.4)) {
  sim[[length(sim) + 1]] <- run_experiment("Table 9", paste0("H=(0.1;", H2, ")"), c(0.1, H2), 0.4,
                                           biv, 1:2, cfg$seed + 100 + 10 * H2)
}

# table 10: more components, forecast of B1 only
rho_sigma1 <- matrix(c(1, 0.4, 0.4,
                       0.4, 1, 0,
                       0.4, 0, 1), 3, 3)
rho_sigma2 <- diag(4)
rho_sigma2[1, 2:4] <- rho_sigma2[2:4, 1] <- 0.4
sim[[length(sim) + 1]] <- run_experiment(
  "Table 10", "Sigma1", c(0.1, 0.4, 0.4), rho_sigma1,
  list(fBm = 1, bfBm = 1:2, mfBm3 = 1:3), 1, cfg$seed + 200)
sim[[length(sim) + 1]] <- run_experiment(
  "Table 10", "Sigma2", c(0.1, 0.4, 0.4, 0.4), rho_sigma2,
  list(fBm = 1, bfBm = 1:2, mfBm3 = 1:3, mfBm4 = 1:4), 1, cfg$seed + 300)

# table 11: H = (0.2, 0.4), parameters known vs estimated on every path
for (r in c(0, 0.4)) {
  sim[[length(sim) + 1]] <- run_experiment("Table 11", paste0("rho=", r), c(0.2, 0.4), r,
                                           biv, 1:2, cfg$seed + 400 + 10 * r, estimate = TRUE)
}

sim <- do.call(rbind, sim)
sim$reps <- cfg$sim_reps

# comparison with the paper. monte carlo error of an RMSFE around 0.48 with
# 10000 reps is about 0.48 / sqrt(2 * 10000) = 0.0034
paper <- utils::read.csv("docs/paper_values/code_b_simulation_tables.csv")
cmp <- merge(sim, paper, by = c("table", "setting", "model", "i", "kind", "h"), all.x = TRUE)
cmp$diff_mc <- cmp$rmsfe - cmp$paper_value
cmp$diff_theory <- cmp$theory - cmp$paper_theory
cmp$mc_se <- cmp$rmsfe / sqrt(2 * cmp$reps)
cmp <- cmp[order(cmp$table, cmp$setting, cmp$i, cmp$model, cmp$kind, cmp$h), ]
code_b_write(sim, "code_b_tables08_11_simulations.csv")
code_b_write(cmp, "code_b_tables08_11_vs_paper.csv")

code_b_log("theory vs paper: max |diff| = ", signif(max(abs(cmp$diff_theory), na.rm = TRUE), 3))
code_b_log("monte carlo vs paper: max |diff| / mc se = ",
           signif(max(abs(cmp$diff_mc) / cmp$mc_se, na.rm = TRUE), 3))
code_b_log("05 done")
