#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("smoke", "paper")) else "smoke"

run <- function(label, script, script_args = character()) {
  cat("\n==", label, "==\n")
  status <- system2(file.path(R.home("bin"), "Rscript"), c(script, script_args))
  if (!identical(status, 0L)) stop(label, " failed with exit status ", status)
}

if (!file.exists("R/code_a_config.R") || !file.exists("R/code_b_config.R")) {
  stop("Run this script from the repository root")
}

run("A/B tests", "scripts/run_tests.R")
run("Code C tests", "code_c/scripts/run_tests.R")

if (mode == "smoke") {
  run("Code A estimates", "scripts/03_estimate_parameters.R", "smoke")
  run("Code A Monte Carlo", "scripts/04_run_monte_carlo.R", "smoke")
  run("Code B forecast simulations", "scripts/05_run_forecast_simulations.R", "smoke")
  run("Code B empirical forecasts", "scripts/06_run_empirical_forecasts.R", "smoke")
  run("Code B robustness", "scripts/07_run_robustness.R", "smoke")
  run("Code C end-to-end", "code_c/scripts/run_experiments.R", c("smoke", "all"))
} else {
  run("Code A author-code check", "scripts/00_check_author_code.R")
  run("Code A estimates", "scripts/03_estimate_parameters.R", "full")
  run("Code A Table 1 SE check", "scripts/03b_check_table01_se.R")
  run("Code A missing-day sensitivity", "scripts/03c_missing_day_sensitivity.R")
  run("Code A saved Monte Carlo comparison", "scripts/04b_compare_mc_with_paper.R", "full")
  run("Code B forecast simulations", "scripts/05_run_forecast_simulations.R", "full")
  run("Code B empirical forecasts", "scripts/06_run_empirical_forecasts.R", "full")
  run("Code B robustness", "scripts/07_run_robustness.R", "full")
  cat("\nPaper pipeline complete. Full Code A Monte Carlo and Code C studies are separate\n",
      "because they are substantially more expensive; see README.md.\n", sep = "")
}
