# Code A configuration: data, estimation and Monte Carlo settings.
# smoke and full runs share all logic; only the numbers below differ.
# Code C can merge this list into R/config.R without changing call sites.

code_a_config <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)

  base <- list(
    mode = mode,
    seed = 20251003L,

    # --- data -------------------------------------------------------------
    risklab_url = "https://dachxiu.chicagobooth.edu/data.php?ticker=",
    # Paper tickers (Table 13) use the historical Risk Lab names of each PERMNO.
    dj30 = c(
      AAPL = 14593, ALD = 10145, AMGN = 14008, AXP = 59176, BA = 19561,
      BEL = 65875, CAT = 18542, CHV = 14541, CRM = 90215, CSCO = 76076,
      DIS = 26403, GS = 86868, HD = 66181, IBM = 12490, INTC = 59328,
      JNJ = 22111, JPM = 47896, KO = 11308, MCD = 43449, MMM = 22592,
      MRK = 22752, MSFT = 10107, NIKE = 57665, PG = 18163, SPC = 59459,
      UNH = 92655, V = 92611, WAG = 19502, WMT = 55976, XOM = 11850
    ),
    mag7 = c(AAPL = 14593, AMZN = 84788, FB = 13407, GOOG = 14542, MSFT = 10107),
    # Column of the Risk Lab file used as daily volatility.
    # "qmle_trade" = Risk Lab headline estimator (annualized volatility).
    rv_measure = "qmle_trade",
    dj30_start = as.Date("2005-01-05"),
    dj30_end = as.Date("2025-01-14"),
    mag7_start = as.Date("2014-03-27"),
    mag7_end = as.Date("2025-01-14"),
    delta_daily = 1 / 252,
    table1_tickers = c("AAPL", "ALD", "AMGN", "AXP", "BA"),
    rolling_window = 504L, # two trading years (Figures 7-8)
    rolling_step = 5L,

    # --- asymptotic series truncation --------------------------------------
    series_terms = 20000L,

    # --- Monte Carlo --------------------------------------------------------
    mc = list(
      # Tables 4-5 (Appendix E.1)
      e1_n = c(500L, 1000L),
      e1_delta = c(1 / 52, 1 / 250),
      e1_rho = c(0, 0.4),
      e1_reps = 1000L,
      # Table 6 (BYZ vs Amblard-Coeurjolly)
      t6_n = c(500L, 1000L),
      t6_grid = expand.grid(rho = c(0.2, 0.4), eta = c(0, 0.5)),
      t6_reps = 1000L,
      ac_dilations = 1:5,
      # Table 7 (size/power of the time-reversibility test)
      t7_n = c(500L, 1000L),
      t7_eta = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65),
      t7_reps = 5000L,
      # Table 12 (measurement error)
      t12_n = c(500L, 1000L),
      t12_intraday = 72L,
      t12_reps = 1000L,
      # Appendix G (Table 20, Figures 9-10)
      g_n = 1000L,
      g_lags = 10L,
      g_reps = 1000L,
      g_boot = 1000L,
      g_reps_cheap = 10000L # extra run for rho_lag / rho^(2) variances only
    )
  )

  if (mode == "smoke") {
    base$mc$e1_reps <- 50L
    base$mc$t6_reps <- 50L
    base$mc$t7_reps <- 200L
    base$mc$t12_reps <- 50L
    base$mc$g_reps <- 50L
    base$mc$g_boot <- 50L
    base$mc$g_reps_cheap <- 200L
    base$series_terms <- 5000L
    base$rolling_step <- 21L
  }
  base
}

# Parse "smoke"/"full" from the command line (default: smoke).
code_a_mode <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) return("smoke")
  match.arg(args[1], c("smoke", "full"))
}

# Number of worker processes for Monte Carlo (forked; 1 on Windows).
code_a_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  env <- Sys.getenv("MFBM_CORES", "")
  if (nzchar(env)) return(as.integer(env))
  max(1L, parallel::detectCores() - 1L)
}

# Source the Code A library files (run scripts from the project root).
source_code_a <- function(root = ".") {
  for (f in c("R/code_a_config.R", "R/covariance.R", "R/estimators.R",
              "R/time_reversibility.R", "R/data_processing.R", "R/monte_carlo.R")) {
    source(file.path(root, f))
  }
  invisible(TRUE)
}

# Append a timestamped line to results/logs/code_a.log and echo it.
code_a_log <- function(...) {
  dir.create(file.path("results", "logs"), recursive = TRUE, showWarnings = FALSE)
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path("results", "logs", "code_a.log"), append = TRUE)
}

# Output root: results/ for full runs, results/smoke/ for smoke runs, so a smoke
# run never overwrites full results. Set once per script via code_a_set_mode().
code_a_set_mode <- function(mode) {
  options(code_a.results = if (mode == "smoke") file.path("results", "smoke") else "results")
  invisible(getOption("code_a.results"))
}

code_a_results <- function(...) {
  file.path(getOption("code_a.results", "results"), ...)
}

# Path of a figure under <output root>/figures/ (directory created on demand).
code_a_fig <- function(name) {
  dir.create(code_a_results("figures"), recursive = TRUE, showWarnings = FALSE)
  code_a_results("figures", name)
}

# Write a results table as CSV under <output root>/tables/.
code_a_write <- function(df, name) {
  dir.create(code_a_results("tables"), recursive = TRUE, showWarnings = FALSE)
  path <- code_a_results("tables", name)
  utils::write.csv(df, path, row.names = FALSE)
  code_a_log("wrote ", path)
  invisible(path)
}
