# code B configuration: forecasting (section 4, section 5.2-5.3, appendix E.3 and F).
# smoke and full runs share the same code, only the numbers below differ.

code_b_config <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)

  cfg <- list(
    mode = mode,
    seed = 20251004L,

    # --- empirical forecasts (tables 2, 3, 15-19) ---------------------------
    delta = 1 / 252, # same daily step as code A
    # "two-year" rolling window. 500 days because both samples have exactly 500
    # common days before the first forecast date of the paper (2007-01-03 and 2016-03-29).
    # on the trading calendar the same 500 rows are 500 trading days (2 x 252 = 504)
    window = 500L,
    # time convention of origins, windows and h (review B1, see code_b_calendars()):
    #   "trading"    primary: h = trading days of the panel calendar, targets fixed
    #                at the origin, a missing target is dropped for all models
    #   "common_obs" historical replication mode: rows = days on which all assets of
    #                the set are observed (a future missing day moves the target)
    calendar = "trading",
    calendars = c("trading", "common_obs"), # conventions computed by scripts 06 and 07
    min_window_coverage = 0.9, # trading calendar: observed share of the window
    horizons = c(1L, 2L, 3L, 4L, 5L, 10L, 15L, 20L),
    origin_step = 1L, # forecast every day. smoke mode skips days
    dj30_assets = c("AAPL", "ALD", "AMGN", "AXP", "BA"),
    mag7_assets = c("AAPL", "AMZN", "FB", "GOOG", "MSFT"),
    # forecast periods are defined by the date being forecast (target date)
    dj30_periods = list(
      full = c("2007-01-03", "2025-01-14"),
      period1 = c("2007-01-03", "2017-04-11"),
      period2 = c("2017-04-12", "2021-07-30")
    ),
    mag7_periods = list(full = c("2016-03-29", "2025-01-14")),

    # --- model confidence set (appendix F.4) --------------------------------
    mcs_boot = 5000L,
    mcs_block = 20L,

    # --- simulations (appendix E.3, tables 8-11) ----------------------------
    sim_n = 500L,
    sim_delta = 1 / 250,
    sim_horizons = 1:5,
    sim_reps = 10000L,

    # --- figures 4-6 ---------------------------------------------------------
    fig6_dims = 1:100
  )

  if (mode == "smoke") {
    cfg$origin_step <- 500L
    cfg$mcs_boot <- 500L
    cfg$sim_reps <- 50L
    cfg$fig6_dims <- c(1:10, 20, 50, 100)
  }
  cfg
}

# "smoke" or "full" from the command line (default smoke)
code_b_mode <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) return("smoke")
  match.arg(args[1], c("smoke", "full"))
}

# number of forked workers (1 on windows). MFBM_CORES overrides it like in code A
code_b_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  env <- Sys.getenv("MFBM_CORES", "")
  if (nzchar(env)) return(as.integer(env))
  detected <- parallel::detectCores()
  if (is.na(detected)) return(1L)
  max(1L, detected - 1L)
}

# forked workers plus multithreaded openblas overload the cpu (each window became
# ~7x slower). blas reads its thread count at start so we restart the script once
# with single-threaded blas. same trick as scripts/04_run_monte_carlo.R
code_b_relaunch_single_thread <- function(script) {
  if (Sys.getenv("OPENBLAS_NUM_THREADS") != "") return(invisible(FALSE))
  status <- system2(file.path(R.home("bin"), "Rscript"), c(script, commandArgs(trailingOnly = TRUE)),
                    env = c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "VECLIB_MAXIMUM_THREADS=1"))
  quit(save = "no", status = status)
}

# load code A (estimators, covariances, simulator, data) and the code B files
source_code_b <- function(root = ".") {
  options(code_b.source_root = normalizePath(root))
  source(file.path(root, "R", "code_a_config.R"))
  source_code_a(root)
  for (f in c("R/forecast_mfbm.R", "R/forecast_har.R", "R/rolling_window.R", "R/metrics.R")) {
    source(file.path(root, f))
  }
  invisible(TRUE)
}

# smoke results go to results/smoke/ so they never overwrite full results
code_b_set_mode <- function(mode) {
  options(code_b.results = if (mode == "smoke") file.path("results", "smoke") else "results")
  invisible(getOption("code_b.results"))
}

code_b_results <- function(...) {
  file.path(getOption("code_b.results", "results"), ...)
}

code_b_log <- function(...) {
  dir.create(file.path("results", "logs"), recursive = TRUE, showWarnings = FALSE)
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path("results", "logs", "code_b.log"), append = TRUE)
}

code_b_write <- function(df, name) {
  dir.create(code_b_results("tables"), recursive = TRUE, showWarnings = FALSE)
  path <- code_b_results("tables", name)
  utils::write.csv(df, path, row.names = FALSE)
  code_b_log("wrote ", path)
  invisible(path)
}

code_b_fig <- function(name) {
  dir.create(code_b_results("figures"), recursive = TRUE, showWarnings = FALSE)
  code_b_results("figures", name)
}
