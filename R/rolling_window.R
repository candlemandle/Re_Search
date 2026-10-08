# rolling-window forecasting pipeline: one function runs every model on the
# same forecast origins, so all models have the same test dates, target and h.

# all models of the paper for one asset set (assets[1] is the series we forecast).
# class is the group used for the MCS p-values (tables 2, 3, 15, 16, 17).
# calendar (review B1): "common_obs" rows are the common observation days of all
# assets (historical replication convention); "trading" rows are trading days of
# the panel and NA marks a day on which an asset is not observed
code_b_models <- function(assets, window, delta, calendar = "common_obs") {
  calendar <- match.arg(calendar, code_b_calendars())
  na_avg <- calendar == "trading"
  mfbm_fun <- if (na_avg) {
    function(L, h) forecast_mfbm_calendar(utils::tail(L, window), h, delta)
  } else {
    function(L, h) forecast_mfbm_window(utils::tail(L, window), h, delta)
  }
  model <- function(name, class, k, fun) list(name = name, class = class, assets = assets[1:k],
                                              window = window, fun = fun)
  suffix <- c("", "2", "3", "4", "5")
  out <- list()
  for (k in seq_len(min(length(assets), 5))) {
    # local() so every closure keeps its own k
    out <- c(out, local({
      k <- k
      list(
        model(c("fBm", "bfBm", "mfBm3", "mfBm4", "mfBm5")[k], "mfbm", k, mfbm_fun),
        model(paste0(if (k == 1) "HAR" else "VHAR", suffix[k]), "vhar", k,
              function(L, h) forecast_har_window(L, h, window, "vhar", na_avg = na_avg)),
        model(paste0(if (k == 1) "HAR" else "VHARF", suffix[k]), "vharf", k,
              function(L, h) forecast_har_window(L, h, window, "vharf", na_avg = na_avg)),
        model(paste0(if (k == 1) "LHAR" else "VLHAR", suffix[k]), "vlhar", k,
              function(L, h) forecast_har_window(L, h, window, "vhar", log = TRUE, na_avg = na_avg))
      )
    }))
  }
  # HAR sits in two classes (tables 3 and 15). it is the same model so we keep one copy
  names(out) <- vapply(out, `[[`, "", "name")
  out[!duplicated(names(out))]
}

# the two time conventions for origins, training windows and h (review B1)
#   common_obs : historical replication mode. rows with a missing value in any asset
#                of the set are deleted first, so h counts common observation days and
#                a missing future day of another asset moves the target
#   trading    : h counts trading days of the panel calendar. the target of an origin
#                is fixed in advance; a missing target removes the forecast for all
#                models, it is never moved to another day
code_b_calendars <- function() c("common_obs", "trading")

# convention of a configuration. NULL (older configurations) is the historical one
code_b_calendar <- function(cfg) {
  if (is.null(cfg$calendar)) return("common_obs")
  match.arg(cfg$calendar, code_b_calendars())
}

# trading-calendar origins that can be used. decided only from rows 1..t:
#   1. a full window of `window` trading days ends at t
#   2. every asset of the set is observed at t (the conditioning vector of the
#      largest model and the daily HAR terms); applied to all models together so
#      that every model has the same forecast keys
#   3. at least min_coverage of the window's days are observed for all assets
# returns a data.frame with the origin rows and the reason when it is skipped
code_b_eligible_origins <- function(L, origins, window, min_coverage = 0.9) {
  complete <- stats::complete.cases(L)
  cov <- vapply(origins, function(t) if (t < window) NA_real_ else mean(complete[(t - window + 1):t]), 0)
  reason <- ifelse(origins < window, "short_window",
            ifelse(!complete[origins], "asset_missing_at_origin",
            ifelse(cov < min_coverage, "low_window_coverage", "")))
  data.frame(origin = origins, coverage = cov, eligible = reason == "", reason = reason)
}

# model classes of the tables, univariate model first
code_b_classes <- function() {
  list(
    mfbm = c("fBm", "bfBm", "mfBm3", "mfBm4", "mfBm5"),
    vhar = c("HAR", "VHAR2", "VHAR3", "VHAR4", "VHAR5"),
    vharf = c("HAR", "VHARF2", "VHARF3", "VHARF4", "VHARF5"),
    vlhar = c("LHAR", "VLHAR2", "VLHAR3", "VLHAR4", "VLHAR5")
  )
}

# forecast origins: row numbers t of the panel. the first one has `window` rows of
# history and we stop when even h = 1 has no realised value
forecast_origins <- function(n_rows, window, step = 1L) {
  seq(window, n_rows - 1L, by = step)
}

# run one model over all origins.
#   L     : log vol matrix (rows = days of the chosen calendar, column 1 = target)
#   dates : dates of the rows of L
# the model only ever sees L[1:t, ] so it can not use the future.
# the target of origin t and horizon h is row t + h; it is dropped (for every model,
# because the actual is shared) when it is outside the panel or not observed, never
# moved. an origin without any evaluable horizon is skipped before the model is
# fitted, so it is omitted for every model alike. returns a typed zero-row table
# (empty_forecast_table()) when no origin has an evaluable horizon
rolling_forecasts <- function(L, dates, model, origins, horizons, cores = 1L) {
  window <- if (is.null(model$window)) NA_integer_ else model$window
  one_origin <- function(t) {
    ok <- t + horizons <= nrow(L)
    ok[ok] <- is.finite(L[t + horizons[ok], 1])
    if (!any(ok)) return(NULL) # nothing to evaluate at this origin
    past <- L[1:t, model$assets, drop = FALSE]
    fc <- model$fun(past, horizons)
    data.frame(
      model = model$name,
      class = model$class,
      h = horizons[ok],
      origin_date = dates[t],
      target_date = dates[t + horizons[ok]],
      forecast = fc[ok],
      actual = exp(L[t + horizons[ok], 1]),
      train_start = if (is.na(window)) dates[NA_integer_] else dates[max(1L, t - window + 1L)],
      n_train_days = if (is.na(window)) NA_integer_ else
        sum(stats::complete.cases(past[max(1L, t - window + 1L):t, , drop = FALSE]))
    )
  }
  res <- parallel::mclapply(origins, one_origin, mc.cores = cores)
  failed <- vapply(res, inherits, TRUE, "try-error")
  if (any(failed)) stop(model$name, ": ", res[[which(failed)[1]]])
  res <- res[!vapply(res, is.null, TRUE)]
  if (!length(res)) return(empty_forecast_table(dates))
  do.call(rbind, res)
}

# zero-row forecast table with the same columns and types as rolling_forecasts()
empty_forecast_table <- function(dates) {
  data.frame(model = character(), class = character(), h = integer(),
             origin_date = dates[0], target_date = dates[0], forecast = numeric(),
             actual = numeric(), train_start = dates[0], n_train_days = integer())
}

# run every model on the same origins and stack the results.
# cfg$calendar picks the time convention (see code_b_calendars()):
#   common_obs : only the days where all assets are observed are used, for every model
#   trading    : all days of the panel; origins from code_b_eligible_origins()
# model_names picks a subset of models (NULL = all).
# checkpoint: if a name is given every finished model is saved to
# <results>/forecasts/checkpoints/<name>_<fingerprint>_<model>.rds and a rerun loads
# it instead of computing it again. the fingerprint covers data, configuration
# (calendar included) and the source of every model function
run_all_models <- function(panel, assets, cfg, cores = 1L, model_names = NULL,
                           checkpoint = NULL) {
  calendar <- code_b_calendar(cfg)
  if (anyNA(panel$date) || anyDuplicated(panel$date) || is.unsorted(panel$date, strictly = TRUE))
    stop("panel dates must be unique and sorted")
  if (!is.null(checkpoint)) checkpoint <- paste0(checkpoint, "_", forecast_fingerprint(panel, assets, cfg))
  sub <- if (calendar == "common_obs") {
    panel[stats::complete.cases(panel[, assets]), c("date", assets)]
  } else {
    panel[, c("date", assets)]
  }
  L <- as.matrix(sub[, assets])
  origins <- forecast_origins(nrow(L), cfg$window, cfg$origin_step)
  if (!is.null(cfg$origin_dates)) {
    origins <- match(as.Date(cfg$origin_dates), sub$date)
    if (!length(origins) || anyNA(origins) || anyDuplicated(origins) ||
        any(origins < cfg$window | origins >= nrow(L))) stop("invalid explicit origin dates")
  }
  if (calendar == "trading") {
    el <- code_b_eligible_origins(L, origins, cfg$window, code_b_min_coverage(cfg))
    if (!is.null(cfg$origin_dates) && !all(el$eligible)) stop("explicit origin dates are not eligible")
    if (any(!el$eligible)) {
      code_b_log("  trading calendar: ", sum(!el$eligible), " of ", nrow(el), " origins skipped (",
                 paste(names(table(el$reason[!el$eligible])), table(el$reason[!el$eligible]),
                       sep = " = ", collapse = ", "), ")")
    }
    origins <- origins[el$eligible]
  }
  if (!length(origins)) stop("no eligible forecast origins")
  models <- code_b_models(assets, cfg$window, cfg$delta, calendar)
  if (!is.null(model_names)) models <- models[model_names]

  out <- list()
  for (m in models) {
    file <- NULL
    if (!is.null(checkpoint)) {
      dir.create(code_b_results("forecasts", "checkpoints"), recursive = TRUE, showWarnings = FALSE)
      file <- code_b_results("forecasts", "checkpoints", paste0(checkpoint, "_", m$name, ".rds"))
      if (file.exists(file)) {
        out[[m$name]] <- readRDS(file)
        code_b_log("  ", m$name, " loaded from checkpoint")
        next
      }
    }
    t0 <- Sys.time()
    out[[m$name]] <- rolling_forecasts(L, sub$date, m, origins, cfg$horizons, cores)
    if (!is.null(file)) saveRDS(out[[m$name]], file)
    code_b_log(sprintf("  %-7s %d origins, %.1f s", m$name, length(origins),
                       as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  fc <- do.call(rbind, out)
  rownames(fc) <- NULL
  if (!nrow(fc)) {
    # a run in which no origin has an observed target inside the panel for any h
    stop(structure(class = c("code_b_no_evaluable_forecasts", "error", "condition"), list(
      message = paste0("no evaluable forecasts: none of the ", length(origins), " origins has an observed ",
                       "target within the panel for horizons ", paste(cfg$horizons, collapse = ", "),
                       " (calendar = ", calendar, ")"), call = NULL)))
  }
  fc$calendar <- calendar
  fc
}

# share of the trading window that must be observed (trading calendar only)
code_b_min_coverage <- function(cfg) if (is.null(cfg$min_window_coverage)) 0.9 else cfg$min_window_coverage

# checks asked for in the task. returns a one-row data.frame and stops on a hard error
check_forecasts <- function(fc) {
  if (!nrow(fc)) stop("forecast checks failed: empty forecasts")
  key <- paste(fc$h, fc$origin_date, fc$target_date)
  same_dates <- all(tapply(key, fc$model, function(k) setequal(k, key[fc$model == fc$model[1]])))
  res <- data.frame(
    n_models = length(unique(fc$model)),
    n_rows = nrow(fc),
    same_test_dates = same_dates,
    no_lookahead = all(fc$origin_date < fc$target_date),
    no_na_inf = all(is.finite(fc$forecast) & is.finite(fc$actual)),
    positive_forecasts = all(is.finite(fc$forecast) & fc$forecast > 0),
    unique_keys = !anyDuplicated(paste(fc$model, fc$h, fc$origin_date, fc$target_date)),
    positive_actuals = all(is.finite(fc$actual) & fc$actual > 0),
    same_actuals = all(vapply(split(fc$actual, key), function(x) length(unique(x)) == 1L, TRUE))
  )
  if (!isTRUE(res$same_test_dates) || !isTRUE(res$no_lookahead) ||
      !isTRUE(res$no_na_inf) || !isTRUE(res$positive_forecasts) ||
      !isTRUE(res$unique_keys) || !isTRUE(res$positive_actuals) ||
      !isTRUE(res$same_actuals)) stop("forecast checks failed")
  res
}

# Checkpoints are specific to data, configuration and the loaded model functions.
forecast_fingerprint <- function(panel, assets, cfg) {
  code_b_fingerprint(list(panel = panel[, c("date", assets)], assets = assets, cfg = cfg,
                          calendar = code_b_calendar(cfg), min_coverage = code_b_min_coverage(cfg)))
}

# Shared identity for forecast and Monte Carlo checkpoints.
code_b_fingerprint <- function(input) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  files <- file.path(getOption("code_b.source_root", "."), "R", c(
    "code_a_config.R", "covariance.R", "estimators.R", "time_reversibility.R",
    "data_processing.R", "monte_carlo.R", "code_b_config.R", "forecast_mfbm.R",
    "forecast_har.R", "rolling_window.R", "metrics.R"))
  if (!all(file.exists(files))) stop("cannot fingerprint model source files")
  function_names <- unlist(lapply(files, function(f) {
    expr <- as.list(parse(f))
    unlist(lapply(expr, function(e) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && length(e) == 3L &&
          is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) as.character(e[[2]]) else NULL
    }))
  }))
  # Deparse omits source-reference environments, keeping identities stable across
  # clean sessions while still including runtime function-body/formal changes.
  functions <- lapply(sort(unique(function_names)), function(nm) {
    fun <- get(nm, mode = "function")
    list(formals = deparse(formals(fun)), body = deparse(body(fun)))
  })
  saveRDS(list(input = input, functions = functions,
    source_md5 = unname(tools::md5sum(files))), tmp, version = 2)
  unname(tools::md5sum(tmp))
}

# Match entire origin/target/h keys before comparing different window lengths.
align_forecast_keys <- function(reference, comparison) {
  check_forecasts(reference); check_forecasts(comparison)
  key <- function(x) paste(x$h, x$origin_date, x$target_date)
  common <- intersect(key(reference), key(comparison))
  if (!length(common)) stop("no common forecast keys")
  a <- reference[key(reference) %in% common, , drop = FALSE]
  b <- comparison[key(comparison) %in% common, , drop = FALSE]
  rownames(a) <- rownames(b) <- NULL
  list(reference = a, comparison = b)
}
