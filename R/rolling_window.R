# rolling-window forecasting pipeline: one function runs every model on the
# same forecast origins, so all models have the same test dates, target and h.

# all models of the paper for one asset set (assets[1] is the series we forecast).
# class is the group used for the MCS p-values (tables 2, 3, 15, 16, 17)
code_b_models <- function(assets, window, delta) {
  model <- function(name, class, k, fun) list(name = name, class = class, assets = assets[1:k], fun = fun)
  suffix <- c("", "2", "3", "4", "5")
  out <- list()
  for (k in seq_len(min(length(assets), 5))) {
    # local() so every closure keeps its own k
    out <- c(out, local({
      k <- k
      list(
        model(c("fBm", "bfBm", "mfBm3", "mfBm4", "mfBm5")[k], "mfbm", k,
              function(L, h) forecast_mfbm_window(utils::tail(L, window), h, delta)),
        model(paste0(if (k == 1) "HAR" else "VHAR", suffix[k]), "vhar", k,
              function(L, h) forecast_har_window(L, h, window, "vhar")),
        model(paste0(if (k == 1) "HAR" else "VHARF", suffix[k]), "vharf", k,
              function(L, h) forecast_har_window(L, h, window, "vharf")),
        model(paste0(if (k == 1) "LHAR" else "VLHAR", suffix[k]), "vlhar", k,
              function(L, h) forecast_har_window(L, h, window, "vhar", log = TRUE))
      )
    }))
  }
  # HAR sits in two classes (tables 3 and 15). it is the same model so we keep one copy
  names(out) <- vapply(out, `[[`, "", "name")
  out[!duplicated(names(out))]
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
#   L     : log vol matrix (rows = common dates of all assets)
#   dates : dates of the rows of L
# the model only ever sees L[1:t, ] so it can not use the future
rolling_forecasts <- function(L, dates, model, origins, horizons, cores = 1L) {
  one_origin <- function(t) {
    past <- L[1:t, model$assets, drop = FALSE]
    fc <- model$fun(past, horizons)
    ok <- t + horizons <= nrow(L)
    data.frame(
      model = model$name,
      class = model$class,
      h = horizons[ok],
      origin_date = dates[t],
      target_date = dates[t + horizons[ok]],
      forecast = fc[ok],
      actual = exp(L[t + horizons[ok], 1])
    )
  }
  res <- parallel::mclapply(origins, one_origin, mc.cores = cores)
  failed <- vapply(res, inherits, TRUE, "try-error")
  if (any(failed)) stop(model$name, ": ", res[[which(failed)[1]]])
  do.call(rbind, res)
}

# run every model on the same origins and stack the results.
# only the days where all assets are observed are used, for every model.
# model_names picks a subset of models (NULL = all).
# checkpoint: if a name is given every finished model is saved to
# <results>/forecasts/checkpoints/<name>_<model>.rds and a rerun loads it
# instead of computing it again (useful if the laptop sleeps or the run dies).
# delete that folder after changing the model code
run_all_models <- function(panel, assets, cfg, cores = 1L, model_names = NULL,
                           checkpoint = NULL) {
  sub <- panel[stats::complete.cases(panel[, assets]), c("date", assets)]
  L <- as.matrix(sub[, assets])
  origins <- forecast_origins(nrow(L), cfg$window, cfg$origin_step)
  models <- code_b_models(assets, cfg$window, cfg$delta)
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
  fc
}

# checks asked for in the task. returns a one-row data.frame and stops on a hard error
check_forecasts <- function(fc) {
  key <- paste(fc$h, fc$target_date)
  same_dates <- all(tapply(key, fc$model, function(k) setequal(k, key[fc$model == fc$model[1]])))
  res <- data.frame(
    n_models = length(unique(fc$model)),
    n_rows = nrow(fc),
    same_test_dates = same_dates,
    no_lookahead = all(fc$origin_date < fc$target_date),
    no_na_inf = all(is.finite(fc$forecast) & is.finite(fc$actual)),
    positive_forecasts = all(fc$forecast > 0)
  )
  if (!res$same_test_dates || !res$no_lookahead) stop("forecast checks failed")
  res
}
