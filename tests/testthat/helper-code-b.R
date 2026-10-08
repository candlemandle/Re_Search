# load code B from the project root (code A is loaded by helper-code-a.R before this file)
code_b_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
options(code_b.source_root = code_b_root)
for (f in c("code_b_config.R", "forecast_mfbm.R", "forecast_har.R", "rolling_window.R", "metrics.R")) {
  source(file.path(code_b_root, "R", f))
}
