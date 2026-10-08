if (dir.exists(".Rlib")) .libPaths(c(normalizePath(".Rlib"), .libPaths()))
if (!requireNamespace("testthat", quietly = TRUE)) stop("Run Rscript scripts/setup_environment.R first")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("smoke", "full")) else "smoke"
dir <- file.path(if (mode == "smoke") "results/smoke" else "results", "logs")
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MFBM_CORES = "1")
res <- testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
tab <- as.data.frame(res)
utils::write.csv(tab[, setdiff(names(tab), "result"), drop = FALSE],
  file.path(dir, "test_results.csv"), row.names = FALSE)
if (any(tab$failed > 0 | tab$error)) stop("Tests failed: see results logs")
cat("All", nrow(tab), "test blocks passed;", sum(tab$passed), "expectations.\n")
