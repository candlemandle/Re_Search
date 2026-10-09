#!/usr/bin/env Rscript

# The scientific pipeline uses base R. testthat is the only package required by
# the automated test suites.
required <- c("testthat")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  repos <- getOption("repos")
  if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
    repos["CRAN"] <- "https://cloud.r-project.org"
  }
  install.packages(missing, repos = repos)
}

still_missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Failed to install required packages: ", paste(still_missing, collapse = ", "))
}

cat("Environment ready. R:", R.version.string, "\n")
cat("testthat:", as.character(utils::packageVersion("testthat")), "\n")
