# Download Risk Lab daily volatilities for the DJ30 (paper tickers) and the
# Magnificent 7 subset used in Appendix F.3.
# Usage (from project root): Rscript scripts/01_download_data.R [--overwrite]

source("R/code_a_config.R")
source_code_a()
cfg <- code_a_config("full")
overwrite <- "--overwrite" %in% commandArgs(trailingOnly = TRUE)

ids <- c(cfg$dj30, cfg$mag7[setdiff(names(cfg$mag7), names(cfg$dj30))])
code_a_log("01_download_data: ", length(ids), " tickers from ", cfg$risklab_url)
for (tk in names(ids)) {
  path <- download_risklab(tk, ids[[tk]], cfg$risklab_url, overwrite = overwrite)
  code_a_log(sprintf("  %-5s permno %-6s -> %s (%.0f KB)", tk, ids[[tk]], path,
                     file.size(path) / 1024))
}
code_a_log("01_download_data: done")
