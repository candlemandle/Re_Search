# Research_Seminars_DSBA2027

## Rough volatility with multivariate fBm — replication

Replication of Bibinger, Yu and Zhang (2026), *Modeling and Forecasting Realized Volatility
with Multivariate Fractional Brownian Motion*, JBES (arXiv:2504.15985).

> This README currently documents **Code A** (data, mfBm covariance structure, estimators,
> time-reversibility test, Monte Carlo). Code B (forecasting) and Code C (extension,
> `run_all`, notebooks) add their sections.

## Requirements

* R ≥ 4.3 (developed with R 4.6.1)
* CRAN package `testthat` (tests only). Code A uses base R otherwise (`stats`, `parallel`,
  `grDevices`, `graphics`, `utils`).

```bash
Rscript -e 'install.packages("testthat")'
```

## Code A: commands (run from the project root)

| step | command | runtime (8-core laptop) |
|---|---|---|
| 0. compare with the authors' code | `Rscript scripts/00_check_author_code.R` | ~10 s |
| 1. download Risk Lab data (~20 MB) | `Rscript scripts/01_download_data.R` | ~2 min |
| 2. build log-volatility panels | `Rscript scripts/02_prepare_data.R` | ~2 s |
| 3. estimates, test, Figures 1–3, 7–8 | `Rscript scripts/03_estimate_parameters.R full` | ~15 s |
| 4. Monte Carlo, smoke (outputs in `results/smoke/`) | `Rscript scripts/04_run_monte_carlo.R smoke` | ~1.5 min |
| 4. Monte Carlo, full | `Rscript scripts/04_run_monte_carlo.R full` | ~15–25 min |
| tests | `Rscript -e 'testthat::test_dir("tests/testthat")'` | ~10 s |

`scripts/04_run_monte_carlo.R full t7` runs a single experiment (`e1`, `t6`, `t7`, `t12`, `g`).
The number of worker processes is `detectCores() - 1`, or `MFBM_CORES`. Results do not depend
on it (fixed seed per chunk of 50 replications). The script relaunches itself with
`OPENBLAS_NUM_THREADS=1` to avoid BLAS thread oversubscription in forked workers.

## Code A: layout

```text
R/code_a_config.R        settings (smoke/full), seeds, tickers, logging helpers
R/covariance.R           mfBm covariances (1), (4), (5), (19), rho_max, admissibility,
                         exact simulation by multivariate circulant embedding
R/estimators.R           MM estimators (6)-(8b), asymptotic variances (9), (10), (14), (27),
                         Isserlis cross-checks, Appendix G GMM, Amblard-Coeurjolly estimator
R/time_reversibility.R   test of eta = 0, eq. (11)
R/data_processing.R      Risk Lab download/parse, panels, panel estimation, rolling H
R/monte_carlo.R          Appendix E.1, E.2, E.4 and G experiments
scripts/00-04            pipeline above
author_code/BYZ/         authors' R files (BYZ.zip), unchanged
docs/paper_values/       numbers transcribed from the paper for automatic comparison
docs/paper_coverage.csv  every table/figure of the paper -> script, output, status
data/README.md           data source, ticker mapping, sample rules
```

Outputs: `results/tables/code_a_*.csv`, `results/figures/code_a_*.pdf`, log `results/logs/code_a.log`.
Every Monte Carlo table has a `*_vs_paper.csv` companion with the paper value and the difference.

## Interfaces for Code B / C

```r
source("R/code_a_config.R"); source_code_a()
panel <- load_logvol_panel("dj30")                      # date + 30 log-vol columns
inc   <- panel_increments(panel, c("AAPL", "ALD"))      # increments on common dates
est   <- estimate_mfbm(inc$X, delta = 1/252)            # H, sigma2, rho, eta
par   <- mfbm_params(est$H, est$rho, 0, sqrt(est$sigma2))
S     <- mfbm_level_cov_matrix(par, times = (1:500) / 252)  # Sigma_{t,d} for forecast (15)
X     <- simulate_mfbm(505, c(0.1, 0.4), 0.4, delta = 1/250, nsim = 100)  # increments
```
