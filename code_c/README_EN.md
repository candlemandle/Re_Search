# Code C: additional realized-volatility forecasting experiments

Code C extends reviewed A/B v2 with stationary mfOU forecasting, training-window sensitivity, known/estimated-parameter simulations and asset-selection experiments. Its purpose is to measure possible improvements and their limits. A better result is never assumed by the implementation.

Read [the simple explanation](EXPLAINED_EN.md), then [executed outcomes](RESULTS_EN.md). [Methods and references](docs/METHODS_EN.md) contain the mathematical specification. [The branch guide](COMMIT_GUIDE_EN.md) explains which files to upload.

## Integration with the team project

B v2 contains A v2. C sources `R/code_b_config.R` and calls `source_code_b(root)`. The baseline is a direct call to `forecast_mfbm_calendar()`, with B's parameter utilities, admissibility checks, trading calendar, minimum coverage and loss functions. The reviewed baseline is not reimplemented.

All additions are under `code_c/`, and functions use the `c2_` prefix. C writes its own outputs, does not edit A/B files, and has a separate entry point. No modification of a shared team runner is needed. In the merged repository:

```text
repository/
  R/                       # reviewed A/B
  data/processed/           # colleagues' existing data
  scripts/                 # colleagues' entry points
  code_c/                  # C-owned overlay
  results/code_c_v2/        # generated C outcomes
```

Locally, the loader selects `part B version 2` when present; otherwise it uses the repository root. Override this with `--ab-root=PATH` or `CODE_C_AB_ROOT`. An older B without the trading-calendar API is rejected with an explanatory error.

## Research questions

| Experiment | Question | Outputs and basis |
|---|---|---|
| `mfou` | Does stationary mfOU improve on mfBm with identical information? | Paired MSFE/RMSFE/QLIKE, gains, parameter diagnostics, cumulative losses; article §6 |
| 250/500/1000-day sensitivity | Does the amount of training history change the answer? | Matched-key window comparisons; supplied research task, extending the article's two-year window |
| `simulation` | What changes when model parameters have to be estimated? | Paired losses/MCSE, known/estimated forecasts, parameter bias/RMSE; extension of Appendix E.3/Table 11 |
| `selection` | Does choosing other additional assets help? | Fixed, correlation and H-gap-guided rules, with fBm/HAR benchmarks; motivated by the article's dependence/Hurst-gap theory |

The score `abs(rho)*abs(H_target-H_asset)` is a heuristic, not a proved optimal selector. Exact definitions, departures from published estimators, and source anchors are in the methods document.

## Input data and units

C uses `data/processed/dj30_logvol.csv` in the selected A/B project. Its first column is `date` (`YYYY-MM-DD`); the remaining 30 columns contain numeric **log RV**, with `NA` for missing observations. The supplied snapshot contains 5019 trading rows, 2005-01-05 through 2025-01-14. Historical identifiers such as `ALD` are preserved.

The original measure is annualized `qmle_trade` volatility. Return to its units with `exp(log_RV)`; do not square it. The team also supplies `mag7_logvol.csv`, with 2698 rows and five available series (AAPL/AMZN/FB/GOOG/MSFT); these C experiments use DJ30. No new download is required.

The main target is AAPL. Nested sets add ALD, AMGN, AXP and BA, using a 500-trading-row window and horizons `1,2,3,4,5,10,15,20`. Empirical time uses `1/252` years per trading day; simulation explicitly uses `1/250` for its separate setup.

Missing days do not compress the time axis. All active nested sets share the largest-set eligibility calendar within each training window, as in B: assets must be observed at the origin and window coverage must be at least 90%. A missing target removes the corresponding outcome jointly. Future values never enter fitting, ranking or calendar eligibility.

## Quick start

Use R and the team's existing environment. Verified with R 4.5.3 and `testthat` 3.3.2. Experiments use base R and A/B functions; tests additionally require `testthat`. Install it in your usual R library if absent. C adds no Quarto, Pandoc or Python runtime dependency.

Run **from the repository root**:

```sh
Rscript code_c/scripts/run_tests.R
Rscript code_c/scripts/run_experiments.R smoke all
```

For an explicit local path:

```sh
Rscript code_c/scripts/run_experiments.R smoke all "--ab-root=part B version 2"
```

Preview full-workload counts and memory before computing forecasts:

```sh
Rscript code_c/scripts/run_experiments.R full mfou --plan
```

Settings used for the recorded studies:

```sh
# Approximately quarterly origins; all nested sets and three windows
Rscript code_c/scripts/run_experiments.R study mfou --study-step=63 --output=results/code_c_v2/study

# Daily forecasts, one/two assets, 2019–2021
Rscript code_c/scripts/run_experiments.R study mfou --dimensions=1,2 --sensitivity=none --study-step=1 --output=results/code_c_v2/dense_study

# Thirty independent replications per generator
Rscript code_c/scripts/run_experiments.R study simulation --output=results/code_c_v2/study_simulation

# Asset selection for AAPL/MSFT/JPM on the quarterly grid
Rscript code_c/scripts/run_experiments.R study selection --study-step=63 --output=results/code_c_v2/study_selection

Rscript code_c/scripts/compare_blocks.R results/code_c_v2/dense_study/mfou/empirical_run.rds
Rscript code_c/scripts/verify_numerics.R results/code_c_v2/study/mfou/empirical_run.rds
```

Default `study` origins span 2019-01-02 through 2021-12-31 at a 21-trading-day step. Explicit options: `--study-start`, `--study-end`, `--study-step`, `--dimensions=1,2`, `--sensitivity=none` or `250,1000`, `--max-origins=1`, `--output=PATH`. In `all`, the origin settings also apply to selection. The preflight plan is produced before the `--max-origins` execution limit.

`smoke` uses three sparse origins and four simulation replications. `full` uses the configured full empirical period and 200 replications. Full daily empirical execution can take days. A five-asset, 1000-row covariance has up to 5000×5000 entries: about 191 MiB for the matrix, with a conservative working-memory estimate of 1144 MiB (1.12 GiB). Hardware affects runtime. The entry point uses one BLAS worker. The full daily all-window study has not been executed here.

## Checkpoints and provenance

Cells are saved under `checkpoints/`. Repeating exactly the same command and output path resumes completed cells. Data, configuration or R-source changes invalidate the fingerprint and prevent inappropriate reuse. Tables/reports are written after a stage finishes; a partial cache is not a completed research result.

`reanalyze.R` refreshes tables and figures without refitting forecasts. For older per-set C calendars it applies B's joint past-only eligibility rule, first preserving the original RDS with suffix `.original_per_set.rds`. Retained forecasts remain identical. The analysis fingerprint is saved separately from original forecasting provenance; reanalysis does not claim that the final source version refitted the entire historical study.

## Output guide

| Output | Meaning |
|---|---|
| `mfou/REPORT_EN.md`, `REPORT_RU.md` | Run scope, actual dates, summary and limits |
| `mfou/forecasts.csv` | Model, dimension, window, origin/target dates, forecast/actual |
| `mfou/metrics.csv`, `comparisons.csv` | Absolute losses, paired gains, uncertainty and coverage |
| `mfou/window_sensitivity.csv` | Alternative windows on matched outcomes |
| `mfou/parameters.csv` | H, κ, σ, µ, optimization bounds/starts, identification and numerical diagnostics |
| `mfou/origin_audit.csv`, `origin_coverage.csv` | Eligibility and exclusions; origin and horizon counts kept distinct |
| `mfou/estimation_failures.csv` | Failed fits; header-only file indicates zero recorded failures |
| `mfou/cumulative_loss.csv`, `figures/` | Loss accumulation and rolling parameter plots |
| `selection/selected_assets.csv`, `metrics.csv` | Training-only selections and descriptive losses |
| `simulation/comparisons.csv`, `parameter_summary.csv` | Paired losses/MCSE and parameter bias/RMSE |
| `run_manifest.csv`, `configuration.rds`, `session_info.txt` | Completed stages, settings and environment |
| `source_identity.csv`, `data_identity.csv` | R-source and processed-panel identities |

`period=full` denotes the **entire executed sample of that run**, not proof of a full-paper run. Paper subperiods use target dates: through 2017-04-11, and 2017-04-12 through 2021-07-30. Missing subperiods receive no fabricated rows. Origins ending in December 2021 can have 20-day targets in January 2022.

## Small procedural architecture

| Module | Responsibility |
|---|---|
| `R/config.R` | A/B loader, settings, origin/asset scenarios |
| `R/mfou.R` | Stationary covariance, moment fitting, Gaussian forecast |
| `R/experiments.R` | Empirical comparison and joint calendar |
| `R/selection.R` | Additional-asset rules using B predictors |
| `R/simulation.R` | Exact Gaussian data generation, known/estimated outcomes |
| `R/metrics.R` | Losses, pairing, calendar bootstrap and window comparisons |
| `R/outputs.R` | Tables, plots, reports and workload plan |

No service framework, database or class hierarchy is needed. For direct use:

```r
options(code_c.root = normalizePath("code_c"))
for (f in sort(list.files("code_c/R", pattern="[.]R$", full.names=TRUE))) source(f)
ab_root <- c2_load_ab()
panel <- load_logvol_panel("dj30", file.path(ab_root,"data/processed"))
cfg <- c2_config("smoke")
x <- c2_run_empirical(panel, cfg)
c2_write_empirical(x, "results/my_c_experiment")
```

## Interpretation and limits

MSFE is mean squared forecast error; RMSFE is its square root; QLIKE is `r-log(r)-1` for `r=actual/forecast`. Lower is better. `gain_pct=100*(1-loss_mfOU/loss_mfBm)` uses mean losses on the exact same outcomes; positive favors mfOU.

Daily paired series receive exploratory 95% pointwise calendar-block intervals when coverage and sample-size guards pass. Missing pairs retain their trading positions. There is no multiple-testing adjustment. Sparse grids receive no intervals. Simulation MCSE quantifies finite-replication error, not real-data forecast uncertainty.

The mfOU estimator is a custom reduced moment method, not the literature's full joint GMM implementation. True parameters appear only in simulations. `ln(2)/κ` is a deterministic drift half-life, not a fractional-autocorrelation half-life. Plug-in conditional variance excludes parameter-estimation uncertainty. The empirical contrast changes model, estimator and initialization together; gains cannot be causally attributed to mean reversion alone.

A/B were accepted with reproduction limitations; C preserves these caveats and does not claim exact reproduction of published numbers. Final independent acceptance of C is not asserted.
