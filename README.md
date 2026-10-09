# Rough-volatility replication and mfOU extension

This repository reproduces and extends Bibinger, Yu and Zhang (2026),
*Modeling and Forecasting Realized Volatility with Multivariate Fractional Brownian
Motion* (JBES; arXiv:2504.15985).

- **Code A** reproduces the mfBm covariance structure, moment estimators,
  time-reversibility test, empirical parameter estimates and Monte Carlo exercises.
- **Code B** reproduces theoretical forecasting figures and implements rolling
  mfBm/HAR forecasts on the authors' DJ30 and Mag7 panels.
- **Code C** is an explicitly exploratory extension: stationary mfOU forecasting,
  window-length sensitivity, known-versus-estimated simulation, and asset selection.

Russian overview: [README_RU.md](README_RU.md). Detailed audit:
[docs/FINAL_AUDIT.md](docs/FINAL_AUDIT.md).

## What is reproduced, and what is not

| Part | Status | Evidence |
|---|---|---|
| Table 1 point estimates | exact to four decimals, 25/25 cells | `results/tables/code_a_table01_estimates.csv` |
| Table 1 standard errors | 23/25 exact to four decimals; two differ by 0.0001 because of series truncation | `results/tables/code_a_table01_se_check.csv` |
| Appendix Monte Carlo | 290/316 cells within two combined Monte Carlo standard errors | `results/tables/code_a_mc_agreement_summary.csv` |
| Figures 1–7 and forecast formulas | reproduced/implemented | `results/figures/`, tests, `docs/paper_coverage.csv` |
| Table 2 forecasting pattern | close for the full sample and period 1 | `results/tables/code_b_dj30_tables_vs_paper*.csv` |
| Table 2 period-2 levels | **not an exact replication**; most losses remain about 6–9% above the paper although model ratios agree | `docs/response_to_review_B.md` |
| Table 3 / vector HAR levels | implementation is reproducible, but positivity filtering makes some QLIKE results unstable | `docs/response_to_review_B.md` |
| mfOU extension | new research result, not a result reported in the paper | `code_c/RESULTS_EN.md` |

The primary missing-day convention is the **trading calendar**: missing observations
retain their positions in time. `common_obs` outputs are preserved as a historical
sensitivity calculation. The paper does not completely specify this rule, so both
conventions are reported rather than silently choosing the one closest to a table.

## Headline findings

1. Adding related assets improves mfBm forecasts gradually, especially at longer
   horizons, but it does not make mfBm the best model in every comparison. In the
   full DJ30 sample, mfBm5 lowers MSFE relative to univariate fBm by about 0.48% at
   `h=1`, 1.23% at `h=5`, and 2.24% at `h=20`.
2. A 250-day window is materially worse than the 500-day baseline for the five-asset
   mfBm at long horizons; the matched `h=20` MSFE ratio is about 1.32. Estimating a
   5×5 dependence structure needs enough history.
3. On Mag7, five-asset mfBm is worse at very short horizons but improves at longer
   horizons; the crossover occurs near `h=4` in this implementation.
4. The mfOU extension is promising only at longer horizons in the executed daily
   study. For AAPL, 2019–2021, `h=20`, mfOU reduces MSFE by about 12.2% and QLIKE by
   about 9.2% versus mfBm (750 paired forecasts; exploratory block-bootstrap intervals
   exclude zero). At `h=1`, the result is mixed.
5. The H-gap asset-selection heuristic is not a universal winner: for AAPL it loses
   at `h=1` and gains about 5–6% at `h=20`, based on only 12 quarterly origins.

These are forecasting results, not trading-profit or causal claims. Non-rejection of
time reversibility is not proof that the asymmetry parameter is zero, and a positive
fitted mfOU mean-reversion rate is not by itself evidence of true mean reversion.

## Requirements and setup

- R 4.3 or newer (the saved runs used R 4.5–4.6)
- CRAN package `testthat` for tests
- base R is sufficient for the estimation and forecasting pipelines

From the repository root:

```bash
Rscript scripts/setup_environment.R
Rscript scripts/run_tests.R
Rscript code_c/scripts/run_tests.R
```

The processed data snapshot is committed. Its expected MD5 values are recorded in
`data/processed/MD5SUMS`. To rebuild it from Risk Lab instead, run
`Rscript scripts/01_download_data.R` and `Rscript scripts/02_prepare_data.R`.

## Reproduction commands

Quick end-to-end check:

```bash
Rscript scripts/run_all.R smoke
```

Paper pipeline using the committed data:

```bash
Rscript scripts/run_all.R paper
```

This runs the non-Monte-Carlo Code A checks and full Code B forecasts. The full Code A
Monte Carlo is deliberately separate because it is expensive:

```bash
Rscript scripts/04_run_monte_carlo.R full
Rscript scripts/04b_compare_mc_with_paper.R full
```

Code C examples and executed-study settings:

```bash
Rscript code_c/scripts/run_experiments.R smoke all
Rscript code_c/scripts/run_experiments.R study mfou --dimensions=1,2 --sensitivity=none --study-step=1 --output=results/code_c_v2/dense_study
Rscript code_c/scripts/run_experiments.R study mfou --study-step=63 --output=results/code_c_v2/window_study
Rscript code_c/scripts/run_experiments.R study selection --study-step=63 --output=results/code_c_v2/selection_study
Rscript code_c/scripts/run_experiments.R study simulation --output=results/code_c_v2/simulation_study
```

The complete daily Code C grid over every dimension and every window was **not** run;
the included outputs distinguish the dense two-dimensional study from sparse
window/selection studies. Do not present the 12-origin results as confirmatory tests.

## Repository map

```text
R/                         Code A/B implementation
scripts/                   A/B pipelines and integrated runners
tests/testthat/             A/B tests
code_c/                    additive mfOU and selection extension
data/processed/            pinned log-volatility panels and checksums
results/tables/             generated numerical results
results/figures/            generated paper/forecast figures
docs/paper_values/          values transcribed from the article
docs/paper_coverage.csv     table/figure-by-table/figure coverage
docs/FINAL_AUDIT.md         claims, discrepancies and limitations
docs/PRESENTATION_GUIDE.md  which real figures to use in the talk
```

## Reproducibility notes

- All rolling forecasts are trained only on observations available at the origin.
- Forecast comparisons use identical target keys and realized outcomes.
- Code B records calendar, positivity, uniqueness and no-look-ahead checks.
- Code C records source/data hashes, failures, parameter boundaries and numerical
  diagnostics; its saved baseline hashes match the integrated A/B sources.
- `eta` has two sign conventions in the printed paper/code. The implementation exposes
  both; the time-reversibility test is unaffected by the sign choice.
- Model Confidence Set p-values are random-bootstrap quantities. Exact table values
  depend on the seed; the implementation was independently compared with NumPy and
  the `arch` package in `docs/review_B/`.

For presentation-ready assets and the interpretation to put next to each one, see
[docs/PRESENTATION_GUIDE.md](docs/PRESENTATION_GUIDE.md).
