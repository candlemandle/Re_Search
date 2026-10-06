# PR: Code B: optimal forecasts, rolling forecasts, HAR family, DJ30 and Mag7 results

Reviewer: participant 5. Branch: `Amir` (built on top of `code-a`).

## Команды запуска

From the project root, in a clean R session, after Code A steps 01–02 (the processed panels
`data/processed/{dj30,mag7}_logvol.csv` are already in the repository):

```bash
Rscript scripts/05_run_forecast_simulations.R smoke   # ~2 min, writes to results/smoke/
Rscript scripts/06_run_empirical_forecasts.R smoke    # ~20 s, every 40th day
Rscript scripts/07_run_robustness.R smoke             # ~15 s, needs 06 smoke first

Rscript scripts/06_run_empirical_forecasts.R full     # ~13 min on 11 cores
Rscript scripts/07_run_robustness.R full              # ~9 min, needs 06 full first
Rscript scripts/05_run_forecast_simulations.R full    # ~14 min
Rscript -e 'testthat::test_dir("tests/testthat")'     # ~1 min, Code A + Code B tests
```

Workers: `MFBM_CORES` (default `detectCores() - 1`), results do not depend on it. Like Code A,
every script restarts itself once with single-threaded BLAS (forked workers + multithreaded
OpenBLAS made each window ~7x slower).

Checkpoints: every finished model (06, 07) and every simulation experiment (05) is saved under
`results/forecasts/checkpoints/` and `results/simulations/checkpoints/`. A rerun loads them, so an
interrupted run continues where it stopped. **Delete these folders after changing model code.**

## Модели и периоды

Target: RV of AAPL = Risk Lab `qmle_trade` annualized volatility (the series Code A identified
by matching Table 1 exactly). Horizons h = 1, 2, 3, 4, 5, 10, 15, 20 days, direct forecasts of
the single day t + h. All models use the same days (where all 5 assets are observed), the same
forecast origins and the same rolling window of **500 days**.

| class | models (assets added in the paper order) | how |
|---|---|---|
| mfBm | fBm, bfBm, mfBm3, mfBm4, mfBm5 | estimate H, sigma, rho on the window (Code A MM estimators, eta = 0), conditional mean (15) of log RV given all levels of the window, RV forecast exp(mean + var / 2) |
| VHAR | HAR, VHAR2–5 | OLS of RV_{t+h} on daily / weekly / monthly RV of all assets, (33)–(34) |
| VHARF | HAR, VHARF2–5 | own HAR terms + HAR terms of the average of the other assets, (35) |
| VLHAR | LHAR, VLHAR2–5 | (34) in log RV, forecast exp(fit + s² / 2) |

DJ30 assets: AAPL, ALD, AMGN, AXP, BA. Mag7 assets: AAPL, AMZN, FB, GOOG, MSFT.

Periods (by target date): DJ30 full 2007-01-03 – 2025-01-14 (4508 forecast days),
period 1 2007-01-03 – 2017-04-11, period 2 2017-04-12 – 2021-07-30; Mag7 2016-03-29 – 2025-01-14 (2188 days).

Metrics: MSFE × 100, QLIKE = RV/F − log(RV/F) − 1, MCS p-values (Hansen, Lunde and Nason 2011,
T_max statistic, circular block bootstrap, block 20, 5000 replications), one MCS per model class
and horizon, as in the paper.

## Созданные outputs

Tables (`results/tables/`):

| file | paper |
|---|---|
| `code_b_fig04_05_weights_msfe.csv` | data of Figures 4–5 |
| `code_b_fig06_dimension.csv`, `code_b_fig06_limits.csv` | data of Figure 6 and the limits d → ∞ |
| `code_b_tables08_11_simulations.csv`, `code_b_tables08_11_vs_paper.csv` | Tables 8–11 |
| `code_b_table02_dj30.csv` | Table 2 (3 periods) |
| `code_b_table03_dj30.csv`, `code_b_table15_dj30.csv`, `code_b_table16_dj30.csv` | Tables 3, 15, 16 |
| `code_b_table18_dj30.csv` | Table 18 (QLIKE) |
| `code_b_dj30_tables_vs_paper.csv`, `code_b_dj30_original_vs_ours.csv` | all DJ30 tables vs paper |
| `code_b_mag7_tables_vs_paper.csv`, `code_b_mag7_original_vs_ours.csv` | Tables 17, 19 vs paper |
| `code_b_robust_dj30_subperiods.csv`, `code_b_robust_dj30_window250.csv` | extra robustness (not in the paper) |
| `code_b_checks_{dj30,mag7,dj30_window250}.csv` | automatic checks (same dates, look-ahead, NA/Inf) |

`*_vs_paper.csv` has our value, the paper value, the difference and the MCS p-values next to each
other. `*_original_vs_ours.csv` is the short version: loss of the 5-asset model / univariate model
for every table, panel and h, ours vs paper, and whether the direction agrees.

Figures (`results/figures/`): `code_b_fig04_weights.pdf`, `code_b_fig05_rel_msfe.pdf`,
`code_b_fig06_dimension.pdf`, `code_b_fig_dj30_msfe_ratio.pdf`, `code_b_fig_dj30_class_ratio.pdf`,
`code_b_fig_dj30_aapl_forecasts.pdf`, `code_b_fig_mag7_class_ratio.pdf`.

Raw forecasts: `results/forecasts/code_b_forecasts_{dj30,mag7}.rds` (one row per model, origin, h).
Log: `results/logs/code_b.log`. Paper numbers used for the comparison:
`docs/paper_values/code_b_forecast_tables.csv`, `docs/paper_values/code_b_simulation_tables.csv`.
Coverage rows: `docs/paper_coverage_code_b.csv`.

## Runtime

11 workers on an M-series MacBook: 06 full ~13 min (mfBm5 alone 6.5 min, 0.36 s per window),
07 full ~9 min, 05 full ~14 min (Table 11 dominates: parameters re-estimated and a
1000 × 1000 system solved on each of 10 000 paths), MCS for all tables ~1 min. Dependencies:
R ≥ 4.3, `testthat`; everything else is base R (`stats`, `parallel`).

## Что совпало с оригиналом

* **Table 2 (DJ30, mfBm), the main result.** Full period and period 1: at h = 1 all five MSFEs
  equal the paper to 4 decimals (e.g. fBm 0.5003, mfBm5 0.4979; period 1 fBm 0.6349, mfBm5 0.6301).
  For h ≥ 2 within 0.7%. The gain of mfBm5 over fBm equals the paper within 0.01 percentage points
  at every h in all three periods (full h = 5: 1.28% vs 1.27%; period 1 h = 5: 1.79% vs 1.79%;
  h = 20: 2.30% vs 2.30% and 3.21% vs 3.20%; period 2: about 0, as in the paper).
  MCS: mfBm5 has p = 1 at every h, as in the paper; correlation of all our MCS p-values with the
  paper's 0.94 (e.g. fBm h = 1 0.385 vs 0.402).
* **Tables 3, 15, 16:** HAR, VHAR, VHARF and VLHAR mean |relative difference| 1.8%, 1.2%, 0.6%.
  Same conclusion: every vector HAR is worse than the univariate HAR and gets worse with every
  added asset; VHARF2 = VHAR2 exactly (as in the paper).
* **Tables 18 (QLIKE):** mfBm within 0.3%, VLHAR within 2.4%.
* **Tables 17, 19 (Mag7):** the 5-asset vs univariate comparison has the same direction as the
  paper in 62 of 64 cells; mfBm improves on fBm at long horizons (h = 20: −4.8% vs −5.4% in the
  paper) and all vector HAR classes lose.
* **Section 4:** Figures 4–5 from the general conditional mean equal the closed form of
  Prop. 4.1 (unit test); Figure 6 limits 0.891 and 0.774 vs "89%" and "77.4%" in the text.
* **Tables 8–11, theoretical RMSFEs (brackets): exact.** All 195 compared bracketed values equal ours to
  4 decimals (max |diff| 5e-5 = rounding), including mfBm3/mfBm4 of Table 10.
* **Tables 8–11, Monte Carlo (10 000 reps):** our simulated RMSFEs are within 2 MC standard errors
  of the theory in 96% of cells (mean |rel. diff| 0.53%); the paper's own Monte Carlo is as close
  to the theory (96%, 0.58%). All conclusions hold: no gain for rho = 0 or H1 = H2, gains grow
  with |rho| (−12% at rho = 0.8, h = 1), with H2 − H1 and with more components (Table 10:
  0.4802 → 0.4756 → 0.4686 → 0.4563); estimation error raises RMSFE by ~0.2% and bfBm still beats
  fBm in all 10 cells with rho = 0.4 (Table 11).

## Что отличается и почему

1. **Table 2, period 2 levels.** Our MSFEs are −0.5% to +9% off (fBm h = 1 0.4100 vs 0.3842),
   while the model ratios still match (mfBm5 / fBm 1.0027 vs 1.0034 at h = 1). Full period and
   period 1 match exactly with the same code, so the paper's period-2 sample must differ from the
   printed dates; no window of the same length reproduces all three printed values at once.
2. **Window = 500 days, not 504.** The paper says "two-year rolling window". Both samples have
   exactly 500 common days before the paper's first forecast date (2007-01-03 and 2016-03-29), and
   with 500 the h = 1 numbers of Table 2 are exact.
3. **VHAR / VHARF at long horizons and their QLIKE.** With 3 × 5 regressors on 500 days the
   direct 15–20 day regressions are unstable; small differences in the fit give large differences
   in MSFE and especially in QLIKE (which explodes when a forecast is close to 0). A linear HAR
   forecast ≤ 0 (3 of 684 228 DJ30 forecasts) is replaced by the smallest RV of the window,
   otherwise QLIKE is undefined. The paper does not say what it does.
4. **Mag7 levels ~4% higher.** Same code as DJ30; the Mag7 common sample has 10 days with a missing
   asset (AMZN 2, GOOG 7, MSFT 1) which we drop for every model. Conclusions unchanged.
5. **Figure 6 parameters.** The text says H1 = 0.1, H = 0.4, the caption H1 = 0.4, H = 0.1; only
   the caption reproduces the quoted 89% and 77.4% (text order gives 91.4% and 77.9%).
6. **Admissibility.** If the estimated rho on a window are not admissible together with the
   estimated H (no mfBm exists, Amblard et al. 2013 condition), rho is shrunk by 5% steps until it
   is. The paper does not discuss this.
7. **MCS implementation.** Our own ~30 lines (T_max, circular block bootstrap) instead of a
   package; Appendix F.4's description of T_max as a distance between empirical distribution
   functions does not match Hansen et al. (2011), we follow Hansen et al.
8. **Simulations vs the paper's Monte Carlo.** Different random numbers (Code A's exact circulant
   embedding, not the authors' simMFBM), so cells differ by MC noise: 97% within 3 MC standard
   errors of the paper's value; the largest gaps (Table 9, H = (0.1, 0.1), h = 3: 0.5198 vs 0.5348)
   are cells where the paper's value is itself ~2 SE from the theory. Note the MC errors of different
   h and of fBm/bfBm on the same paths are correlated.

## Проверка look-ahead

* A model receives only `L[1:t, ]` (all data up to the forecast origin t); HAR regressions use
  pairs (s, s + h) with s + h ≤ t. Test `forecasts do not change when the future is changed`
  shifts all data after t by +5 and checks the forecasts are identical for every model of all
  four classes (11 models on a 3-asset test panel).
* `check_forecasts()` runs on every full run and stops on failure: `origin_date < target_date`
  for every forecast, the same (h, target date) set for every model, no NA/Inf, all forecasts > 0.
  Results: `results/tables/code_b_checks_*.csv` (all TRUE).
* Parameters of each window are estimated on that window only.
* RMSFE / MSFE checked by hand on a small example (test-metrics.R).

## Какие тесты проходят

`tests/testthat/test-forecast-mfbm.R`, `test-metrics.R`, `test-rolling-window.R`: 69 expectations,
all pass, together with the 113 of Code A. They cover: theoretical RMSFEs = Table 8 brackets,
general forecast = closed forms of Prop. 4.1 and 4.3, no gain for equal H (Prop. B.1) or rho = 0,
admissibility vs rho_max, shrinking to admissible parameters, the window forecast written out
by hand from eq. (15), simulated vs theoretical RMSFE, QLIKE and MSFE by hand, MCS behaviour and
reproducibility, block bootstrap, HAR features, VHARF2 = VHAR2, no look-ahead, same test dates.

## Merge notes

* The branch needs Code A (`R/code_a_config.R`, `R/covariance.R`, `R/estimators.R`,
  `R/data_processing.R`, data panels). Merge `code-a` first.
* After merging, append to `.gitignore`:
  `results/forecasts/checkpoints/` and `results/simulations/`.
* Add the rows of `docs/paper_coverage_code_b.csv` to `docs/paper_coverage.csv` in place of the
  `B,not_started` rows (kept in a separate file to avoid a merge conflict).
* README section: see below.

### README section (to paste after the Code A section)

```markdown
## Code B: forecasts (run from the project root)

| step | command | runtime (11 workers) |
|---|---|---|
| 5. Figures 4–6, Tables 8–11 | `Rscript scripts/05_run_forecast_simulations.R full` | ~14 min |
| 6. DJ30 forecasts, Tables 2, 3, 15, 16, 18 | `Rscript scripts/06_run_empirical_forecasts.R full` | ~13 min |
| 7. Mag7 (Tables 17, 19) and robustness | `Rscript scripts/07_run_robustness.R full` | ~9 min |

`smoke` instead of `full` runs the same code on fewer days / replications into `results/smoke/`.
Finished models are checkpointed; delete `results/*/checkpoints/` after changing model code.

R/code_b_config.R     settings (smoke/full), logging, single-thread BLAS restart
R/forecast_mfbm.R     optimal forecast (15), Prop. 4.1-4.3, rolling-window mfBm forecast, E.3 Monte Carlo
R/forecast_har.R      HAR, VHAR, VHARF, LHAR, VLHAR (Appendix F.2)
R/rolling_window.R    all models on the same forecast origins, checks
R/metrics.R           MSFE, QLIKE, model confidence set, tables vs paper
```
