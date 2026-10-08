# Code B: response to review (Block B, findings 1 and 3)

Author response for independent re-review. It does not change the reviewer acceptance status:
Block B remains `changes_requested` until the reviewer decides.

**Provenance labels used below**

- **NEW**: computed on 2026-10-08 in the integrated root, from code in this state. This includes
  full runs of scripts 06 and 07 under both time conventions.
- **LEGACY**: saved artifacts of the earlier part-B full run (`Re_Search-Amir-part_b` and the
  team's copy of its `results/forecasts/*.rds`). These were inspected, not rerun. The NEW
  `common_obs` run reproduces them bit for bit (see "Root integration" below), so LEGACY numbers
  equal NEW `common_obs` numbers.

**Inputs not available to the author.** `docs/review_report.md`, `tmp/reviewer/adversarial_review.R`,
`results/logs/reviewer_adversarial.log` and `results/tables/original_vs_reproduced.csv` were not
present in the delivered tree. The reviewer probe was therefore re-implemented from the
instructions in `docs/review_B/b1_probe.R`. It uses a different synthetic panel, so its numbers
differ from 0.2677517 / 0.2232834 / 0.2448805. The mechanism is the same.

## Summary

| Item | Verdict | One line |
|---|---|---|
| B1 | `confirmed` | Future missingness moved h-step targets in the historical pipeline. A trading-calendar convention is now primary, and the historical one is kept as the labelled mode `common_obs`. |
| B2 | `partly_applicable` | Transcription and keys are correct, and the calendar does not explain the gaps. Two causes are identified, filter fragility and one dominant forecast. Period 2 is `unresolved`. |
| B3 | `confirmed` | Smoke MCS was degenerate (block 20 ≥ n ≤ 10, every p-value = 1). It is now reported as unavailable. Full-sample MCS is validated against an independent implementation and against `arch`. |
| B4 | `confirmed` | The handoff overclaimed: cross-class "mfBm wins", "best at every horizon", unqualified "look-ahead free", period 2 "harmless". The narrative is corrected and the numerical code is unchanged. |

## B1. Future missingness and the evaluation calendar — `confirmed`

### Estimand

- **What the paper says.** The paper defines forecasts as "h-day-ahead" (captions of Tables 2–3
  and 15–19, p. 22 and pp. 55–59). Forecast periods are calendar dates, for example
  "January 3, 2007 – January 14, 2025".
- **What it does not say.** The text never mentions missing days in the forecasting section.
  The supplied author code (`author_code/BYZ`) covers estimation and inference only; it has no
  forecasting pipeline. Code A uses the common sample for Table 1 estimation, but that is an
  estimation choice, not a definition of h.
- **What supports common observation steps.** Only numerical agreement: the config comment
  "exactly 500 common days before the first forecast date", and the evidence below.
- **Intended estimand (primary).** The RV of AAPL h trading days after the origin, on the panel
  calendar. The panel calendar is AAPL's own: AAPL is observed on all 5019 DJ30 and 2698 Mag7
  panel dates. The panel follows NYSE trading days; for example, 2007-01-02, the closure for
  President Ford's funeral, is absent.

### Evidence (NEW unless labelled)

**Probe** (`docs/review_B/b1_probe.R`, output `docs/review_B/b1_probe.csv`). Full pipeline
`run_all_models()`, origin row 70, h = 5, on both a daily calendar and an explicit
Monday–Friday calendar. Weekday-calendar results:

| mode | scenario | h = 5 target | fBm forecast | actual |
|---|---|---|---:|---:|
| common_obs | baseline | 2010-04-16 (Fri) | 0.1729426 | 0.1862095 |
| common_obs | NA in asset B after the origin | **2010-04-19 (Mon, 6th trading day)** | 0.1729426 | 0.1665561 |
| common_obs | NA in the target itself at t+5 | **2010-04-19 (moved, not dropped)** | 0.1729426 | 0.1665561 |
| trading | NA in asset B after the origin | 2010-04-16 | 0.1729426 | 0.1862095 |
| trading | future values of B set to 9 | 2010-04-16 | 0.1729426 | 0.1862095 |
| trading | NA in the target itself at t+5 | row removed for every model | — | — |

The forecast never changes, so no future *values* leak into estimates. Only the mapping from
forecast to target depends on future availability. The same holds for mfBm5 and VHAR5
(see the CSV).

**Real samples.** LEGACY forecasts, `docs/review_B/b1_legacy_shifted_targets.csv`.

| sample | h | forecasts | target ≠ h trading days | max extra days |
|---|---:|---:|---:|---:|
| DJ30 | 1 | 4508 | 10 | 1 |
| DJ30 | 5 | 4504 | 48 | 2 |
| DJ30 | 20 | 4489 | 168 (3.7%) | 3 |
| Mag7 | 1 | 2188 | 6 | 1 |
| Mag7 | 20 | 2169 | 120 (5.5%) | 1 |

The NEW trading run has 0 shifted targets (`docs/review_B/b1_trading_shifted_targets.csv`).

**Training windows.** `results/tables/code_b_calendar_keys_{dj30,mag7}.csv` shows that among
origins used by both conventions, 1204 of 4508 (DJ30) and 1144 of 2188 (Mag7) 500-row windows
start on a different day. A complete-case window of 500 rows spans more than 500 trading days
whenever it contains a deleted day.

Missing values in the five assets: DJ30 has 11 days (ALD, AMGN and AXP have 4 each; 9 of these
days fall in 2020–2024). Mag7 has 10 days (GOOG 7, AMZN 2, MSFT 1).

### Is common_obs a replication convention?

Partly, yes. Table 2 at h = 1 (full period and period 1, 5 models each) equals the paper to 4
decimals in **10/10** cells under `common_obs`, against 3/10 under `trading`. At h > 1 neither
convention is exact. The median |relative difference| for the full period plus period 1 is:

| h | trading | common_obs |
|---:|---:|---:|
| 1 | 0.011% | 0.005% |
| 5 | 0.10% | 0.12% |
| 20 | 0.44% | 0.53% |

The honest claim is therefore: "the authors' h = 1 numbers are consistent with windows built on
common observation days". It is not "the paper defines h in common steps".

### Changes

| File | Change |
|---|---|
| `R/rolling_window.R` | `cfg$calendar` with `code_b_calendars()` = `common_obs` (historical, unchanged arithmetic) or `trading`. Adds `code_b_eligible_origins()`. Missing targets are dropped jointly. New output columns `train_start`, `n_train_days`, `calendar`. The checkpoint fingerprint includes the resolved calendar and the coverage rule. |
| `R/forecast_mfbm.R` | Adds `forecast_mfbm_calendar()` and `window_mfbm_params_calendar()` (gap-aware moment estimators). `admissible_window_params()` is factored out with identical arithmetic. |
| `R/forecast_har.R` | Observed-day averages (`na_avg`), `rolling_mean_observed()`, and `har_fit_window()` diagnostics (raw prediction, filter flag, rank). Numerically identical without NA. |
| `R/code_b_config.R` | `calendar = "trading"`, `calendars = c("trading", "common_obs")`, `min_window_coverage = 0.9`. |
| `scripts/06`, `scripts/07` | Compute both conventions. The primary convention keeps the historical file names; `common_obs` files get a suffix. They also write `code_b_calendar_sensitivity_*.csv`, `code_b_calendar_keys_*.csv` and `code_b_provenance_0{6,7}.csv` (input md5, code fingerprint, full config). |

**Rules of the trading convention.** Every decision uses rows 1..t only.

- **Origin t is eligible if** (i) a full window of 500 trading days ends at t; (ii) all five assets
  are observed at t; (iii) at least 90% of the window's days are observed for all five assets.
  Rule (ii) is applied to every model jointly, so all compared models share the same forecast
  keys. Skipped origins: DJ30 10 of 4519, Mag7 6 of 2198, all for reason (ii).
- **Target:** row t + h of the calendar. If the target is missing, the forecast is removed for
  every model and never moved.
- **Gaps inside the training window.** The calendar is neither filled in nor compressed:
  - mfBm: H, σ and ρ come from the moment estimators on genuine 1-day and 2-day increments only
    (increments spanning a missing day are excluded). The conditional mean is exact Gaussian
    conditioning on the observed days at their true times (irregular-time covariance from
    eq. (5)).
  - HAR: the daily term must be observed. The weekly and monthly terms average the observed days
    among the last 5 / 22 trading days. A training pair needs regressors at s and the response
    at s + h.
  - Without missing days both forecasters reduce exactly to the historical code (unit tests).
- **Not adopted: the strict rule "skip if any training day is missing"** (Code C's choice). For
  the fixed DJ30 set it would delete almost all of 2020–2025, because the 9 missing days are less
  than 500 days apart.

### Acceptance tests

All of these go through `run_all_models()`, in `tests/testthat/test-calendar.R`:

- future values and NA patterns of non-target assets leave keys, forecasts and actuals unchanged
  for all 19 models, and targets equal origin + h trading days;
- a missing future target removes that (origin, h) for every model, other horizons are unchanged,
  and no forecast lands on t + 6;
- `common_obs` documents the historical move to t + 6;
- an earlier missing training day keeps the origin eligible: fBm, HAR and LHAR are unchanged
  (60 days used), the multivariate models use 59 days, bfBm equals `forecast_mfbm_calendar()` and
  differs from calendar compression;
- a missing asset at the origin, or <90% window coverage, skips the origin for every model;
- with no missing data both conventions give identical forecasts;
- the gap-aware estimators equal Code A's when nothing is missing;
- the HAR observed-day averages behave as specified;
- the fingerprint differs by calendar and coverage rule.

### Numerical effect of the calendar (NEW, both conventions)

`results/tables/code_b_calendar_sensitivity_{dj30,mag7}.csv` gives every cell under both
conventions. MSFE levels move by at most 0.6% (period 2: 1.3%). Linear-HAR QLIKE moves by up to
11% (Table 18 VHAR4 h = 15: −5.9%). None of the five B2 discrepancies is explained by the
calendar; see B2.

### Root integration

The NEW `common_obs` full forecasts equal the LEGACY forecasts **bit for bit**:

- DJ30: 684,228 rows; Mag7: 331,588 rows; identical keys, forecasts and actuals;
- every table value and MCS p-value is identical (max |diff| = 0).

## B2. Large discrepancies — `partly_applicable` (period 2: `unresolved`)

**Keys and transcription.** `docs/review_B/verify_paper_cells.py` extracts Tables 2, 3 and 15–19
directly from the PDF text layer. All **720/720** cells (loss and MCS p-value) of
`docs/paper_values/code_b_forecast_tables.csv` match the PDF. All five reviewer cells are
confirmed, together with three directly checked reference cells:

| cell | PDF page | paper |
|---|---:|---:|
| Table 17 VHAR5 h = 20 | 56 | 2.5475 |
| Table 19 VHAR3 h = 10 | 59 | 0.1487 |
| Table 18 VHAR4 h = 15 | 58 | 0.2896 |
| Table 3 VHAR4 h = 20 | 22 | 1.9388 |
| Table 2 period 2 fBm h = 1 | 22 | 0.3842 |
| Table 2 full mfBm5 h = 20 | 22 | 1.3370 |
| Table 3 HAR h = 20 | 22 | 1.3072 |
| Table 16 LHAR h = 20 | 55 | 1.2721 |

**Revised comparison (NEW full runs).**

| cell | paper | trading | common_obs | rel. trading | rel. common_obs | n (tr / co) |
|---|---:|---:|---:|---:|---:|---|
| T17 Mag7 VHAR5 h = 20 MSFE×100 | 2.5475 | 1.2826 | 1.2832 | −49.65% | −49.63% | 2173 / 2169 |
| T19 Mag7 VHAR3 h = 10 QLIKE | 0.1487 | 0.3390 | 0.3398 | +127.98% | +128.54% | 2183 / 2179 |
| T18 DJ30 VHAR4 h = 15 QLIKE | 0.2896 | 0.1383 | 0.1470 | −52.24% | −49.25% | 4495 / 4494 |
| T3 DJ30 VHAR4 h = 20 MSFE×100 | 1.9388 | 1.6981 | 1.7016 | −12.41% | −12.24% | 4490 / 4489 |
| T2 DJ30 period 2 fBm h = 1 | 0.3842 | 0.4099 | 0.4100 | +6.68% | +6.72% | 1075 / 1075 |
| T2 DJ30 full fBm h = 1 | 0.5003 | 0.5001 | 0.5003 | −0.03% | 0.00% | 4508 / 4508 |
| T2 DJ30 period 1 mfBm5 h = 1 | 0.6301 | 0.6302 | 0.6301 | +0.01% | +0.01% | 2587 / 2587 |
| T2 DJ30 full mfBm5 h = 20 | 1.3370 | 1.3391 | 1.3416 | +0.16% | +0.34% | 4490 / 4489 |

**Date ranges** (NEW trading, target dates):

- DJ30 full: 2007-01-03 (h = 1) or 2007-01-30 (h = 20) to 2025-01-14;
- period 1: 2007-01-03 to 2017-04-11; period 2: 2017-04-12 to 2021-07-30;
- Mag7: 2016-03-29 (h = 1) or 2016-04-19 (h = 20) to 2025-01-14.

### Investigation, in the requested order

1. **Sample and dates.** Counts and ranges are listed above. Calendar sensitivity is separated
   from everything else: the five cells move by at most 5.9% between conventions, and the gaps
   persist under both.
2. **Units and transformations.** RV = annualised Risk Lab QMLE volatility, and the paper reports
   MSFE×100. The h = 1 equality to 4 decimals (full period and period 1) rules out unit or scale
   errors.
3. **Window endpoints.** Variant: 500 HAR regression pairs instead of 500 − h. It moves the MSFE
   cells by −0.9% and −3.8% and Table 18 by −7.0%, and closes no gap (the gaps become −13.2%,
   −51.6% and −55.6%). Table 19 falls by 74.6% only because the single near-zero forecast
   (item 8) changes, which illustrates its fragility, not a cause; the gap becomes −42.1%
   (`docs/review_B/b2_trading_variants.csv`).
4. **h-step alignment.** The calendar comparison is above. In both conventions the response of a
   training pair satisfies s + h ≤ origin.
5. **HAR / VHAR / VHARF definitions.** Equations (33)–(35): VHARk has 1 + 3k regressors and VHARF
   has 7. VHARF2 = VHAR2 exactly, as in the paper (unit test).
6. **Lognormal correction.** It applies only to LHAR and mfBm, which agree within 2.6% (mfBm
   within 0.44% in Table 18). It is not involved in the linear-VHAR cells.
7. **Rank and conditioning.** Zero rank-deficient fits in all four HAR cells
   (`b2_trading_cells.csv`).
8. **Positivity filter.** This is where the linear-HAR gaps sit:
   - **Mag7 VHAR3 QLIKE h = 10.** One forecast, origin 2020-04-08, has a *positive* OLS
     prediction of 0.000485 against an actual of 0.270. Its QLIKE is 549.5, which is **74% of the
     total loss**. Without the top 10 rows the mean is 0.071. The ≤ 0 filter does not catch
     near-zero positive forecasts, so this level is not robustly estimable. The paper's 0.1487
     could come from a slightly different near-zero forecast; we cannot attribute it further.
   - **DJ30 VHAR4 QLIKE h = 15.** The top 10 rows give 25% of the loss (minimum forecast
     0.0063). The level is likewise dominated by near-zero forecasts, and the cause is not
     identified.
   - **Linear-VHAR MSFE (EXPLORATORY, not adopted).** Hypothesis: the paper's MSFE used raw OLS
     forecasts without our filter. It was tested on **all** 144 linear (V)HAR(F) MSFE cells of
     Tables 3, 15 and 17, not fitted to one cell (`docs/review_B/b2_filter_hypothesis_*.csv`):

     | convention | cells with filter events | within 1% with filter | within 1% without | median abs. error with / without |
     |---|---:|---:|---:|---:|
     | common_obs | 33 | 1 | 19 | 3.5% / 0.8% |
     | trading | 33 | 2 | 12 | 3.6% / 1.1% |

     Examples: Table 3 VHAR4 h = 20 goes from 1.702 to 1.929 (paper 1.939); Mag7 VHAR5 h = 20 goes
     from 1.283 to 2.209 (paper 2.548, still −13%). This is evidence about the paper's convention,
     not proof. The pipeline keeps its documented filter, because QLIKE is undefined for
     RV forecasts ≤ 0.
9. **Loss aggregation and MCS.** Losses are averaged over target dates within the period, one
   row per (target date, model, h). MCS: see B3.

**Table 2 period 2 — `unresolved`.** For fBm, bfBm and mfBm5 at h = 1, our full-period and
period-1 values equal the paper to 4 decimals, yet period 2 does not (0.4100 vs 0.3842). Under
our daily losses, the paper's own full and period-1 numbers imply a 2021-07-31 – 2025-01-14 mean
of 0.236 against our 0.203 (`docs/review_B/b2_table2_consistency_legacy.csv`). So the paper's
period-2 cells are not consistent with the stated period on these forecasts.

- **April 11 vs April 12.** Using the prose date instead of the caption date changes the gap only
  from +6.68% to +6.76%. This confirms that the boundary alone is no explanation.
- **Boundary scan (EXPLORATORY).** Fitting one boundary on fBm h = 1 (end 2020-08-31)
  reproduces 9–10 of 40 period-2 cells within 0.1%, with up to 3.4% error elsewhere. No single
  boundary reproduces the table (`b2_period2_scan_*.csv`).
- The cause is not established. We do not call it a paper error, and we do not use the matching
  ratios as an explanation of the level.

**Scope of "about 4% higher" (Mag7).** It holds for the mfBm class (+2.4% to +5.3%, median
+4.1% in Table 17; +2.1% to +5.3% in Table 19) and for VHARF/VLHAR MSFE (+0.4% to +4.4%). It does
not hold for linear VHAR at long h (down to −49.7%) or for HAR-type QLIKE (−47% to +128%). The
coverage file now states this scope.

**Changes.** No numerical change was made to the B2 pipeline. No periods, assets, filters or
parameters were tuned. All variants are labelled EXPLORATORY and live in `docs/review_B/`. The
documentation no longer claims exact or uniformly close replication: `docs/paper_coverage_code_b.csv`,
`docs/code_b_handoff.md` and `docs/replication_limitations.md` were updated.

## B3. Smoke MCS and independent validation — `confirmed`

**Evidence.** With block ≥ n, `block_bootstrap_index()` draws a single start, and the first n
indices are one full rotation of the sample. Every bootstrap mean then equals the sample mean,
the bootstrap variance is 0, and every p-value is 1 by construction. Smoke uses origin step 500,
which gives n = 2–10 losses per cell against block 20 in all 720 smoke cells. The earlier smoke
p-values were therefore uninformative, not merely noisy.

**Change** (`R/metrics.R`):

- `mcs_status(n, block)` returns `unavailable_block_ge_n` or `unavailable_n_lt_2`, and in those
  cases `mcs_pvalues()` returns NA;
- `evaluate_class()` reports `mcs_status`, `mcs_B`, `mcs_block`, `mcs_seed`, `first_target`
  and `last_target`.

**NEW results:**

- smoke: 720/720 cells are `unavailable_block_ge_n`;
- full: 720/720 are `ok`, with n = 1075–4508, B = 5000, block 20 and seed = 20251004 + h;
- unit test added in `tests/testthat/test-metrics.R`.

The block length was not changed, and there is no reduction for significance.

**Specification kept.** The statistic is Hansen–Lunde–Nason T_max: each model's loss relative to
the average of the models still in the set, with variance from the circular block bootstrap. The
model with the largest t is eliminated, and its p-value is the running maximum of the test
p-values. Note that the paper's Appendix F.4 names T_max but then describes "the maximum absolute
difference between the empirical distribution functions", which is a different statistic. We
implement T_max. The R-type statistic is not used.

**Independent validation (NEW).** `docs/review_B/mcs_export.R` and `mcs_validate.py` were run on
8 full-sample cells, including the four B2 HAR cells:

1. **Exact check.** A separate NumPy implementation of the T_max elimination, fed the identical
   resamples exported from R (B = 1000, seed 7). Max |p_numpy − p_R| = **0** for both the LEGACY
   and the NEW trading losses.
2. **External check.** `arch` 8.0.0 `MCS(method="max", bootstrap="circular", block_size=20,
   reps=5000)`, 5 seeds, using its own resamples. Max |mean p_R − mean p_arch| =
   **0.0067 (LEGACY) / 0.0077 (NEW)**. Our own range across 5 seeds is up to 0.026, so the two
   implementations agree within bootstrap noise. `arch` was installed in a scratch virtualenv
   outside the project.

Bootstrap variability, from 5 seeds of our implementation:

| cell | p range |
|---|---|
| DJ30 Table 2 h = 5 fBm | 0.093 – 0.106 |
| Table 2 h = 20 fBm | 0.023 – 0.028 |
| Table 3 h = 20 VHAR5 | 0.002 – 0.005 |

Full results: `docs/review_B/mcs_{legacy,trading}/mcs_validation.csv`. The exported resample
matrices (`*_idx.bin`, ~140 MB per run) were deleted after the check. `mcs_export.R` regenerates
them deterministically (seed 7).

**Limitation.** All MCS sets are within-class. Smoke p-values are execution checks only.

## B4. Claims — `confirmed`; narrative corrected, no numerical code changed for it

Overclaims found in `docs/code_b_handoff.md`, now corrected:

- **"Why mfBm wins".** This implied cross-class superiority. It is now "Why extra stocks can help
  mfBm (theory)". The handoff states explicitly that at h = 20 HAR and LHAR beat mfBm5, both in
  the paper (1.3072 and 1.2721 vs 1.3370) and in our results (trading: 1.3045 and 1.2705 vs
  1.3391). At h = 10 LHAR and mfBm5 are equal (1.0568 vs 1.0570). No cross-class MCS was run.
- **"mfBm5 is the best model of its class (MCS p = 1) at every horizon".** False even for the
  LEGACY results: at h = 10, p = 0.768 (common_obs) and 0.404 (trading).
- **"The empirical result of the paper is reproduced".** Now qualified: period 2 is unresolved,
  h > 1 differs by up to 0.7%, and h = 1 is exact only under `common_obs`.
- **"Is it look-ahead free? Yes".** Now distinguishes future values (never used) from future
  missingness, which decided the target date in `common_obs`.
- **"Period 2 differs … but the comparison between models is the same".** Now unresolved, with
  the explicit note that ranking agreement does not explain the level gap.
- **Oracle vs estimated.** Props. 4.1–4.3, Figures 4–6 and the known-parameter rows of Tables
  8–11 are now labelled oracle results. The empirical 0.5–3% gains are finite-sample results with
  estimated parameters.
- **Subperiods.** The split chosen from the rolling-H plot (Figure 7) is labelled exploratory.
- **Smaller corrections.** The Mag7 numbers are updated: mfBm5 beats fBm for h ≥ 4 (paper h ≥ 3),
  by 4.6% at h = 20 (paper 5.4%). The 250-day window loses at h ≥ 4, not h ≥ 5. The "about 4%"
  scope is clarified.

The mfBm forecasting identities, which were accepted separately, are unchanged.

## Commands executed (NEW, integrated root, R 4.6.1, macOS, 12 cores)

```sh
Rscript scripts/run_tests.R smoke                      # 74 test blocks, 293 expectations, 0 failures
Rscript scripts/05_run_forecast_simulations.R smoke    # exit 0
Rscript scripts/06_run_empirical_forecasts.R smoke     # exit 0 (both calendars)
Rscript scripts/07_run_robustness.R smoke              # exit 0 (both calendars)
Rscript scripts/06_run_empirical_forecasts.R full      # exit 0; ~11 min per calendar on 11 workers
Rscript scripts/07_run_robustness.R full               # exit 0; Mag7 both calendars + 250-day window
Rscript docs/review_B/b1_probe.R both
Rscript docs/review_B/b1_legacy_real.R <legacy forecasts dir>
python3 -I docs/review_B/verify_paper_cells.py
Rscript docs/review_B/b2_diagnose.R results/forecasts "" trading docs/review_B/b2_trading
Rscript docs/review_B/b2_filter_hypothesis.R results/forecasts "" trading docs/review_B/b2_filter_hypothesis_trading.csv
Rscript docs/review_B/b2_period2_scan.R results/forecasts/code_b_forecasts_dj30.rds docs/review_B/b2_period2_scan_trading.csv
Rscript docs/review_B/mcs_export.R results/forecasts "" docs/review_B/mcs_trading
<venv with arch 8.0.0>/python -I docs/review_B/mcs_validate.py docs/review_B/mcs_trading
```

The same B2 and MCS scripts were also run on the LEGACY forecasts (`*_legacy*` outputs). Logs are
in `docs/review_B/*.log`, and run identity is in `results/tables/code_b_provenance_0{6,7}.csv`.

**Not rerun in full:** script 05. Its simulation code is unchanged, and only its checkpoint
fingerprint changes. The unified `docs/paper_coverage.csv` is regenerated by script 09 (Code C),
which was not rerun; `docs/paper_coverage_code_b.csv` is updated.

**Checkpoints.** Forecast checkpoints are named `<sample>_<calendar>_<fingerprint>_<model>.rds`.
The fingerprint covers the data, the full configuration (including the resolved calendar and the
coverage rule) and the source of every function in the model files. Older checkpoints therefore
cannot be reused for new settings; stale ones remain on disk unused.

**Environment note.** In the delivered tree, the project library `.Rlib` lacked
`brio/libs/brio.so`. Its other compiled packages were R 4.5.2 builds under a macOS quarantine flag,
and macOS refused to load `testthat.so`. Intermediate test runs used the system testthat 3.3.2.
The packages were then rebuilt from source at the `renv.lock` versions: brio 1.1.5, ps 1.9.3,
processx 3.9.0, diffobj 0.3.9, renv 1.3.1, testthat 3.3.2.

Homebrew R 4.6.1 requests `-std=gnu23`, which Apple clang 15 rejects, so the build used a one-off
`R_MAKEVARS_USER` with `CC = clang -std=gnu2x`. The final `Rscript scripts/run_tests.R smoke` ran
from `.Rlib`: 74 test blocks, 293 expectations, 0 failures.

## Remaining limitations

- The trading convention is the stated estimand but does not reproduce the paper at h = 1. The
  paper itself is consistent only with common-observation windows. Both are reported; neither
  reproduces h > 1 exactly.
- Gap handling (irregular times for mfBm, observed-day averages for HAR, 90% coverage) is a
  documented choice. With ≤ 11 missing days per sample its numerical effect is below 0.6% of MSFE.
- The linear-HAR QLIKE levels and the long-horizon linear-VHAR MSFE levels are not robust
  replications. The period-2 gap is unresolved.
- MCS is within-class. B = 5000 bootstrap p-values carry Monte Carlo noise of up to ±0.013.
