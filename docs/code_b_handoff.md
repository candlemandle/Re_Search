# Code B: handoff for the presenter (via participants 4–5)

Only results accepted by the reviewer should go on slides. All numbers below come from
`results/tables/code_b_*.csv` of the full run of 2026-10-08 (trading convention unless stated).
Reviewer status: B block pending independent re-review.

## The metrics

* **MSFE** (mean squared forecast error) = average of (realized RV − forecast)² over all forecast
  days. Smaller is better. The paper prints MSFE × 100. RV here is annualized volatility (≈ 0.25),
  so MSFE × 100 = 0.50 means a typical error of about √0.005 ≈ 0.07, i.e. 7 vol points.
* **RMSFE** = √MSFE, used in the simulations (Tables 8–11), same units as the forecast.
* **QLIKE** = RV/F − log(RV/F) − 1. Zero for a perfect forecast; penalizes forecasts that are too
  low more than forecasts that are too high. Standard robust loss for volatility (Patton 2011).
* **MCS p-value** (model confidence set, Hansen et al. 2011): the set of models that cannot be
  distinguished from the best one. p = 1 means "the best model"; a model with p < 0.10 is
  statistically worse than the best at the 10% level. Computed within each model class and horizon.
* **Ratio to the univariate model** (our summary plots): MSFE of the 5-asset model / MSFE of the
  1-asset model of the same class. Below 1: the extra assets help.

## The models, from the simplest to the most complex

1. **HAR**: regression of RV in h days on today's RV, the last week's mean and the last month's
   mean (Corsi 2009). The standard benchmark.
2. **LHAR**: the same in log RV.
3. **fBm**: log RV of AAPL is a fractional Brownian motion; forecast = conditional mean given the
   last 500 days (H and sigma estimated on these 500 days).
4. **VHAR2–5**: HAR with the daily/weekly/monthly RV of 1–4 other stocks as extra regressors
   (3 more coefficients per stock).
5. **VHARF2–5**: HAR plus the HAR terms of the *average* of the other stocks (only 3 extra
   coefficients in total).
6. **VLHAR2–5**: VHAR in log RV.
7. **bfBm, mfBm3–5**: log RV of AAPL and 1–4 other stocks are a multivariate fBm (each with its own
   H, joint correlations rho); forecast of AAPL = conditional mean given the past of all of them.
   Only d H's, d sigma's and d(d−1)/2 correlations are estimated.

## Time convention (review B1)

Every empirical number exists in two versions (see `docs/response_to_review_B.md`):

* **trading** (primary, files without suffix): h = trading days of the panel calendar; the
  target of an origin is fixed in advance and a missing target removes the forecast for all
  models. This is the "h-day-ahead" forecast the paper describes.
* **common_obs** (historical replication mode, files `*_common_obs.csv`): days on which one of
  the five stocks is missing are deleted first, so a future missing day of another stock moves
  the target by a day. This is the convention of the earlier results and the only one that gives
  the paper's Table 2 at h = 1 to 4 decimals (10 of 10 cells vs 3 of 10 for trading).

Between the two, MSFE levels move by at most 0.6% (1.3% in period 2); linear-HAR QLIKE moves by
up to 11% because it is dominated by a few near-zero forecasts.

## Two main results

1. **Table 2 is closely reproduced for the full period and period 1, not for period 2.** With the
   historical common_obs convention the h = 1 MSFEs equal the paper to 4 decimals; at h > 1 the
   levels differ by up to 0.7% under either convention. The MSFE gain of mfBm5 over fBm is
   0.48% / 1.23% / 2.24% at h = 1 / 5 / 20 over 2007–2025 (paper 0.48 / 1.27 / 2.30) and
   0.75% / 1.80% / 3.21% in 2007–2017 (paper 0.76 / 1.79 / 3.20). In 2017–2021 there is no gain
   (−0.2% to 0.0%). Period 2 levels are mostly 6–9% above the paper (range −1.4% to +9.1%) and the
   cause is unresolved (see below). Within the mfBm class mfBm5 has MCS p = 1 at every horizon of the full period except
   h = 10 (p = 0.40 trading, 0.77 common_obs; the earlier handoff said "every horizon", which was
   not true even for the historical results).
2. **Within each class, more stocks help mfBm and hurt the vector HAR models (Tables 3, 15–17).**
   Every vector HAR variant has a larger MSFE than the plain HAR (VHAR5 is 8% worse at h = 1 and
   38% worse at h = 20). On the Mag7 stocks mfBm5 beats fBm for h >= 4 (paper: h >= 3), by 4.6% at
   h = 20 (paper 5.4%), and is worse at h = 1 (−1.7%; paper −1.0%). These are within-class
   comparisons; the MCS is run separately in every class and says nothing about which class is
   best. Mag7 levels are 2–5% above the paper for the mfBm class and for VHARF/VLHAR MSFE; this
   "about 4%" statement does not hold for linear VHAR at long horizons or for HAR-type QLIKE.

## What this does and does not show

* **Not "mfBm beats everything".** Across classes mfBm is not uniformly best, already in the
  paper: at h = 20 the paper has MSFE × 100 of 1.3370 for mfBm5, 1.3072 for HAR and 1.2721 for
  LHAR (ours, trading: 1.3391, 1.3045, 1.2705). At h = 10 LHAR and mfBm5 are equal (1.0568 vs
  1.0570). No cross-class MCS was run, so no cross-class superiority is claimed.
* **Theory vs estimated forecasts.** Prop. 4.1–4.3, Figures 4–6 and the "known parameter" rows of
  Tables 8–11 are oracle results: the multivariate forecast cannot be worse than the univariate one
  when the true parameters are used. The empirical gains above are finite-sample results with
  estimated H, sigma and rho on 500 days, and they are small (0.5–3%).
* **Subperiods are exploratory.** The 2007–2017 / 2017–2021 split was chosen after looking at the
  rolling H estimates (Figure 7), so the larger gain in period 1 is suggestive, not a pre-registered
  test.

## Why extra stocks can help mfBm (theory, two sentences)

With known parameters the optimal mfBm forecast uses the other stocks only through the
*difference in roughness*: when H differ and the correlation is non-zero, the other path tells
something about AAPL's future that its own past does not (Prop. 4.1; with equal H the extra weight
is exactly zero). Each added stock costs a vector HAR 3 extra coefficients on 500 days and only one
H, one sigma and the correlations for mfBm, which is consistent with, but does not prove, the
within-class pattern above.

## Main limitation

The gains are small (0.5–3% of MSFE) and depend on the period: in 2017–2021, when the H of the five
stocks are close, mfBm and fBm are the same. With a one-year window (250 days) mfBm5 loses to fBm
at h ≥ 4 (by 19% at h = 20), so the gain needs enough data to estimate the 5 × 5 correlation matrix.

## Section 4 and the simulations

* Figure 4: with one observation of each stock, the weight on the other stock is negative when
  H2 < H1, zero when H2 = H1 and positive when H2 > H1.
* Figure 5: the bivariate forecast always has MSFE ≤ the univariate one; where an mfBm exists the
  gain is up to 5% with rho = 0.5 and up to 19% with rho = 0.9 (H2 = 0.71).
* Figure 6: adding equally correlated components helps less and less; the MSFE ratio settles at
  0.891 (t = 1) and 0.774 (t = 10), as quoted in the text.
* Tables 8–11 (simulated mfBm, parameters known): the theoretical RMSFEs equal the paper exactly;
  the bivariate forecast helps only when rho ≠ 0 and H differ (−12% RMSFE at rho = 0.8), each
  extra correlated component helps more (Table 10), and with estimated parameters the gain
  remains (Table 11).

## Things worth saying if asked

* **Which RV and which window?** Risk Lab QMLE volatility (identified by Code A) and a 500-day
  rolling window; with the common_obs convention they give the paper's Table 2 numbers exactly
  at h = 1 (full period and period 1).
* **Is it look-ahead free?** No future *values* enter any estimate or forecast (unit tests change
  all future data). In the historical common_obs mode, however, future *missingness* of another
  stock decides which day is the h-step target (review B1); the primary trading convention fixes
  the target at the origin and is tested against future values and missing-data patterns.
* **Lognormal correction.** mfBm and log-HAR forecast log RV; RV forecast = exp(mean + var / 2).
* **Period 2 numbers differ** from the paper by up to 9% in level (median +6.6%). Our full-period and period-1
  numbers equal the paper at h = 1, so the paper's period-2 cells are not consistent with the
  stated period on our forecasts; no single alternative boundary reproduces them. That the model
  ratios agree does not explain the level difference; the cause is unresolved.
* **Large HAR differences.** Mag7 VHAR3 QLIKE at h = 10 is dominated by one near-zero forecast
  (74% of the loss); linear VHAR MSFE at long horizons depends on the positivity filter. These
  cells are not robust replications; see `docs/response_to_review_B.md`.

## Figures (axes and lines)

| file | what it shows |
|---|---|
| `code_b_fig04_weights.pdf` | x = H2, H1 = 0.4, rho = 0.5, t = h = 1. Left/middle: relative weights on own / other stock; right: own weight relative to the univariate weight. Red dotted = no mfBm exists |
| `code_b_fig05_rel_msfe.pdf` | x = H2; y = MSFE of the bivariate forecast / univariate MSFE; left rho = 0.5, right rho = 0.9 |
| `code_b_fig06_dimension.pdf` | x = number of components d; y = MSFE ratio to d = 1; dashed = limit d → ∞; left t = 1, right t = 10 |
| `code_b_fig_dj30_msfe_ratio.pdf` | x = horizon h; y = MSFE of bfBm/mfBm3/4/5 divided by fBm; panels: 2007–2025, 2007–2017, 2017–2021; below the grey line = better than fBm |
| `code_b_fig_dj30_class_ratio.pdf` | x = h; y = 5-stock model / 1-stock model for each class; solid = ours, dotted = paper; black (mfBm) below 1, the HAR classes above |
| `code_b_fig_mag7_class_ratio.pdf` | the same for the Mag7 stocks |
| `code_b_fig_dj30_aapl_forecasts.pdf` | realized RV of AAPL (grey) and the fBm (black), mfBm5 (red), HAR (blue dotted) forecasts for h = 1 (top) and h = 20 (bottom), log scale; dashed lines = the two subperiods |
