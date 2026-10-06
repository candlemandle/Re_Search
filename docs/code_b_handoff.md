# Code B: handoff for the presenter (via participants 4–5)

Only results accepted by the reviewer should go on slides. All numbers below come from
`results/tables/code_b_*.csv`.

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

## Two main results

1. **The empirical result of the paper is reproduced (Table 2).** For AAPL with ALD, AMGN, AXP, BA
   the MSFEs at h = 1 match the paper to 4 decimals, and the gain of mfBm5 over fBm matches to
   0.01 percentage points: 0.5% at h = 1, 1.3% at h = 5, 2.3% at h = 20 over 2007–2025, larger in
   2007–2017 (1.8% at h = 5, 3.2% at h = 20) when the H of the stocks differ, and none in
   2017–2021 when they converge (Figure 7 of Code A). mfBm5 is the best model of its class (MCS
   p = 1) at every horizon.
2. **More information helps mfBm but hurts HAR (Tables 3, 15–17).** Every vector HAR variant is
   worse than the plain HAR, and worse with every added stock (VHAR5 is 8% worse at h = 1 and 38%
   at h = 20), while mfBm gets better. Same on the Mag7 stocks (mfBm5 −4.8% vs fBm at h = 20).

## Why mfBm wins (two sentences)

The optimal forecast of an mfBm uses the other stocks only through the *difference in roughness*:
when H differ and the correlation is non-zero, the other path tells something about AAPL's
future that its own past does not (Prop. 4.1; with equal H the extra weight is exactly zero).
HAR-type models estimate 3 extra coefficients per added stock on 500 days and pay more in
estimation noise than they gain, while mfBm adds only one H, one sigma and the correlations.

## Main limitation

The gains are small (0.5–3% of MSFE) and depend on the period: in 2017–2021, when the H of the five
stocks are close, mfBm and fBm are the same. With a one-year window (250 days) mfBm5 loses to fBm
at h ≥ 5, so the gain needs enough data to estimate the 5 × 5 correlation matrix.

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
  rolling window; both give the paper's Table 2 numbers exactly at h = 1.
* **Is it look-ahead free?** Yes: every forecast uses only data up to its origin; this is checked
  automatically on every run and in a unit test that changes all future data.
* **Lognormal correction.** mfBm and log-HAR forecast log RV; RV forecast = exp(mean + var / 2).
* **Period 2 numbers differ** from the paper by up to 9% in level, but the comparison between
  models is the same as in the paper.

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
