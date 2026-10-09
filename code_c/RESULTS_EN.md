# Executed outcomes

Computed on 9 October 2026 using the team's processed panel. **In the daily AAPL study for 2019–2021, mfOU helps more at 20 trading days than at one day.** This is an exploratory result for a specific sample. Full daily execution of the entire paper period, all five assets and all windows has not been performed.

Small tables and figures are in [examples](examples/README.md); full local forecasts remain in `results/code_c_v2/`. Reproduction commands are in the [README](README_EN.md).

## Executed scope

| Study | Actual size | Interpretation |
|---|---|---|
| Daily mfBm/mfOU | 751 scheduled origins; one jointly excluded; 750 pairs per horizon; one/two assets; window 500 | Exploratory paired calendar-block comparison |
| All nested sets/windows | 12 origins, step 63 trading days; nine asset/window scenarios; 1728 forecast rows | Descriptive sensitivity; no inferential intervals |
| Simulation | Three generators × 30 independent replications; 500 training levels; eight horizons | Losses, estimation diagnostics and Monte Carlo error |
| Asset selection | AAPL/MSFT/JPM; 12 origins each; five strategies; 1440 forecast rows | Descriptive heuristic comparisons |
| Smoke | All three stages, including 250/500/1000 windows | Execution verification only |

Daily origins span 2019-01-02 through 2021-12-31; actual target dates span 2019-01-03 through 2022-01-31. The quarterly grid ends at origin 2021-10-07, with last target 2021-11-05. Empirical mfOU studies recorded zero estimation/forecast failures. This does not imply strong identification.

## Daily comparison

Positive gains indicate lower mean mfOU loss relative to mfBm. These are exploratory 95% pointwise intervals from 2000 paired calendar-block resamples, block length 20 trading days. Each row uses 750 pairs; the missing origin retains its original calendar position.

| Assets | Horizon | MSFE gain %, interval | QLIKE gain %, interval |
|---|---:|---:|---:|
| AAPL | 1 | −0.39 [−2.83, 3.37] | 0.89 [−0.72, 2.42] |
| AAPL | 10 | 3.97 [−5.96, 22.49] | 3.60 [−4.77, 13.07] |
| AAPL | 20 | 12.18 [2.16, 34.14] | 9.21 [1.43, 20.99] |
| AAPL + ALD | 1 | −0.53 [−3.31, 3.63] | 0.88 [−0.77, 2.53] |
| AAPL + ALD | 10 | 3.89 [−6.07, 22.27] | 3.60 [−4.20, 13.16] |
| AAPL + ALD | 20 | 12.24 [2.44, 34.85] | 9.36 [1.25, 21.48] |

For two assets and horizon 20, MSFE falls from **0.01804745** to **0.01583913**, and QLIKE from **0.11034045** to **0.10000767**. All absolute losses and horizons are in `examples/dense_study/metrics.csv`.

At block lengths 40 and 60, the 20-day gain intervals still have positive lower bounds. For two assets and block 60: MSFE [3.62%, 31.11%], QLIKE [1.38%, 20.53%]. This demonstrates sensitivity to the tested blocks, not control for searching across models/horizons/periods. The one- and ten-day intervals include zero. Adding ALD changes these mean outcomes only slightly.

## Training-window sensitivity

On twelve matched origins, five-asset mfOU with window 250 has 9.14% higher 20-day MSFE and 8.10% higher QLIKE than window 500. Window 1000 is also worse: MSFE rises 16.69%, QLIKE 20.13%. At one day, window 1000 lowers MSFE by 1.98% but raises QLIKE by 1.48%.

More history therefore does not automatically improve this study. Twelve dates do not prove that 500 is universally best, and no window is selected using future losses. Mean AAPL κ estimates in the five-asset cases are approximately 6.68/4.49/3.21 per year for 250/500/1000 rows. This reveals window sensitivity rather than a constant identified true rate.

With five assets and the main 500-row window, quarterly 20-day gains are 13.82% MSFE and 13.73% QLIKE. These descriptive numbers do not inherit the daily two-asset intervals.

## Controlled simulations

Generators: mfBm; stationary mfOU with κ=(0.2,0.4); stationary mfOU with κ=(2,4). All have H=(0.2,0.4), σ=(1,1), driver correlation 0.4. True parameters are available only here.

| Generator | 20-day QLIKE improvement: estimated mfOU vs estimated mfBm | Paired loss difference ± MCSE |
|---|---:|---:|
| mfBm | −5.74% | +0.00724 ± 0.00911 |
| Slow mfOU | +11.42% | −0.01523 ± 0.02130 |
| Fast mfOU | +10.83% | −0.01363 ± 0.00881 |

The difference here is `loss_mfOU-loss_mfBm`; negative favors mfOU. MCSE is the standard error across independent replications, not a 95% interval. Thirty replications leave substantial uncertainty and do not establish a stable model ranking.

A useful extra outcome is weak-rate estimation difficulty. For slow mfOU's first component, true κ=0.2, estimated bias is approximately +2.81 and RMSE 5.53. For its second component, true κ=0.4, bias is +1.02 and RMSE 1.99. Starting-value disagreement and bound solutions are recorded. The run has 2640 successful forecast rows and 240 `not_applicable` rows, with zero registered failures.

`mfOU_known` is not applicable under an mfBm generator: no positive true stationary κ exists. `mfBm_known_driver` under an mfOU generator knows noise parameters but uses a misspecified covariance; it is not the correct oracle. A small-simulation estimated-parameter win over a correct known-parameter forecast does not mean that estimation improves the population-optimal forecast; inspect MCSE.

## Asset selection

The H-gap heuristic is not a universal winner. For AAPL at horizon one it worsens MSFE by 3.24% and QLIKE by 2.99% relative to the fixed set; at horizon 20 it reduces them by 5.79% and 5.23%. Each comparison has only twelve outcomes and is descriptive. MSFT/JPM, HAR/fBm results and losses are also preserved.

## Verification evidence

- C: **23 blocks / 145 expectations**, no errors or warnings; repeated successfully with merged root layout.
- Unmodified A/B: **66 blocks / 275 expectations** in a disposable copy.
- The complete C smoke and a root-layout execution both exit 0.
- Direct B parity on real 500-day one/two/five-asset windows: maximum mfBm forecast difference 0. A separate baseline audit includes real gaps, exact keys/actuals and future non-target missingness invariance.
- Refining mfOU quadrature from 128 to 256 changes real-window forecasts by approximately 10⁻¹⁴ relative to their scale, below the prespecified 10⁻⁴ tolerance.
- The mfOU covariance is invariant to translating the time origin. B's local window anchor prevents absolute calendar labels becoming a baseline parameter.
- All 281 initially fingerprinted A/B snapshot files retain their SHA-256 identities.

The daily study initially used per-set eligibility. When B's common calendar rule was identified, the original RDS was preserved and a past-only shared mask removed sixteen one-asset forecast rows from one origin without refitting retained forecasts. Original forecasting provenance and current analysis identities remain separate. Final smoke/tests use the corrected common-calendar runner.

## Research hypotheses and limits

- **H1:** positive exploratory evidence at horizon 20; no universal improvement established.
- **H2:** the observed long-horizon comparison is more favorable than the short-horizon comparison. A dedicated test of the difference between those gains was not performed.
- **H3:** stationarity is a model property; better fit to true long-run dynamics is not established. Window sensitivity and weak κ identification are demonstrated.
- **H4:** estimation adds uncertainty, with both gains and losses across simulated generators; thirty replications are not conclusive.

A fitted positive κ under a positive-κ constraint is not a statistical proof of mean reversion. Cumulative losses reveal contributions of changing volatility, but a separate causal crisis analysis was not performed. Trading/risk-management economic utility was not measured. Model, estimator and initial distribution change together, so the entire gain cannot be attributed to mean reversion alone.

## Five statements for one presentation slide

1. **Limitation:** mfBm captures roughness and asset dependence but lacks a stationary mean-reverting level.
2. **Modification:** genuine stationary mfOU with component-specific rates, retaining B's baseline.
3. **Experiment:** identical past windows, assets, trading calendar and outcomes; extra window/simulation studies.
4. **Result:** AAPL, 2019–2021, 750 pairs: approximately 12% lower MSFE and 9% lower QLIKE at 20 days; no convincing one-day gain.
5. **Conclusion:** promising for this sample and horizon, with weak identification and exploratory uncertainty requiring further verification.
