# What Code C does, in plain language

We check whether allowing volatility to return toward a usual level helps predict it.

The article's mfBm model captures rough fluctuations and relationships between assets. Its levels do not have a stationary long-run distribution. The authors explicitly propose mfOU forecasting as future work. MfOU adds a long-run level and a speed of return, while keeping fractional noise.

More parameters may help, but estimating them from limited history can also make forecasts worse. We therefore compare actual errors rather than assume an improvement.

1. Predict AAPL with the reviewed mfBm code and with stationary mfOU, using the same history, assets and future outcomes.
2. Add ALD, AMGN, AXP and BA in the article's order; compare 250, 500 and 1000 days of history.
3. Generate artificial data with known parameters, then compare known-parameter and estimated-parameter forecasts. This reveals estimation effects.
4. Retain the earlier question of how to select additional assets: fixed sets, correlation, and correlation combined with Hurst differences. The last rule is a heuristic.

The program saves forecasts, observed values, errors, parameter estimates and failures. MSFE, RMSFE and QLIKE are error measures: lower is better. Positive percentage gains favor mfOU. A cumulative plot shows when gains or losses appeared.

C lives in its own folder and directly loads A/B's existing functions. Missing days keep their original positions; no future value enters fitting or asset choice. This makes the mathematical extension compatible with the team's project without editing their files.

The recorded daily 2019–2021 study has a promising 20-day result, around 12% lower squared error and 9% lower QLIKE. One-day performance is mixed. This is an exploratory result for one target and one period, not proof that mfOU always wins. Simulations also show that weak mean reversion can be difficult to estimate.

For a presentation: “Code C tests a future-work direction explicitly proposed by the article. It compares the reviewed mfBm forecast with stationary mfOU on identical trading-calendar outcomes, then studies window length and parameter-estimation effects. The usefulness of mean reversion depends on observed losses and uncertainty.”

Detailed usage: [README](README_EN.md). Actual outcomes: [results](RESULTS_EN.md). Formula definitions: [methods](docs/METHODS_EN.md).
