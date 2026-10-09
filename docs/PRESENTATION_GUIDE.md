# Presentation asset guide

Use real article figures for the article narrative and real repository outputs for the
team-results section. Do not redraw a table as a decorative chart when a focused crop
of the original table is clearer.

## Main-story figures

| Story point | Asset | What to say |
|---|---|---|
| Roughness changes through time | `results/figures/code_a_fig07_rolling_H_dj30.pdf` | H is not constant; the paper motivates subperiod analysis |
| Why cross-assets can help | article Figure 4/5 or `results/figures/code_b_fig04_weights.pdf` and `code_b_fig05_rel_msfe.pdf` | the extra series changes conditional weights and can reduce MSFE |
| Dimension effect | article Figure 6 or `results/figures/code_b_fig06_dimension.pdf` | gains can grow with dimension in theory, subject to estimation cost |
| Actual AAPL forecast path | `results/figures/code_b_fig_dj30_aapl_forecasts.pdf` | model forecasts track a common target; compare errors, not line aesthetics |
| mfBm gain by horizon | `results/figures/code_b_fig_dj30_msfe_ratio.pdf` | gains are small at h=1 and larger at long horizons |
| Across model classes | `results/figures/code_b_fig_dj30_class_ratio.pdf` | mfBm improves within its class; HAR can still have lower raw loss |
| Mag7 robustness | `results/figures/code_b_fig_mag7_class_ratio.pdf` | short-horizon cost, longer-horizon benefit; vector HAR deteriorates |

## New Code C slides

### 1. Why mfOU?

Use one equation and one verbal contrast:

- mfBm: rough and cross-correlated, but the level is nonstationary.
- mfOU: keeps fractional noise and adds a long-run level plus mean-reversion rate
  `kappa`.

Label this “Extension proposed in Section 6,” not “paper result.”

### 2. Main additional result

Use `code_c/examples/dense_study/figures/cumulative_loss.png` and a compact two-row
table:

| horizon | MSFE gain | QLIKE gain |
|---:|---:|---:|
| 1 | -0.39% | +0.89% |
| 20 | +12.18% | +9.21% |

The chart's y-axis is cumulative `mfBm loss - mfOU loss`; positive values favor mfOU.
State “750 paired AAPL forecasts, 2019–2021.”

### 3. Why the result is not automatic

Use `code_c/examples/dense_study/figures/rolling_parameters.png`. Highlight that H and
`kappa` vary across rolling windows. Add three short caveats: plug-in uncertainty,
weak identification of slow `kappa`, and no multiple-testing correction.

### 4. Robustness / selection

Use a native slide chart from the saved CSVs rather than a screenshot:

- `code_c/examples/window_study/window_sensitivity.csv`
- `code_c/examples/selection_study/metrics.csv`

Show only h=1 and h=20. The message is asymmetric: 500 days is a reasonable compromise;
the H-gap selection rule loses at h=1 and gains at h=20.

## What belongs in backup

- complete parameter tables;
- all MCS p-values;
- period-2 discrepancy diagnostics;
- full Monte Carlo tables;
- covariance derivation and numerical quadrature details.

The main deck should still include real empirical tables/graphs. Backup is for dense
diagnostics, not for the central results.
