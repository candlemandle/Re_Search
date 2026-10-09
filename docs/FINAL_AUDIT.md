# Final audit against the article

Audit date: 2026-10-09. Integration base: `origin/Amir/part_b-review-b` at
`d051338`, with the supplied `code_c/` package added without modifying its scientific
source files.

## Bottom line

The repository is suitable for a transparent replication-plus-extension project.
Code A is a close reproduction. Code B reproduces the main qualitative forecasting
patterns and theoretical figures, but not every published loss level. Code C is a
well-instrumented research extension proposed by Section 6 of the article; it must not
be described as a paper replication.

The defensible claim is:

> We reproduce the paper's estimators, theory and main forecast patterns, document an
> unresolved period-2 level discrepancy, and show exploratory evidence that adding
> stationary mean reversion can improve 20-day AAPL forecasts in 2019–2021.

## Article-to-code correspondence

| Article component | Implementation | Assessment |
|---|---|---|
| mfBm covariance, admissibility, simulation | `R/covariance.R` | aligned with equations (1), (4), (5), (19) |
| Moment estimators and asymptotic SE | `R/estimators.R` | aligned with equations (6)–(10), (14), (27) |
| Time-reversibility test | `R/time_reversibility.R` | aligned with equation (11) |
| Conditional Gaussian forecast | `R/forecast_mfbm.R` | aligned with Section 4 / equation (15) |
| HAR/LHAR/VHAR families | `R/forecast_har.R` | implemented as direct horizon-specific regressions |
| Rolling evaluation and MCS | `R/rolling_window.R`, `R/metrics.R` | no-look-ahead, paired targets, calendar diagnostics |
| Empirical tables/figures | `scripts/03*`, `05`, `06`, `07` | complete coverage recorded in `docs/paper_coverage.csv` |
| mfOU future-work direction | `code_c/R/mfou.R` | genuine stationary mfOU, diagonal drift, reversible mfBm driver |
| Window and asset-selection extensions | `code_c/R/experiments.R`, `selection.R` | exploratory extensions, not published results |

## Verified numerical evidence

- Table 1 point estimates: 25/25 cells match to four decimals.
- Table 1 SE: 23/25 match to four decimals; two differ by 0.0001 because the article's
  finite truncation and the implementation's longer series truncation round differently.
- Monte Carlo: 290/316 cells are within two combined Monte Carlo standard errors.
- A/B test log: 74 blocks and 293 expectations passed on the review branch.
- Code C supplied test log: 23 blocks and 145 expectations passed; its manifest has
  92/92 matching SHA-256 entries.
- Code C source hashes for the baseline A/B files match the integrated branch exactly.
- Code C headline gain percentages were independently recomputed from the saved CSVs;
  the maximum formula discrepancy was approximately `3e-12` percentage points.
- Processed data MD5 values match `data/processed/MD5SUMS`:
  `544fbd178444f5bff3c695c5cd69ead0` (DJ30) and
  `b3852bb24f4d3ea2eca73b1a1566bbab` (Mag7).

The current machine did not have R installed during this integration audit. Therefore
the saved R logs/results were checked, hashes and CSV calculations were independently
verified, and the R source was reviewed statically, but the entire R pipeline was not
re-executed on this machine. `scripts/setup_environment.R` and `scripts/run_all.R` were
added so the final branch has an explicit rerun path.

## Discrepancies that remain

### Code A

The paper and authors' code use opposite signs for the printed `eta` estimator in one
place. The repository exposes both `proof` and `printed` conventions. This does not
change the two-sided test of `eta = 0`.

Measurement noise creates large downward bias in `rho` and `eta`, and weakens the
time-reversibility test. A failure to reject is not evidence that the true asymmetry is
exactly zero.

### Code B

The trading-calendar implementation is the primary specification. The historical
`common_obs` version is retained as sensitivity output because the article does not
fully specify missing-day handling.

- Full-sample and period-1 Table 2 values are close to the paper.
- Period-2 losses remain mostly 6–9% above the paper. Calendar choices, target alignment,
  data units, and filters were checked; no defensible code change eliminated the gap.
- Relative model ratios remain close, so the qualitative result is stronger than the
  claim of exact absolute replication.
- Some VHAR QLIKE values are sensitive to the rule that discards nonpositive forecasts.
  Present them as implementation-dependent robustness results.

### Code C

MfOU changes the model, estimator and initial distribution at once. Its empirical gain
cannot be attributed causally to mean reversion alone. The conditional variance is
plug-in and does not include parameter-estimation uncertainty. `kappa` is constrained
positive and can be weakly identified, especially when true mean reversion is slow.

The dense study covers only dimensions 1–2 and 2019–2021. Window, selection and
five-dimensional comparisons use 12 quarterly origins. The complete daily grid over
all dimensions/windows was not executed.

## Additional-experiment conclusions

### Stationary mfOU versus mfBm

For AAPL and AAPL+ALD, with a 500-day training window and 750 paired origins:

| Assets | Horizon | MSFE gain | QLIKE gain | Interpretation |
|---|---:|---:|---:|---|
| AAPL | 1 | -0.39% | +0.89% | no clear short-horizon gain |
| AAPL | 20 | +12.18% | +9.21% | exploratory intervals exclude zero |
| AAPL+ALD | 1 | -0.53% | +0.88% | no clear short-horizon gain |
| AAPL+ALD | 20 | +12.24% | +9.36% | exploratory intervals exclude zero |

The cumulative-loss plot shows that most of the long-horizon advantage accumulates
during particular episodes rather than uniformly. This is evidence for a horizon- and
regime-dependent benefit, not universal dominance.

### Window length

In Code B, the 250-day window raises five-asset mfBm MSFE by about 2.0% at `h=1`,
3.6% at `h=5`, and 31.9% at `h=20` relative to 500 days. In sparse Code C mfOU
comparisons, both 250 and 1000 days can be worse than 500 days at `h=20`. More history
is not automatically better because the local dynamics and the number of parameters
both matter.

### Asset selection

The proposed score `abs(rho) * abs(H_target - H_asset)` behaves differently by horizon.
For AAPL it loses 3.24% MSFE at `h=1` and gains 5.79% at `h=20` relative to the fixed
article order. Correlation-only gains 2.94% at `h=20`. With only 12 outcomes, this is a
useful hypothesis generator, not a model-selection rule.

### Simulation

At `h=20`, estimated mfOU is worse when the data generator is mfBm and better on average
when the generator is mfOU. With only 30 repetitions, Monte Carlo standard errors are
large. Slow `kappa` is strongly overestimated in the saved study, showing that
identification, not just forecast formula choice, is a central limitation.

## Claims to avoid

- “The whole paper is exactly reproduced.”
- “mfBm beats HAR at every horizon.”
- “Adding more assets always improves forecasts.”
- “The time-reversibility test proves `eta = 0`.”
- “mfOU is statistically superior in general.”
- “The H-gap rule finds the optimal assets.”

## Release checklist

- Run `Rscript scripts/setup_environment.R`.
- Run `Rscript scripts/run_tests.R` and `Rscript code_c/scripts/run_tests.R`.
- Run `Rscript scripts/run_all.R smoke` on a machine with R.
- Inspect `git diff --check` and confirm only intended paths are staged.
- Merge this integration branch into `main`; do not force-push over the one-file current
  main branch.
