# PR: Code A — data, mfBm covariance structure, estimators, test, Monte Carlo

> Updated after the first review: see `docs/response_to_review_A.md`.

Reviewer: participant 5. Branch: `code-a`.

## Команды запуска

From the project root, in a clean R session:

```bash
Rscript scripts/00_check_author_code.R          # authors' files vs ours
Rscript scripts/01_download_data.R              # Risk Lab, ~20 MB
Rscript scripts/02_prepare_data.R
Rscript scripts/03_estimate_parameters.R full
Rscript scripts/04_run_monte_carlo.R smoke      # ~1.5 min, same code as full; writes to results/smoke/
Rscript scripts/04_run_monte_carlo.R full       # ~15-25 min on 8 cores
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Single experiment: `Rscript scripts/04_run_monte_carlo.R full t7` (`e1`, `t6`, `t7`, `t12`, `g`).
Workers: `MFBM_CORES` (default `detectCores() - 1`); results do not depend on it.

## Использованные данные

* Risk Lab (Dacheng Xiu), daily volatility, column `qmle_trade` (annualized QMLE from all
  trades), `B = log(vol)`. Endpoint, file format and ticker → PERMNO mapping: `data/README.md`.
* DJ30: 30 paper tickers, 2005-01-05 – 2025-01-14 (5019 dates). Mag7 subset: AAPL, AMZN, FB,
  GOOG, MSFT, 2014-03-27 – 2025-01-14.
* Authors' R files `author_code/BYZ/` from <https://fba.um.edu.mo/wp-content/uploads/2025/04/BYZ.zip>
  (link in arXiv v1), unchanged.

## Созданные outputs

Tables (`results/tables/`):

| file | paper |
|---|---|
| `code_a_data_summary.csv` | data/period summary |
| `code_a_table01_estimates.csv` | Table 1 (+ paper values and SEs) |
| `code_a_table13_14_dj30.csv` | Tables 13–14 (+ paper values, SEs) |
| `code_a_time_reversibility_{dj30,mag7,summary}.csv` | test (11) on all pairs |
| `code_a_mag7_estimates.csv` | Mag7 estimates (supplementary) |
| `code_a_mc_table04_05_estimators.csv`, `..._vs_paper.csv` | Tables 4–5 |
| `code_a_mc_table06_byz_vs_ac.csv`, `..._vs_paper.csv` | Table 6 |
| `code_a_mc_table07_test_size_power.csv`, `..._vs_paper.csv` | Table 7 |
| `code_a_mc_table12_measurement_error.csv`, `..._vs_paper.csv` | Table 12 |
| `code_a_mc_table20_gmm.csv`, `..._vs_paper.csv`, `code_a_mc_gmm_optimal_weights_true.csv` | Table 20 |
| `code_a_author_code_comparison.csv` | authors' code vs ours |
| `code_a_original_vs_ours_empirical.csv` | short original-vs-ours (empirical) |
| `code_a_fig07_rolling_H_dj30.csv`, `code_a_fig08_rolling_H_mag7.csv` | data of Figures 7–8 |

Figures (`results/figures/`): `code_a_fig01_sample_paths.pdf`, `code_a_fig02_rho_max.pdf`,
`code_a_fig03_admissible.pdf`, `code_a_fig07_rolling_H_dj30.pdf`, `code_a_fig08_rolling_H_mag7.pdf`,
`code_a_mc_power_curve.pdf`, `code_a_fig09_gmm_estimates.pdf`, `code_a_fig10_gmm_variances.pdf`.

Log: `results/logs/code_a.log`. Coverage: `docs/paper_coverage.csv` (19 Code A rows).

## Runtime

8-core MacBook (machine shared with other apps): download ~2 min; 02 ~2 s; 03 ~15 s;
04 smoke ~1–2 min; 04 full: Tables 4–6, 12 < 10 s each, Table 7 ~1 min, Appendix G
(1000 reps × two-step GMM + 1000-draw bootstrap) dominates, ~5–25 min depending on load.
Tests ~10 s. Dependencies: R ≥ 4.3, `testthat`; everything else is base R.

## Что совпало с оригиналом

* **Table 1:** all 25 point estimates equal the paper at 4 decimals; 23/25 SEs do, two differ by 1
  in the 4th decimal (authors' finite-sum SE formula reproduces 25/25; `scripts/03b_check_table01_se.R`).
* **Asymptotic SEs of Tables 4–5 (parentheses):** 45/48 equal after rounding; the closed forms
  (14), (27) agree with an independent Isserlis computation and with the authors'
  `BYZinference.R` to 1e-5 (relative 1e-4).
* **Author code:** our MM estimators equal `BYZestimator.R` on the same path (H, sigma^2, rho to
  machine precision); our AC estimator equals their `estMFBM(i2, M = 1..5, w = (1, 0, 0))`.
* **Monte Carlo (Tables 4–7, 12, 20):** 290 of 316 cells within 2 combined MC standard errors
  of the paper (criterion and exceptions: `results/tables/code_a_mc_agreement*.csv`,
  `docs/response_to_review_A.md`); e.g. Table 12 noise bias rho −0.2436 vs −0.2430;
  Table 6 AC eta bias at eta = 0: 0.304 vs 0.308; Table 20 n·Var(rho_opt) 0.66 vs 0.65.
* Figures 1–3, 7–8 reproduce the paper's figures (Figure 7 incl. the 2017–2021 convergence).
* rho_max(0.2, 0.8) = 0.662 and rho_max(0.1, 0.9) = 0.383 as quoted in Section 2.

## Что отличается и почему

1. **Tables 13–14 (30 stocks):** 70% of H, 57% of rho and 55% of eta equal to 4 decimals,
   max |diff| 0.0105 (H JNJ) / 0.0038 (rho) / 0.0091 (eta), all < 1 SE. Cause: sample rules for
   missing days. The paper copies Table 1's H for the first five stocks and seems to use the
   five-stock common sample for pairs involving them (this hybrid rule gives 83% exact matches);
   we keep a single transparent rule (pairwise common dates), documented in `data/README.md`.
2. **Sign of eta.** Printed (8b) and the authors' code estimate eta in the Amblard–Coeurjolly
   convention (simulator uses `rho - eta sign(h)`), opposite to the paper's definition (4)/(26)
   and the Appendix C.2 proof. We default to the eq. (4) convention (`eta_mm(..., "proof")`);
   `"printed"` reproduces (8b). Table 14 entry (row i, col j) = our eta(col j, row i).
   No effect on the test (two-sided) or any magnitude.
3. **Table 7, H_hat < 0.** With H1 = 0.1, n = 500 about 1.2% of replications give H_hat1 < 0,
   where AVAR is undefined. Main columns: rates over valid replications; `*_all` columns
   count them as non-rejections. The paper does not say which it uses.
4. **Table 7, eta = 0.65** is 0.04 from the admissibility boundary (0.693, Amblard et al. 2013,
   Prop. 9; verified with the exact covariance matrix). Exact circulant embedding fails there,
   so the standard approximate Wood–Chan method is used (relative clipped eigenvalue 1.6e-5),
   flagged `approximate_simulation = TRUE`.
5. **Table 20, n·Var of rho_1 / rho^(2)** from the 1000 GMM replications are 0.85 / 1.31
   (paper 0.77 / 1.16, asymptotic 0.777 for rho_1); with 10 000 independent replications
   (`n_variance_large_mc`) they are 0.78 / 1.20, i.e. the paper's values. The 1000-replication variance has ~5% MC error.
6. **sigma^2 SE (10)** uses log(1/Delta), as in the proof (the theorem writes log n with Delta = 1/n);
   this is what reproduces Tables 4–5.
7. Our simulator: exact multivariate circulant embedding (validated against Cholesky in tests)
   instead of the authors' `simMFBM`; different random numbers, same distribution.

## Какие тесты проходят

`tests/testthat/test-estimators.R`: 30 test blocks, all pass (3 added after review). They cover
covariance symmetry/PD, consistency of level and increment covariances, rho_max values,
admissibility boundary, simulator vs target covariance and vs Cholesky, near-boundary fallback,
seed reproducibility, estimator recovery and admissible ranges, scale invariance, eta sign
convention, asymptotic SEs vs Tables 4–5, closed forms vs Isserlis, GMM weights, test statistic,
invalid-H handling, Risk Lab parser and cleaning, panel alignment, smoke/full config structure,
core-independence of Monte Carlo, absence of absolute paths, and Table 1 on real data
(skipped if the data are not downloaded).

## For the presenter

`docs/code_a_handoff.md`: method in 7 sentences, two main results, the main limitation,
legend of every figure.
