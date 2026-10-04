# Code A — handoff for the presenter (via participants 4–5)

Only results accepted by the reviewer should go on slides. All numbers below come from
`results/tables/code_a_*.csv`.

## The method in 7 sentences

1. Each log-RV series is modelled as a fractional Brownian motion; the Hurst exponent H
   measures roughness (H < 1/2: rough, anti-persistent increments).
2. A multivariate fBm (mfBm) lets every asset keep its own H and scale sigma, and links
   assets through a correlation rho and an asymmetry (lead–lag) parameter eta.
3. Not every combination is possible: the more two Hurst exponents differ, the smaller the
   maximal correlation rho_max(H1, H2) (Figures 2–3).
4. The paper estimates everything by simple moments of increments: H from the ratio of
   lag-2 to lag-1 realized variances (6), sigma^2 from the realized variance (7), rho as the
   realized correlation of increments (8a), eta from lagged cross-products (8b).
5. The estimators are consistent and asymptotically normal for H < 3/4, with closed-form
   standard errors (Theorems 3.1–3.3); so we can test time-reversibility, eta = 0 (11).
6. Monte Carlo shows small bias and that the asymptotic SEs match the simulated ones; the new
   eta estimator is far more accurate than the older Amblard–Coeurjolly estimator.
7. On 20 years of DJ30 data the estimated H are 0.14–0.25, correlations of increments
   0.2–0.77, and eta is not significantly different from 0 for most pairs, which justifies
   the time-reversible mfBm used for forecasting (Code B).

## Two main results

1. **Exact empirical replication of Table 1.** With Risk Lab `qmle_trade` volatility and the
   five stocks' common sample, all 25 estimates and all 25 standard errors match the paper
   to the 4th decimal (e.g. H_AAPL = 0.2609 (0.0123), rho_AAPL,ALD = 0.3902 (0.0131),
   eta_AAPL,ALD = 0.0950 (0.0322)). Tables 13–14 (30 stocks) match within < 1 SE.
2. **The simulation evidence of the paper is reproduced.** Asymptotic SEs (Tables 4–5,
   parentheses) are reproduced exactly; MC bias/SD/RMSE, the BYZ-vs-AC comparison (Table 6),
   the size and power of the test (Table 7) and the noise bias (Table 12: rho biased by
   −0.24, eta by −0.37) agree within Monte Carlo error.

Time-reversibility on real data: 27 of 435 DJ30 pairs are rejected at 1% (6.2%), 66 at 5%;
1 of 10 Mag7 pairs at 1%. Most pairs are consistent with eta = 0, as the paper states.

## Main limitation

The empirical estimates describe **realized** volatility, not latent volatility: measurement
noise biases rho and eta towards zero (Table 12, rho 0.4 → 0.16) and H downwards. So
"eta ≈ 0" in the data does not prove the latent process is time-reversible (paper, footnote 7).

## Things worth saying if asked

* **Sign of eta.** The printed eq. (8b) and the authors' code estimate eta in the
  Amblard–Coeurjolly convention; the paper's own definition (eq. (4), Appendix C.2) has the
  opposite sign. Same numbers, opposite labels. The test is unaffected.
* **Which RV?** The paper does not say which Risk Lab series it uses; we identified it
  (`qmle_trade`, the noise-robust QMLE estimator) by matching Table 1 exactly. 5-minute RV
  gives clearly lower H and rho.
* **Author code.** The authors' three R files (`author_code/BYZ`) run; their estimators and
  standard errors coincide with ours to machine precision / 1e-5.

## Figures (axes and lines)

| file | what it shows |
|---|---|
| `code_a_fig01_sample_paths.pdf` | two correlated fBm paths (rho = 0.8): top H = 0.1 (rough), bottom H = 0.4 (smoother); x = time in years, step 1/250 |
| `code_a_fig02_rho_max.pdf` | grey level = maximal admissible correlation for each (H1, H2); red contours at 0.25/0.5/0.75/0.9; white diagonal = no restriction when H1 = H2 |
| `code_a_fig03_admissible.pdf` | grey = (H1, H2) pairs where a time-reversible bfBm with rho = 0.5 / 0.75 / 0.9 exists |
| `code_a_fig07_rolling_H_dj30.pdf` | H estimated on rolling two-year windows of log RV for AAPL, ALD, AMGN, AXP, BA; dashed lines at 2017-04-11 and 2021-07-30 (forecast subperiods of Code B) |
| `code_a_fig08_rolling_H_mag7.pdf` | same for AAPL, AMZN, FB, GOOG, MSFT |
| `code_a_mc_power_curve.pdf` | rejection rate of the test vs true eta; solid 5%, dashed 1%; blue n = 500, red n = 1000; grey lines at nominal sizes |
| `code_a_fig09_gmm_estimates.pdf` | 1000 MC estimates of rho (true 0.4): black rho_1 (8a), green second-order rho^(2), blue optimal GMM rho^opt |
| `code_a_fig10_gmm_variances.pdf` | n × variance of the lag-l correlation estimator for l = 1..10 (solid MC, dotted asymptotic), dashed = rho^opt, diamond = rho^(2) |
