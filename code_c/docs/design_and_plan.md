# Code C research extension — design and implementation plan

Date: 2026-10-09. The user authorizes choosing scientifically justified extra experiments and implementing them, while preserving reviewed A/B source.

## Goal and architecture

Add a self-contained `code_c/` overlay to the reviewed B v2 snapshot (which already contains reviewed A v2). No A/B-owned file is changed. Source B's public loader and reuse its estimators, admissibility, trading calendar and exact mfBm forecaster. All new functions use `c2_` names; outputs live in `results/code_c_v2/<mode>/`. The local reviewed snapshot is an input dependency, never copied into the commit payload.

## Scientific choices

1. Genuine stationary, diagonal-kappa mfOU driven by eta=0 mfBm versus B's mfBm, on identical information sets (AAPL through AAPL/ALD/AMGN/AXP/BA), panel trading dates and horizons 1,2,3,4,5,10,15,20. The primary window is 500 days. Section 6, PDF pp.23–24 explicitly proposes mfOU forecasting.
2. Prespecified window sensitivity: 250 and 1000 days for dimensions 1 and 5. This is a sensitivity experiment motivated by the task brief and the paper's two-year baseline, not an experiment reproduced from the paper. Comparisons intersect forecast keys across windows.
3. Known-versus-estimated parameter simulation, extending Appendix E.3/Table 11. Simulate both mfBm and stationary mfOU; never give oracle parameters to real-data forecasts. Save Monte Carlo standard errors, parameter errors and failures.
4. Retain the asset-selection question with fixed/correlation/H-gap strategies, but use the reviewed trading-calendar predictors and gap-aware training estimates. Its score remains a disclosed heuristic, not a theorem.

An AR(1)/Euler-filter surrogate would be easier but would not answer the requested mfOU question. Joint dense multivariate maximum likelihood would increase cost and complexity. We choose exact Gaussian conditioning with a non-optimal finite-lag moment estimator. The mfOU process with different kappas need not itself be time reversible even though its noise is.

## Mathematics and estimation

The stationary solution is X_i(t)=mu_i+sigma_i integral(-infinity,t) exp[-kappa_i(t-s)] dB_i^H(s). Covariance is evaluated using a one-dimensional integration-by-parts identity valid for every 0<H_i<1, including rough H<1/2. Validate against analytic Brownian OU, exact variance, independent adaptive integration, positive definiteness and numerical quadrature refinement.

Use finite-lag squared differences and expected sample-centered variance to estimate H,kappa,sigma jointly per component. Missing observations retain their real trading positions. Use several training-only optimizer starts, finite bounds and recorded local moment-Jacobian diagnostics. Estimate driver rho from genuine one-day cross-increments, then use B's noise admissibility/shrinkage rule. Estimate means by GLS using the same conditioning matrix. Conditional log variance excludes parameter uncertainty; simulation reports its effect separately.

## Implementation tasks and checks

1. Write failing covariance/forecast tests; implement `R/mfou.R`; check OU limit, rough covariance, signs, translation, log-normal conversion and population moments.
2. Write failing reviewed-API/calendar tests; implement `R/config.R`, `R/experiments.R` and C-only entrypoints. Check direct B equality, future mutation invariance, shared outcomes, missing targets and origin attrition.
3. Implement checkpoint identities/resume and artifact-only output; verify a resumed run equals an uninterrupted run and code/config/data changes invalidate caches.
4. Implement descriptive/paired metrics, prespecified periods and honest uncertainty availability. Dense/sparse calendar metadata must accompany bootstrap. Smoke never exports inferential intervals.
5. Implement known/estimated simulations and numerical sensitivity; validate zero correlation/OU anchors and MCSE bookkeeping.
6. Execute the reviewed A/B tests in a disposable merged tree, new C tests, real 500-row smoke and a measured substantive study if runtime permits. Do not launch an unbounded full experiment. Report actual sample and failures without superiority claims from smoke.
7. Create Russian/English READMEs, simple explanation, mathematical/source notes, actual results interpretation and commit manifest. Verify all reviewed snapshot hashes unchanged. Prepare an additive upload payload; do not invent Git history or promise zero conflicts with an unseen remote branch.

No Git repository or remote URL is available initially. User was asked for the remote and branch while independent work continues. Documentation and implementation stay within the additive folder.
