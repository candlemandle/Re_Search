# Response to review — Code A

Scope: branch `code-a` (the reviewer's reference copy `Re_Search-code-a/`). The integrated
`RS 2` tree, `docs/review_report.md` and the reviewer logs were not available to me; this
response works from the revision instructions and reproduces every quoted number here.
I do not change the review status; everything below is evidence for re-review.

Environment of all new runs: R 4.6.1 (aarch64-apple-darwin), testthat 3.x, base R otherwise.
Every script run now appends a provenance row (time, script, mode, R version, git commit,
"code modified" flag, seed, MD5 of the processed panels) to `results/logs/code_a_runs.csv`.

| item | verdict |
|---|---|
| A1 Table 1 standard errors | confirmed (claim wording); numerical cause identified, not a formula defect |
| A2 non-rejection vs latent reversibility | confirmed (handoff wording); test unchanged |
| A3 eta conventions | partly_applicable: already a documented convention difference; labels added |
| A4 MC agreement claims | confirmed: claim was unquantified; criterion added, exceptions listed |
| A5 data preparation / missing days | partly_applicable: rules documented and quantified; snapshot now protected |

---

## A1. Table 1 standard errors

**Verdict: confirmed** for the claim; the cause is a finite-sample convention, not a defect.

**Reproduction.** `Rscript scripts/03b_check_table01_se.R` → `results/tables/code_a_table01_se_check.csv`.

**Evidence.**
* All 25 point estimates equal Table 1 after rounding to 4 decimals (unchanged).
* With our SEs (closed-form limits (9), (14), (27), N = 5007 increments), **23/25** SEs equal the
  paper after rounding — exactly the reviewer's finding. Exceptions:
  eta(AAPL, ALD) 0.03225229 vs 0.0322 and eta(AAPL, AXP) 0.03255567 vs 0.0325.
* Truncation is not the cause: R = 5000, 20000 and 10^6 series terms differ by at most
  2.9e-12 relatively (`convergence_max_rel_diff`). Truncation (instead of rounding) of the
  paper's numbers is also not the cause (only 14/25 would match).
* The authors' own `BYZinference.R` functions (`varHest`, `varrho2`, `vareta`: finite sums with
  (1 − r/n) weights, truncation n = 2000) give 24/25 with N = 5007 and **25/25 with N = 5008**
  (the number of observed levels on the five-stock common sample). This is a post hoc but
  complete explanation; it is consistent with, not proven to be, the authors' computation.
* The same mechanism explains 2 of 3 deterministic SE cells in Tables 4–5 that differ in the 4th
  decimal (sigma2_2, Delta = 1/250, n = 1000: ours 0.274151, authors' formula 0.274113, paper
  0.2741). The third (n = 500: 0.3876 in Table 4) is printed as 0.3877 in Table 5 for the
  identical quantity, i.e. an inconsistency inside the paper.

**Changes.** No formula change: our SEs are the limits stated in Theorems 3.1–3.3; the finite-sum
variant is asymptotically equivalent and is reported side by side. Stale claims of "all 25 SEs
exact" corrected in `docs/pr_code_a.md`, `docs/code_a_handoff.md`, `data/README.md`,
`docs/paper_coverage.csv`.

**Remaining limitation.** The N = 5008 / finite-sum explanation is inferred, not documented by the authors.

## A2. Non-rejection is not proof of latent time reversibility

**Verdict: confirmed** for the narrative in `docs/code_a_handoff.md` ("which justifies the
time-reversible mfBm", "Most pairs are consistent with eta = 0, as the paper states");
the test itself was not questioned and is unchanged.

**Reproduction.** `Rscript scripts/03_estimate_parameters.R full` →
`code_a_time_reversibility_{dj30,mag7,summary}.csv`.

**Evidence (unchanged counts).** DJ30 435 pairs: 27 rejected at 1%, 66 at 5%; Table 1 subset:
3/10 at 1%, 4/10 at 5%; Mag7: 1/10 at both levels; 0 invalid tests (H_hat outside (0, 3/4)).

**Changes.**
* Handoff rewritten: eta = 0 is a *maintained forecasting specification*, partly supported by
  individual non-rejections; these are dependent pairwise tests, not a calibrated global test;
  measurement noise biases eta_hat towards 0 (Table 12: −0.37 at eta = 0.5) and lowers power;
  results concern observed log RV, not latent volatility (paper footnote 7, Appendix E.4).
* Summary table gains an `interpretation` column and separately labelled **exploratory**
  multiplicity columns (Holm / Benjamini–Hochberg at 5%: DJ30 3 / 8 rejections; Table 1 subset
  2 / 3; Mag7 1 / 1). The unadjusted counts of the paper's procedure are kept as the primary output.

## A3. eta conventions

**Verdict: partly_applicable.** The difference was already identified and documented as a
convention difference (not an inconsistency) before the review; what was missing were labels in
the exported tables and an end-to-end orientation test.

**Evidence** (`Rscript scripts/00_check_author_code.R` → `code_a_author_code_comparison.csv`):
* Authors' simulator `simMFBM` with eta = +0.5 and printed (8b) recover +0.507 (50 paths);
  the Appendix C.2 version on the same paths gives −0.507.
* Our simulator, built from eq. (4), with eta = +0.5 and the Appendix C.2 version recovers +0.502.
* Our printed-(8b) implementation equals the authors' code on the same path up to their
  `mean()` normalisation (|diff| 1e-3 at n = 1000).
* So: printed (8b) + simMFBM form one consistent convention, eq. (4) + Appendix C.2 the other;
  `eta_mm(x1, x2, "printed") == eta_mm(x2, x1, "proof")`. The paper's Table 1/14 cell in row i,
  column j equals our eq. (4) eta of (column, row); 25/25 Table 1 signs agree under this mapping.
* Forecasting uses eta = 0 and the test (11) is two-sided, so neither depends on the sign.

**Changes.** `eta_convention` column in `code_a_table01_estimates.csv`,
`code_a_table13_14_dj30.csv` and the reversibility tables; new test "signed eta keeps its
orientation through the empirical pipeline" (simulated eta = +0.5 → `estimate_panel`,
`time_reversibility_table`, printed/proof identity). No sign was changed anywhere.

## A4. Monte Carlo: execution vs validated replication

**Verdict: confirmed.** "Within Monte Carlo error" in the PR/handoff had no stated criterion
and, once a criterion is applied, is too strong.

**Provenance of the full tables.** `code_a_mc_*.csv` in `results/tables/` come from one full
run on 2026-10-04 20:44–20:58 (log `results/logs/code_a_mc_full_stdout.log`) with the original
`code-a` code. They were **not** rerun for this response: the simulation and estimation code
(`R/covariance.R`, `R/estimators.R`, the experiment functions in `R/monte_carlo.R`) is unchanged
since (`git diff` touches only appended reporting functions). Replications: Tables 4–6 and 12:
1000; Table 7: 5000 per cell; Table 20: 1000 (+10 000 independent paths for rho_lag/rho^(2)
variances); seeds per scenario and per chunk of 50 (`mc_seed`, `mc_run`), independent of core count.

**Criterion** (`mc_agreement()` in `R/monte_carlo.R`; `Rscript scripts/04b_compare_mc_with_paper.R full`
→ `code_a_mc_agreement.csv`, `code_a_mc_agreement_summary.csv`): for each paper cell
z = (ours − paper)/sqrt(se_ours² + se_paper²), both MC standard errors estimated from the
respective replication counts (bias: sd/√R; sd: sd/√(2(R−1)); RMSE: delta method; rates:
binomial; Table 20 variance: var·√(2/(R−1))); agreement = |z| ≤ 2 or difference within the
paper's rounding. Cells share paths and are numerous, so about 5% of |z| > 2 is expected even
under exact agreement; this is a screening criterion, not a joint test.

**Results.**

| table | MC cells | \|z\| ≤ 2 | \|z\| > 3 | max \|z\| |
|---|---:|---:|---:|---:|
| 4 | 72 | 66 | 1 | 3.06 |
| 5 | 72 | 67 | 0 | 2.84 |
| 6 | 96 | 90 | 1 | 4.96 |
| 7 | 32 | 26 | 2 | 4.13 |
| 12 | 32 | 30 | 0 | 2.50 |
| 20 | 12 | 11 | 0 | 2.23 |
| **total** | **316** | **290 (91.8%)** | 4 | |

Asymptotic SEs (deterministic): 45/48 equal after rounding (see A1 for the three cells).

Exceptions worth naming:
* **Table 7, eta = 0.3–0.4**: our power is higher by 0.8–2.2 percentage points in 5 cells
  (n = 500 and n = 1000, z = 2.7–4.1); one cell (n = 500, eta = 0.1, 1%) is lower by 1.2 points
  (z = −2.6). Size is fine (5%: 0.052 / 0.053 vs 0.058 / 0.053). Tested hypothesis: counting replications with H_hat1 < 0 (1.0–1.6%) as
  non-rejections (columns `*_all`) — rejected, it worsens agreement at large eta (z up to −8).
  Rates over valid replications remain the primary output; the difference is **unresolved**.
* **Table 6, AC rho bias, rho = 0.2, n = 500**: −0.0059 vs +0.0045 (z = −5.0), although our AC
  estimator equals the authors' `estMFBM` to machine precision on identical paths; the other 47
  AC cells agree. **Unresolved**; possibly a simulator difference in the authors' run.
* Sd/RMSE of sigma^2 (Tables 4–5, Delta = 1/250): |z| 2.0–3.1. The se of an SD uses a normal
  approximation, which understates it for the skewed sigma^2_hat, so these z overstate the gap.
* eta = 0.65 cells use the approximate circulant embedding (`approximate_simulation = TRUE`,
  relative clipped eigenvalue −1.6e-5); both are at power 1.000 and agree.

**Changes.** Criterion, scripts and outputs above; PR/handoff/coverage now cite the cell
counts instead of a blanket claim; new test of the criterion (it also exposed and fixed an
edge case in `mc_agreement` with an empty reference subset).

**Remaining work.** A fresh full Monte Carlo run by the reviewer is still needed to certify
paper-sized execution independently; command `Rscript scripts/04_run_monte_carlo.R full`
(13–25 min on 8 cores).

## A5. Data preparation and missing days

**Verdict: partly_applicable.** The rules were documented, but (i) the raw/processed snapshot
was not protected or checksummed and (ii) the explanation of the paper's rule was stated more
firmly than the evidence supports.

**Facts.**
* Every rule collapses missing days: consecutive *available* observations are one step Delta.
  This is a deliberate replication convention, not an exactly equally spaced model.
  Bridged increments per series: 0–16 over 2005–2025 (`code_a_missing_day_gaps.csv`; V starts 2008-03-19).
* Samples: Table 1 = five-stock common sample (5008 dates); Tables 13–14 = H on own dates,
  rho/eta on pairwise common dates; Mag7 = five-stock common sample.

**Sensitivity** (`Rscript scripts/03c_missing_day_sensitivity.R` → `code_a_missing_day_sensitivity.csv`),
share of cells equal to the paper at 4 decimals / max |diff|:

| rule | Table 13 H | Table 13 rho | Table 14 eta |
|---|---|---|---|
| own dates / pairwise common dates (used) | 70% / 0.0105 | 57% / 0.0038 | 55% / 0.0091 |
| hybrid: Table 1 stocks on five-stock sample | 87% / 0.0105 | 83% / 0.0008 | 83% / 0.0014 |
| common sample of all 30 | 3% / 0.028 | 0% / 0.073 | 1% / 0.067 |
| calendar grid, spanning increments dropped | 17% / 0.012 | 6% / 0.0037 | 3% / 0.010 |

Table 1 H: common sample 5/5 equal (max diff 4.6e-5); own dates 0/5 (max diff 0.0021).
The hybrid rule fits best but still not exactly (residual: JNJ H, some GS/PG/WAG pairs); it is
reported as a **hypothesis**, and the used rule is kept because it is a single transparent rule
(changing it to chase agreement would be tuning on the published table).

**Snapshot.** Processed panels are committed with `data/processed/MD5SUMS`; the raw files used
(downloaded 2026-10-04; ticker, PERMNO, bytes, MD5, first/last day) are recorded in
`data/raw_snapshot.csv`. `scripts/02_prepare_data.R` now refuses to overwrite the processed
snapshot when rebuilt panels differ (Risk Lab extends and may revise history) unless
`--refresh-snapshot` is given; a test checks the processed files against `MD5SUMS`. Raw files
themselves are not committed (third-party data, ~20 MB); with the processed snapshot every
Code A result is reproducible without them.

---

## Commands run for this response (2026-10-08)

```bash
Rscript scripts/00_check_author_code.R
Rscript scripts/02_prepare_data.R                 # snapshot guard: panels identical
Rscript scripts/03_estimate_parameters.R full
Rscript scripts/03b_check_table01_se.R
Rscript scripts/03c_missing_day_sensitivity.R
Rscript scripts/04b_compare_mc_with_paper.R full  # on the saved 2026-10-04 full MC tables
Rscript scripts/04_run_monte_carlo.R smoke        # outputs in results/smoke/
Rscript -e 'testthat::test_dir("tests/testthat")' # 30 test blocks, all pass
```

The integrated project's `scripts/check_data.R` and `scripts/run_tests.R` belong to Code C and do
not exist on this branch.

## Changed files

* New: `scripts/03b_check_table01_se.R`, `scripts/03c_missing_day_sensitivity.R`,
  `scripts/04b_compare_mc_with_paper.R`, `data/raw_snapshot.csv`, `data/processed/MD5SUMS`,
  `docs/response_to_review_A.md`.
* Modified: `R/code_a_config.R` (MD5 helper, run log), `R/monte_carlo.R` (`mc_agreement`,
  `mc_agreement_summary`), `scripts/02_prepare_data.R` (snapshot guard),
  `scripts/03_estimate_parameters.R` (eta convention labels, exploratory multiplicity,
  interpretation), `scripts/04_run_monte_carlo.R` (agreement + run log at the end),
  `scripts/00_check_author_code.R` (run log), `tests/testthat/test-estimators.R` (3 tests),
  `README.md`, `data/README.md`, `docs/pr_code_a.md`, `docs/code_a_handoff.md`,
  `docs/paper_coverage.csv`.
* Regenerated outputs: `code_a_table01_*`, `code_a_table13_14_dj30.csv`,
  `code_a_time_reversibility_*`, `code_a_missing_day_*`, `code_a_mc_agreement*`,
  `code_a_author_code_comparison.csv`, `results/logs/code_a_runs.csv`.
