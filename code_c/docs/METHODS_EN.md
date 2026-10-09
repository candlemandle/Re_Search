# Code C: models, experiments and limits

## What this extension asks

The article forecasts realized volatility from its own history and the histories of other assets. Its mfBm model captures rough paths and cross-asset dependence. Code C asks whether adding a tendency to return toward a long-run level helps forecasting. It also checks how the length of the training window and parameter estimation affect the result. Improvement is an empirical question, not an assumption.

The extension is an additive `code_c/` folder. It uses the reviewed B v2 interfaces, which include A v2. It does not change their source. Mathematical functions may be involved, but the workflow stays simple: load the reviewed functions and data, fit on a past window, forecast, save outcomes, compare losses.

## Which experiments are justified

| Experiment | Question | Basis and interpretation |
|---|---|---|
| mfBm versus genuine stationary mfOU | Does mean reversion help, particularly at 10, 15 and 20 trading days? | Article §6, PDF pp. 23–24 explicitly proposes mfOU modelling and forecasting. A positive result must be demonstrated on paired outcomes. |
| 250/500/1000-day windows | Is the answer sensitive to the amount of history available for estimation? | The article uses approximately two years; the supplied research task asks for one/two/four years. This sensitivity experiment extends the paper; it is not a reproduced published result. |
| Known versus estimated parameters under two simulated processes | How much forecasting accuracy is lost because parameters must be estimated? | Extends Appendix E.3, Table 11, PDF p. 48. Both mfBm and mfOU data-generating processes are used. Known parameters are available only in simulation. |
| Existing fixed/correlation/H-gap asset selection | Does the choice of additional assets matter? | Motivated by the article's correlation/Hurst-gap results, §4 and Appendix E.3. The H-gap score is a heuristic, not a proved optimum. |

The primary mfOU comparison uses 500 trading-calendar rows and the article's nested sets: AAPL; AAPL/ALD; then AMGN, AXP and BA, up to five assets. Forecast horizons are `1, 2, 3, 4, 5, 10, 15, 20` trading days. Window sensitivity uses dimensions one and five at 250 and 1000 days, alongside the primary 500-day result. Window comparisons use the intersection of forecast keys, so a change in sample coverage is not silently treated as a model improvement.

Full-period and prespecified subperiod losses provide additional outcomes. The article's §5.2, PDF p. 21, motivates the boundary `2017-04-11` and period ending `2021-07-30` through convergence of estimated Hurst exponents. These are descriptive regime comparisons, not a new structural-break estimator. Exact sample dates, boundaries and usable counts must accompany every result.

## Baseline and new model

The baseline calls B's reviewed time-reversible mfBm forecast directly. Its component-wise Hurst exponents, scales and noise correlations are estimated from the training data; admissibility and handling of trading-calendar gaps follow B. The window is anchored to its first observed value, as in B's implementation. Absolute calendar labels are not a new baseline parameter.

For mfOU, component $i$ has the stationary solution

$$
X_i(t)=\mu_i+\sigma_i\int_{-\infty}^{t}
e^{-\kappa_i(t-u)}\,dB_i^{H_i}(u),\qquad \kappa_i>0.
$$

The driver $B^H$ is an admissible multivariate fractional Brownian motion with component-wise $H_i\in(0,1)$, driver correlations $\rho_{ij}$, and $\eta_{ij}=0$. The drift matrix is diagonal. The integration over the infinite past defines a stationary Gaussian process; an AR(1) or Euler-filter approximation is not substituted for it.

Different mean-reversion rates can create asymmetric cross-lag covariance even with reversible driving noise. The precise description is **mfOU driven by time-reversible mfBm**. It is generally not itself time reversible when the rates differ. The driver correlation $\rho_{ij}$ is also generally different from the contemporaneous correlation of the mfOU levels.

With time measured in years and a daily step of $1/252$, $\kappa_i$ has units per year. $\log(2)/\kappa_i$ is the deterministic drift half-life. It is not an autocorrelation half-life: fractional OU is not Markovian when $H_i\ne1/2$, and rough-process covariances may become negative.

## Model differences

| Property | Reviewed mfBm baseline | Stationary mfOU extension |
|---|---|---|
| Level dynamics | Nonstationary fractional Gaussian levels, locally anchored | Stationary fractional Gaussian levels with diagonal mean reversion |
| Parameters | H, σ, ρ; η fixed to zero | H, σ, ρ, κ and µ; driver η fixed to zero |
| Starting distribution | First observed window value is the anchor | Stationary covariance from the infinite past |
| Estimation | Reviewed B moment estimators and admissibility | Custom variogram/centered-variance moments, multiple starts, driver admissibility, mean GLS |
| Forecast | B's Gaussian conditional forecast | Gaussian conditional forecast from stationary mfOU covariance |
| Return to RV | exp(m+v/2) | exp(m+v/2), with plug-in variance |
| Reversibility | Reversible specification | Reversible driver; unequal κ can make the resulting process non-reversible |

## Exact covariance and its numerical evaluation

Use the convention $C_{ij}(\tau)=\operatorname{Cov}(X_i(t+\tau),X_j(t))$, and write $a=H_i+H_j$. An integration-by-parts identity gives, for $\tau\ge0$,

$$
C_{ij}(\tau)=\frac{\rho_{ij}\sigma_i\sigma_j}{2}
\left[\frac{\kappa_j}{\kappa_i+\kappa_j}E|\tau+V|^a
+\frac{\kappa_i}{\kappa_i+\kappa_j}E|\tau-U|^a-\tau^a\right],
$$

where $U\sim\mathrm{Exp}(\kappa_i)$ and $V\sim\mathrm{Exp}(\kappa_j)$. The expectations are evaluated through an upper incomplete gamma term and one bounded integral:

$$
\begin{aligned}
P_-&=a\kappa_j^{-a}e^{x_j}\Gamma(a,x_j),\quad x_j=\kappa_j\tau,\\
M_-&=\Gamma(a+1)\kappa_i^{-a}e^{-x_i}
-\tau^a\int_0^1 e^{-x_i(1-z^{1/a})}\,dz,\quad x_i=\kappa_i\tau,\\
C_{ij}(\tau)&=\frac{\rho_{ij}\sigma_i\sigma_j}
{2(\kappa_i+\kappa_j)}(\kappa_jP_-+\kappa_iM_-).
\end{aligned}
$$

At zero use

$$
C_{ij}(0)=\frac{\rho_{ij}\sigma_i\sigma_j\Gamma(a+1)}
{2(\kappa_i+\kappa_j)}
(\kappa_i^{1-a}+\kappa_j^{1-a});
\qquad C_{ij}(-\tau)=C_{ji}(\tau).
$$

This avoids cutting off a slowly converging spectral tail when $H$ is small. Gamma probabilities are evaluated on the log scale. Numerical validation includes the ordinary OU limit $H_i=H_j=1/2$, exact marginal variance, rough $H=0.1$, independent integration, covariance symmetry/positive definiteness, and quadrature refinement. Negative off-diagonal or long-lag covariance is not automatically an error and must not be clipped to zero.

## Training-only finite-lag estimation

The estimator is a custom nonoptimal moment estimator, not maximum likelihood and not a reproduction of the full joint GMM estimator in the mfOU literature.

For each component, compare observed squared differences at prespecified lags with

$$
E[(X_i(t+\ell\Delta)-X_i(t))^2]
=2\{C_{ii}(0)-C_{ii}(\ell\Delta)\}.
$$

The other moment is the sample-centered variance, calculated with denominator $n$: $\widehat V_i=n^{-1}\sum(Y_i-\overline Y_i)^2$. Subtracting a sample mean changes its expectation in a correlated finite window. If $S_i$ is the covariance matrix at the $n$ actually observed training positions,

$$
E[\widehat V_i]=\frac{\operatorname{tr}(S_i)}{n}
-\frac{\mathbf1^\top S_i\mathbf1}{n^2}.
$$

Thus the fit does not equate a finite-window centered variance directly to the population marginal variance. For fixed $H_i,\kappa_i$, the model moments are linear in $\sigma_i^2$, which is profiled out. If $g_k$ is a unit-scale model moment and $e_k$ its positive empirical counterpart, write $z_k=g_k/e_k$. The relative-error objective has the analytic scale estimate $\widehat\sigma_i^2=\sum z_k/\sum z_k^2$.

The default squared-difference lags are `1, 2, 5, 10, 20, 60` trading days. Several starting values optimize $H_i$ and $\log\kappa_i$ within finite bounds. The implementation defaults use $H_i\in[0.02,0.8]$, annual $\kappa_i\in[0.02,50]$ and rate starts `0.1, 1, 10`; these are numerical design choices, not restrictions asserted by the paper. Record convergence, boundary solutions, disagreements between starts, and local moment-Jacobian diagnostics. These checks do not prove global identification; small mean reversion can be weakly identified in a short window.

Driver correlation is calibrated from genuine same-day cross-increments using

$$
E[\Delta X_i\,\Delta X_j]
=2C_{ij}(0)-C_{ij}(\Delta)-C_{ij}(-\Delta).
$$

Its dependence on mfOU parameters is retained. Missing observations are excluded from the relevant moment pairs while their original trading positions remain in the time axis. B's mfBm noise-admissibility check and correlation shrinkage are reused; shrinkage is recorded. No realized outcome or future observation enters parameter fitting, asset ranking or starting-value selection.

## Conditional forecast and parameter uncertainty

Stack training observations in B's time-major order. Let $S$ be their mfOU covariance matrix, $D$ the matrix assigning each observation to its component mean, and $y$ the observation vector. Estimate means jointly by GLS:

$$
\widehat\mu=(D^\top S^{-1}D)^{-1}D^\top S^{-1}y.
$$

For a future target with covariance vector $c$,

$$
m=\widehat\mu_i+c^\top S^{-1}(y-D\widehat\mu),\quad
v=C_{ii}(0)-c^\top S^{-1}c,\quad
\widehat{RV}=\exp(m+v/2).
$$

Use Cholesky solves, not an explicit matrix inverse. This conditions on the available history, not just the last value. The variance $v$ treats the fitted means and covariance parameters as known. It excludes their estimation uncertainty. The log-normal adjustment is therefore a plug-in model adjustment, not a complete predictive-uncertainty calculation. Known/estimated simulation exposes the resulting performance gap; it does not automatically correct real-data intervals.

## Evaluation and interpretation

Compare models on identical target, origin, target-date and horizon keys with the same realized outcomes. The team data use annualized `qmle_trade` volatility, not its square; preserve B's units and loss definitions. Report MSFE, RMSFE and QLIKE, absolute losses and percentage change relative to mfBm. QLIKE is $r-\log r-1$, where $r=RV/\widehat{RV}$. A cumulative loss-difference plot shows when gains or losses accumulate; it is not a significance test.

Bootstrap comparisons retain pairing and acknowledge horizon overlap, calendar gaps and block choice. On a daily trading grid, missing pairs remain at their original positions as NA; the same moving calendar-block indices resample both loss series, and means use only shared finite pairs. Inference requires at least 90% paired calendar coverage, a median origin step of one and at least ten blocks of observed pairs. These are prespecified availability guards, not a guarantee of independence. The requested block is at least the horizon; sensitivity uses 20/40/60 trading days. Sparse grids and smoke receive no intervals. Intervals are exploratory, pointwise and have no multiple-comparison correction. Report coverage and failures beside losses; a failed mfOU forecast never becomes an mfBm prediction.

There is no theorem guaranteeing that mfOU beats mfBm. A negative result may reflect weak mean reversion, limited identification, additional estimation error or a useful local mfBm approximation. A better average loss alone does not establish statistical or economic superiority. Runtime or numerical tolerances affect feasibility and numerical accuracy; they do not by themselves improve empirical forecasting.

## Source map

- Local article: `04_Bibinger_Yu_Zhang_2026_JBES.pdf`, printed/PDF pp. 14–15 (§4, conditional Gaussian forecasting); pp. 20–21 (§5.2, assets, horizons, rolling windows and subperiods); pp. 23–24 (§6, proposed mfOU/estimation/structural-break research); p. 48 (Appendix E.3, Table 11); pp. 56–57 (Appendices F.4–F.5, MCS and QLIKE). The local version is authoritative for these page anchors. [Authors' public preprint](https://arxiv.org/abs/2504.15985).
- Dugo, Giorgio and Pigato, *The multivariate fractional Ornstein–Uhlenbeck process*: [definition, covariance and time-reversibility, §§2.2–2.4](https://arxiv.org/html/2408.03051v2); [published DOI](https://doi.org/10.1016/j.spa.2025.104814). The bounded-integral identity above is derived from the stationary definition, not claimed as a copied published algorithm.
- Dugo, Giorgio and Pigato, *Multivariate Rough Volatility*: [§§2–3, covariance-matching estimation; §4.4, slow mean reversion](https://arxiv.org/html/2412.14353v3). This motivates the estimation family; Code C's reduced estimator differs from their full procedure.
- Wang, Xiao and Yu, *Modeling and forecasting realized volatility with the fractional Ornstein–Uhlenbeck process*: [published paper](https://doi.org/10.1016/j.jeconom.2021.08.001). Univariate background; it does not establish superiority of this multivariate implementation.

All active nested sets share the largest-set past-only eligibility calendar within each window, following B. The empirical contrast changes the model, estimator and initial distribution together; it cannot isolate a causal contribution of mean reversion alone. Known-parameter simulation clarifies estimation effects but does not identify the cause of an empirical gain.

Within each model, conditioning and moment fitting use jointly observed rows across its assets, matching B; partially observed rows are not used as extra conditioning observations. Their positions remain in the trading time axis.
