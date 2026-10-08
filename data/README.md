# Data

## Source

Daily realized volatility estimates from **Risk Lab** (Dacheng Xiu, Chicago Booth),
<https://dachxiu.chicagobooth.edu/#risklab> — the source cited in footnote 5 of the paper.

The website's chart loads one file per security from the public endpoint

```text
https://dachxiu.chicagobooth.edu/data.php?ticker=<PERMNO>
```

`scripts/01_download_data.R` downloads these files into `data/raw/risklab_<TICKER>.txt`
(one request per ticker, 1 s pause, ~0.7 MB per file, ~20 MB total). Raw files are not
committed (`data/raw/` is in `.gitignore`); rerun the script to recreate them.
Risk Lab updates its history daily, so later downloads contain more days; the sample
is cut to the paper's period in `scripts/02_prepare_data.R`.

## File format

Six header lines (ticker, PERMNO, company name, number of days, first day, last day),
then one line per trading day with 12 space-separated fields:

| field | meaning |
|---|---|
| `permno`, `date` | CRSP PERMNO, `YYYYMMDD` |
| `qmle_trade`, `ma_trade`, `ci_trade` | QMLE volatility from all trades, selected MA order of the noise, CI half-width |
| `rv5_trade`, `rv15_trade` | 5- and 15-minute realized volatility from trades |
| `qmle_quote`, `ma_quote`, `ci_quote` | same from mid-quotes |
| `rv5_quote`, `rv15_quote` | 5- and 15-minute realized volatility from mid-quotes |

Values are **annualized volatilities** (square root of annualized variance).

## Which measure the paper uses

We use **`qmle_trade`** (Risk Lab's headline estimator). This was identified empirically:
with `qmle_trade`, all 25 point estimates of Table 1 (H, rho, eta) are reproduced to the
4th decimal (23 of the 25 standard errors as well; see `docs/response_to_review_A.md`, A1); the 5- and 15-minute RVs give H lower by 0.05–0.11 and
rho lower by 0.08–0.17 (see `results/logs/code_a.log` and the PR notes).

The modelled process is `B_t = log(vol_t)`. H, rho, eta and all their standard errors
are invariant to using log variance (`2 log vol`) instead; only sigma^2 changes by a factor 4.

## Tickers

The paper lists DJ30 stocks under the *historical* Risk Lab names of each PERMNO
(Table 13). Mapping used (`R/code_a_config.R`):

| paper | today | PERMNO | | paper | today | PERMNO |
|---|---|---|---|---|---|---|
| AAPL | AAPL | 14593 | | INTC | INTC | 59328 |
| ALD | HON (AlliedSignal → Honeywell) | 10145 | | JNJ | JNJ | 22111 |
| AMGN | AMGN | 14008 | | JPM | JPM | 47896 |
| AXP | AXP | 59176 | | KO | KO | 11308 |
| BA | BA | 19561 | | MCD | MCD | 43449 |
| BEL | VZ (Bell Atlantic → Verizon) | 65875 | | MMM | MMM | 22592 |
| CAT | CAT | 18542 | | MRK | MRK | 22752 |
| CHV | CVX (Chevron) | 14541 | | MSFT | MSFT | 10107 |
| CRM | CRM | 90215 | | NIKE | NKE | 57665 |
| CSCO | CSCO | 76076 | | PG | PG | 18163 |
| DIS | DIS | 26403 | | SPC | TRV (St Paul → Travelers) | 59459 |
| GS | GS | 86868 | | UNH | UNH | 92655 |
| HD | HD | 66181 | | V | V (from 2008-03-19) | 92611 |
| IBM | IBM | 12490 | | WAG | WBA (Walgreens) | 19502 |
| | | | | WMT | WMT | 55976 |
| | | | | XOM | XOM | 11850 |

Magnificent 7 subset (Appendix F.3): AAPL 14593, AMZN 84788, FB (META) 13407,
GOOG (class C, from 2014-04-02) 14542, MSFT 10107.

## Samples and processed files

| file | content |
|---|---|
| `processed/dj30_logvol.csv` | log volatility, 30 stocks, 2005-01-05 – 2025-01-14 (5019 dates; NA when a stock has no value) |
| `processed/mag7_logvol.csv` | log volatility, 5 stocks, 2014-03-27 – 2025-01-14 (two years before the first forecast on 2016-03-29) |
| `../results/tables/code_a_data_summary.csv` | per-series coverage, number of observations, moments |

Cleaning: values below 1e-6 (Risk Lab's own cut-off) and non-finite values are dropped;
no other filtering or outlier treatment.

## How missing days are handled

* **H, sigma^2** of a single series: increments between consecutive available days of that series.
* **rho, eta** of a pair: increments on the dates where both series are observed.
* **Table 1**: the five stocks are treated as one 5-dimensional mfBm on their common dates
  (5008 dates, 5007 increments). This reproduces the Table 1 point estimates; estimating H on each
  stock's own dates changes H by up to 0.002.
* **Tables 13–14** (30 stocks): rules above. 70% of H, 57% of rho and 55% of eta values
  match the paper to 4 decimals; the largest deviations are 0.0105 (H of JNJ), 0.0038 (rho)
  and 0.0091 (eta), all below one standard error. The paper copies the five Table 1 values
  of H into Table 13 and appears to use the Table 1 common sample for pairs involving those
  five stocks (this hybrid rule raises the exact-match rate to 83%); we keep the simpler
  pairwise rule. Which rule the paper used is not stated, so these are hypotheses; the
  comparison of four rules is in `results/tables/code_a_missing_day_sensitivity.csv`
  (`scripts/03c_missing_day_sensitivity.R`).
* All rules collapse missing days: consecutive *available* observations are treated as one
  sampling step Delta (a replication convention, not an exactly equally spaced model).

## Snapshot

The processed panels used for every reported result are fixed by `processed/MD5SUMS`; the raw
files they were built from (downloaded 2026-10-04) are listed with size and MD5 in
`raw_snapshot.csv`. `scripts/02_prepare_data.R` stops instead of overwriting the processed
snapshot if a fresh download changes the panels; use `--refresh-snapshot` only deliberately.

The dates on which a series is missing are listed by `n_obs` / `share_of_panel_dates` in
`code_a_data_summary.csv` (V starts in March 2008; every other stock misses at most 43 days).
