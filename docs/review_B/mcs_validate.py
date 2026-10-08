"""B3: validate Code B's MCS p-values (R, metrics.R::mcs_pvalues) independently.

1. exact check: a separate NumPy implementation of the Hansen-Lunde-Nason (2011)
   T_max elimination (statistic of each model relative to the average of the models
   still in the set, bootstrap variance, p-value of an eliminated model = running
   max of the test p-values) fed with the *same* resamples exported from R
   (B = 1000). Agreement must be exact up to floating point.
2. external check: arch 8.0.0 (K. Sheppard) MCS(method="max", bootstrap="circular",
   block_size=20, reps=5000) with its own random resamples, 5 seeds. Agreement is
   statistical: compare with the spread of our own p-values over 5 seeds.

Usage: <python with arch> -I docs/review_B/mcs_validate.py <export dir>
"""
import csv, sys, os, warnings
import numpy as np

out = sys.argv[1]
rows = list(csv.DictReader(open(os.path.join(out, "r_pvalues.csv"))))
cells = sorted({r["cell"] for r in rows}, key=lambda c: int(c[4:]))


def tmax_mcs(L, idx):
    """HLN T_max MCS given resample indices idx (B x n, 0-based)."""
    n, M = L.shape
    alive = list(range(M))
    pval = np.full(M, np.nan)
    pmax = 0.0
    while len(alive) > 1:
        Lk = L[:, alive]
        d = Lk - Lk.mean(axis=1, keepdims=True)          # d_i. = L_i - mean_j L_j
        dbar = d.mean(axis=0)
        dboot = d[idx].mean(axis=1)                       # B x m bootstrap means
        se = np.sqrt(((dboot - dbar) ** 2).mean(axis=0))
        se[se == 0] = np.inf
        t = dbar / se
        tboot = ((dboot - dbar) / se).max(axis=1)
        p = np.mean(tboot >= t.max())
        pmax = max(pmax, p)
        worst = alive[int(np.argmax(t))]
        pval[worst] = pmax
        alive.remove(worst)
    pval[alive[0]] = 1.0
    return pval


report = []
try:
    from arch.bootstrap import MCS
    import arch
    have_arch = True
except ImportError:
    have_arch = False

for c in cells:
    rr = [r for r in rows if r["cell"] == c]
    with open(os.path.join(out, f"{c}_losses.csv")) as f:
        rd = list(csv.reader(f))
    models = rd[0][1:]
    L = np.array([[float(v) for v in r[1:]] for r in rd[1:]])
    n = L.shape[0]
    idx = np.fromfile(os.path.join(out, f"{c}_idx.bin"), dtype=np.int32).reshape(-1, n) - 1
    p_np = tmax_mcs(L, idx)
    p_r = np.array([float(r["p_B1000_seed7"]) for r in rr])
    exact_diff = float(np.max(np.abs(p_np - p_r)))
    arch_p = np.full((5, len(models)), np.nan)
    if have_arch:
        for s in range(5):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = MCS(L, size=0.10, reps=5000, block_size=20, method="max",
                        bootstrap="circular", seed=20251004 + s)
                m.compute()
            pv = m.pvalues["Pvalue"]
            arch_p[s] = [pv.loc[i] for i in range(len(models))]
    for j, mdl in enumerate(models):
        r = rr[j]
        report.append(dict(cell=c, table=r["table"], sample=r["sample"], metric=r["metric"], h=r["h"],
                           model=mdl, n=n, r_pipeline=float(r["p_pipeline_seed"]),
                           r_seed_min=float(r["p_seed_min"]), r_seed_max=float(r["p_seed_max"]),
                           r_seed_mean=float(r["p_seed_mean"]),
                           arch_seed_mean=float(np.nanmean(arch_p[:, j])),
                           arch_seed_min=float(np.nanmin(arch_p[:, j])),
                           arch_seed_max=float(np.nanmax(arch_p[:, j])),
                           numpy_B1000=float(p_np[j]), r_B1000=float(p_r[j]),
                           exact_abs_diff=exact_diff))

with open(os.path.join(out, "mcs_validation.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(report[0]))
    w.writeheader(); w.writerows(report)
print("arch", arch.__version__ if have_arch else "not installed")
print("max |numpy - R| with identical resamples:", max(r["exact_abs_diff"] for r in report))
d = [abs(r["r_seed_mean"] - r["arch_seed_mean"]) for r in report]
print("max |mean_R - mean_arch| over 5 seeds each:", round(max(d), 4))
for r in report:
    print(f'{r["table"]:8s} {r["sample"]} {r["metric"]:5s} h={r["h"]:>2s} {r["model"]:7s} n={r["n"]} '
          f'R {r["r_seed_mean"]:.4f} [{r["r_seed_min"]:.4f},{r["r_seed_max"]:.4f}]  '
          f'arch {r["arch_seed_mean"]:.4f} [{r["arch_seed_min"]:.4f},{r["arch_seed_max"]:.4f}]')
