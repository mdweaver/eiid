# Rank and conditioning of V-hat used by the projection interval (P3.5(b2)(ii)) in each design of
# phase3_montecarlo.py. Same designs and seeds as the 1000-replication run. Only the regression and V-hat
# are recomputed, not the intervals.
# For each replication: dim(beta), numerical rank q (eigenvalues > 1e-10 * largest, the rule used by the
# script and by the package), the count of eigenvalues below that threshold, and the ratio of the smallest
# to the largest eigenvalue. Design B fits all 3 columns (dim 18), whose structural rank is 12 by adding up;
# it is also reported on the C - 1 = 2 free columns (dim 12), as P3 Remark P3.4.1 prescribes.
# Run: python3 -I phase3_montecarlo_rank.py [reps]
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phase3_montecarlo as P


def stats(V):
    ev = np.linalg.eigvalsh((V + V.T) / 2)
    thr = 1e-10 * ev.max()
    return len(V), int(np.sum(ev > thr)), int(np.sum(ev <= thr)), float(max(ev.min(), 0) / ev.max())


def one(args):
    name, seed = args
    kind, base, spec = P.DESIGNS[name]
    rng = np.random.default_rng(seed)
    if kind == "B":
        X = rng.dirichlet([4, 3, 2.5], 1000)
        M = P.M3(X)
        bu = np.stack([np.stack([rng.dirichlet(40.0 * M[i, r]) for r in range(3)]) for i in range(1000)])
        T = np.einsum("ir,irc->ic", X, bu)
        Z = np.column_stack([X[:, r] * X[:, s] for r, s in P.MONO])
        coef = np.linalg.lstsq(Z, T, rcond=None)[0]
        E = T - Z @ coef
        return {"M0 (3 columns)": stats(P.hc3_cov(Z, E)), "M0 (C-1 columns)": stats(P.hc3_cov(Z, E[:, :2]))}
    simname, simargs, models, _ = spec
    s = P.SIMS[simname](rng, *simargs)
    out = {}
    for model in models:
        dfit, Dfit = (s["d"], s["D"]) if model == "M1" else (np.zeros_like(s["d"]), 1)
        Z = P.design22(dfit, s["x"], Dfit)
        beta = np.linalg.lstsq(Z, s["t"], rcond=None)[0]
        out[model] = stats(P.hc3_cov(Z, (s["t"] - Z @ beta)[:, None]))
    return out


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    print(f"Rank of V-hat, {reps} replications per design (seeds as in phase3_montecarlo.py).")
    print("design / fit: dim(beta); q (rank used): distribution; reps with q < dim; eigenvalues <= 1e-10*max; "
          "min/max eigenvalue ratio (median, minimum)")
    with Pool(os.cpu_count() or 1) as pool:
        for name, (kind, base, spec) in P.DESIGNS.items():
            rows = pool.map(one, [(name, base + i) for i in range(reps)], chunksize=4)
            for fit in rows[0]:
                st = np.array([r[fit] for r in rows])
                qs, cnt = np.unique(st[:, 1], return_counts=True)
                print(f"  {name:13s} {fit:16s} dim {int(st[0, 0]):3d}; q: "
                      + ", ".join(f"{int(a)} x{c}" for a, c in zip(qs, cnt))
                      + f"; q < dim in {int(np.sum(st[:, 1] < st[:, 0]))}; near-zero eigenvalues {int(st[:, 2].min())}-{int(st[:, 2].max())}"
                      + f"; ratio median {np.median(st[:, 3]):.1e}, min {st[:, 3].min():.1e}")
