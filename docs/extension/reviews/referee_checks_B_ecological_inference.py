from itertools import product
from pathlib import Path
import runpy
import numpy as np

reference = Path("verify_R10.py")

def union_experiment():
    r = runpy.run_path(str(reference))
    for kind in ("boundary, ER far from 0", "boundary, ER near 0", "interior"):
        x, p = r["design"](150)
        theta, ab = r["make_truth"](x, p, kind)
        parts = r["linear_parts"](x, p)
        rows = []
        for _ in range(600):
            intervals = r["one_rep"](x, p, theta, ab, 15, parts)["CR-union"]
            if intervals is None:
                continue
            intervals = sorted(intervals)
            raw = sum(hi - lo for lo, hi in intervals)
            overlap = 0
            if len(intervals) == 2:
                (lo1, hi1), (lo2, hi2) = intervals
                overlap = max(0, min(hi1, hi2) - lo2)
            rows.append((raw, raw - overlap, overlap))
        a = np.asarray(rows)
        print(kind, np.median(a[:, :2], axis=0), sum(a[:, 2] > 0))

def misspecification_experiment():
    r = runpy.run_path(str(reference))
    # All reference functions share this global generator.
    r["one_rep"].__globals__["RNG"] = np.random.default_rng(1739)
    for n in (150, 600):
        x, p = r["design"](n)
        parts = r["linear_parts"](x, p)
        for s in (0, .08, .25):
            ab = (.40, -s, .50+s, -s)
            theta = np.array([.40, .45, .50, -.05-s])
            truth = parts[3](theta)["D"]
            results = {name: [] for name in ("JKST", "CR-slab")}
            for _ in range(400):
                out = r["one_rep"](x, p, theta, ab, 15, parts)
                for name, records in results.items():
                    ci = out[name]
                    if ci is None:
                        records.append((1, 0, np.nan))
                    else:
                        lo, hi = ci
                        records.append((0, lo <= truth <= hi, hi-lo))
            for name, records in results.items():
                a = np.asarray(records)
                print(n, s, name, truth, a[:, :2].sum(axis=0),
                      np.nanmedian(a[:, 2]))

def small_tables():
    x = np.array([[.2, .7, .1], [.6, .05, .35]])
    m = np.array([[.4, .8, .1], [.42, .82, .12]])
    p = np.array([.5, .5])
    focal = x[:, 0]
    def contribution(share, rate):
        weights = p * share
        center = np.average(focal, weights=weights)
        return weights @ ((focal - center) * rate)
    rest = np.sum(x[:, 1:] * m[:, 1:], axis=1) / (1-focal)
    w1 = contribution(focal, m[:, 0]) + contribution(1-focal, rest)
    w2 = sum(contribution(x[:, j], m[:, j]) for j in range(3))
    print("Rest rates, W1, W2", rest, w1, w2)
    # Each four-voter precinct has one P vote, which can go to either row.
    gaps = {k/3 - (2-k)/5 for votes in product((0, 1), repeat=2)
            for k in [sum(votes)]}
    print("Integer-feasible gaps", sorted(gaps))

union_experiment()
misspecification_experiment()
small_tables()
