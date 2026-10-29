# Worked 2x2 example: riding dummies (M1) vs one province-wide model (M0).
# Rows: P = voted for party P last time (share x), O = everyone else (1-x).
# Outcome: votes P this time.  loyalty = rate for P-voters, crossover = rate for O-voters.
# Population (noise-free) calculation: polls spread uniformly over each riding's x-range.
# Run: python3 -I example_2x2_ridings.py
import numpy as np

g = 0.10                                    # common contextual slope (truth), both rows
ridings = {                                 # x-range, loyalty intercept, crossover intercept
    "Riding 1 (P stronghold)": dict(rng=(0.5, 0.9), aP=0.89, aO=0.20),
    "Riding 2 (P weak)":       dict(rng=(0.1, 0.5), aP=0.45, aO=0.02),
}

def interval_for_g(c0, c1, xs):
    """Values of the unknown slope g allowed by 0<=crossover<=1 and 0<=loyalty<=1
    at the riding's extreme polls. Identified: T = c0 + c1 x (curvature 0 here).
    Given g: crossover = c0 + g x ; loyalty = (c0 + c1 - g) + g x = c0 + c1 - g (1 - x)."""
    lo, hi = -np.inf, np.inf
    for x in xs:
        lo = max(lo, -c0 / x, (c0 + c1 - 1) / (1 - x))   # crossover>=0, loyalty<=1
        hi = min(hi, (1 - c0) / x, (c0 + c1) / (1 - x))  # crossover<=1, loyalty>=0
    return lo, hi

grid = {k: np.linspace(*v["rng"], 2001) for k, v in ridings.items()}
num = den = 0.0
lo_all, hi_all = -np.inf, np.inf
rows = []
for k, v in ridings.items():
    c0, c1 = v["aO"], v["aP"] - v["aO"] + g             # what the riding's polls reveal
    lo, hi = interval_for_g(c0, c1, v["rng"])
    lo_all, hi_all = max(lo_all, lo), min(hi_all, hi)
    x = grid[k]
    num += np.sum(x * (v["aP"] + g * x)); den += np.sum(x)
    rows.append((k, c0, c1, lo, hi))
B_true = num / den

def B_M1(gv):
    n = 0.0
    for k, v in ridings.items():
        x = grid[k]; c0, c1 = v["aO"], v["aP"] - v["aO"] + g
        n += np.sum(x * (c0 + c1 - gv * (1 - x)))
    return n / den

print("What each riding's polls reveal:  T = c0 + c1*x")
for k, c0, c1, lo, hi in rows:
    print(f"  {k}: c0={c0:.2f} c1={c1:.2f}  ->  g in [{lo:+.3f}, {hi:+.3f}]")
print(f"M1 (common g, riding intercepts): g in [{lo_all:+.3f}, {hi_all:+.3f}]  (truth {g:+.2f})")
print(f"  B (province loyalty) in [{B_M1(hi_all):.3f}, {B_M1(lo_all):.3f}]   truth {B_true:.3f}")

# M0: one model for the whole province, fitted to the pooled polls
x = np.concatenate(list(grid.values()))
t = np.concatenate([v["aO"] + (v["aP"] - v["aO"] + g) * grid[k]
                    for k, v in ridings.items()])
w0, c1, d1 = np.linalg.lstsq(np.c_[np.ones_like(x), x, x**2], t, rcond=None)[0]
lo, hi = -np.inf, np.inf
for v in (x.min(), x.max()):
    lo = max(lo, -w0 / v, (w0 + c1 + d1 * v - 1) / (1 - v))
    hi = min(hi, (1 - w0) / v, (w0 + c1 + d1 * v) / (1 - v))
B0 = lambda w1: np.sum(x * (w0 + c1 + d1 * x - w1 * (1 - x))) / np.sum(x)
print(f"M0 pooled fit: T = {w0:.3f} + {c1:.3f} x + {d1:.3f} x^2 ; crossover slope in [{lo:+.3f}, {hi:+.3f}]")
print(f"  B in [{B0(hi):.3f}, {B0(lo):.3f}]")
dd = (np.sum(np.maximum(0, t - (1 - x))) / np.sum(x), np.sum(np.minimum(t, x)) / np.sum(x))
print(f"DD (no model): B in [{dd[0]:.3f}, {dd[1]:.3f}]")
