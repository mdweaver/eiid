"""
Avenue C, ordered outcomes (e.g. governor vote and a Black-suffrage referendum).
2 groups (veteran share x); each man casts a governor vote (Rep/Dem) and a referendum vote
(yes/no); per county only the one-way margins x, R (Rep share), y (yes share) are observed.
Nesting ("Guttman") assumption: nobody votes Dem for governor and yes on suffrage.
Checks, over 500 random counties with y <= R:
 (a) nesting does not change the Duncan-Davis bound on veterans' yes rate;
 (b) it narrows the bound on the cross-outcome contrast (veterans' yes rate - Rep rate).
No abstention or roll-off in this check.
"""
import numpy as np
from scipy.optimize import linprog
rng = np.random.default_rng(3)
def idx(g, a, b): return g * 4 + a * 2 + b      # a: 0 Rep, 1 Dem; b: 0 yes, 1 no
def bounds(x, R, y, nest, target):
    A, bb = [], []
    for g, s in ((0, x), (1, 1 - x)):
        r = np.zeros(8); r[[idx(g, a, b) for a in (0, 1) for b in (0, 1)]] = 1; A.append(r); bb.append(s)
    r = np.zeros(8); r[[idx(g, 0, b) for g in (0, 1) for b in (0, 1)]] = 1; A.append(r); bb.append(R)
    r = np.zeros(8); r[[idx(g, a, 0) for g in (0, 1) for a in (0, 1)]] = 1; A.append(r); bb.append(y)
    ub = [(0, 0) if (nest and a == 1 and b == 0) else (0, None) for g in (0, 1) for a in (0, 1) for b in (0, 1)]
    c = np.zeros(8)
    for a in (0, 1): c[idx(0, a, 0)] += 1 / x
    if target == "contrast":
        for b in (0, 1): c[idx(0, 0, b)] -= 1 / x
    lo = linprog(c, A_eq=A, b_eq=bb, bounds=ub).fun; hi = -linprog(-c, A_eq=A, b_eq=bb, bounds=ub).fun
    return lo, hi
d1 = 0; w = []
for _ in range(500):
    x = rng.uniform(.1, .6); R = rng.uniform(.3, .8); y = rng.uniform(.05, R)
    a0, a1 = bounds(x, R, y, False, "rate"), bounds(x, R, y, True, "rate")
    d1 = max(d1, abs(a0[0] - a1[0]), abs(a0[1] - a1[1]))
    c0, c1 = bounds(x, R, y, False, "contrast"), bounds(x, R, y, True, "contrast")
    w.append((c1[1] - c1[0]) / (c0[1] - c0[0]))
txt = (f"(a) max change in DD for veterans' yes rate from nesting: {d1:.1e}\n"
       f"(b) contrast (yes - Rep) width, nesting / no nesting: median {np.median(w):.2f}, range {min(w):.2f}-{max(w):.2f}\n")
print(txt, end=""); open(__file__.replace(".py", "_output.md"), "w").write(txt)
