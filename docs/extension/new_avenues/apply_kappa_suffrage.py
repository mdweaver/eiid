"""
Variance-bounded contextual effects (the "kappa" dial) on the veterans / suffrage data.
Usage: python3 apply_kappa_suffrage.py <wi csv> <ia csv>      (data not committed)

Assumption: Var_X(E[beta_r | X]) <= kappa * Var(Y) for each group r (yes column).
Under the pooled linear model (M0), E[beta_r | X] = a_r + b_r X, so this is |b_r| <= s with
s = sqrt(kappa Var(Y) / Var(X)):  a bounded-slope (Lipschitz) restriction whose scale is set by
the data.  With yes/no/other columns: rate_n slope = k_yes, rate_v slope = d + k_yes (d identified).
kappa = 0 is Goodman.  kappa_min = smallest kappa with a non-empty set.
Same data preparation and trimmed constraint points as apply_B_suffrage.py.  Plug-in only.
"""
import sys, os
import numpy as np, pandas as pd
from scipy.optimize import linprog
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apply_B_suffrage import prep_wi, prep_ia

def run(d, label):
    d = d[d.N > 0].copy(); d["x"] = d.V / d.N; d["other"] = d.N - d.yes - d.no
    x, N = d.x.values, d.N.values
    T = np.column_stack([d.yes / d.N, d.no / d.N, d.other / d.N])
    keep = np.where((N >= 100) & (x >= np.quantile(x, .025)) & (x <= np.quantile(x, .975)))[0]
    C = np.linalg.lstsq(np.column_stack([np.ones_like(x), x, x * x]), T, rcond=None)[0]
    varY, varX = T[:, 0].var(), x.var()
    E = [np.array([1., 0.]), np.array([0., 1.]), np.array([-1., -1.])]
    def lp(kappa):
        A, b = [], []
        for c in range(3):
            a_n, c1, dd = C[0, c], C[1, c], C[2, c]
            rv0 = c1 + a_n + dd * x; rn0 = a_n + 0 * x
            for i in keep:
                A.append(-(x[i] - 1) * E[c]); b.append(rv0[i])
                A.append(-x[i] * E[c]); b.append(rn0[i])
        if kappa is not None:
            s = np.sqrt(kappa * varY / varX); dd = C[2, 0]
            A += [E[0], -E[0], E[0], -E[0]]; b += [s, s, s - dd, s + dd]
        A, b = np.array(A), np.array(b)
        a_n, c1, dd = C[0, 0], C[1, 0], C[2, 0]
        wv, wn = N * x, N * (1 - x)
        rv0 = c1 + a_n + dd * x
        tv = (np.sum(wv * rv0) / wv.sum(), np.sum(wv * (x - 1)) / wv.sum())
        tn = (a_n, np.sum(wn * x) / wn.sum())
        res = {}
        for nm, (c0, sl) in (("vet", tv), ("non", tn), ("gap", (tv[0] - tn[0], tv[1] - tn[1]))):
            lo = linprog([sl, 0], A_ub=A, b_ub=b, bounds=[(None, None)] * 2, method="highs")
            hi = linprog([-sl, 0], A_ub=A, b_ub=b, bounds=[(None, None)] * 2, method="highs")
            if lo.status or hi.status: return None
            res[nm] = (c0 + lo.fun, c0 - hi.fun)
        return res
    lo_k, hi_k = 0.0, 50.0
    if lp(hi_k) is None:
        kmin = None
    else:
        for _ in range(40):
            mid = (lo_k + hi_k) / 2
            lo_k, hi_k = (lo_k, mid) if lp(mid) is not None else (mid, hi_k)
        kmin = hi_k
    out = [f"### {label}: {len(d)} units; Var(Y_yes) {varY:.4f}, Var(X) {varX:.4f}; identified slope gap d = {C[2,0]:+.3f}",
           f"kappa floor implied by d alone: d^2 Var(X) / (4 Var(Y)) = {C[2,0]**2 * varX / (4 * varY):.3f}; "
           f"kappa_min with the [0,1] rows: {('%.3f' % kmin) if kmin is not None else 'none (M0 empty)'}",
           "| kappa | veterans' yes rate | non-veterans' yes rate | gap (v - n) |", "|---|---|---|---|"]
    f = lambda t: f"[{t[0]:.3f}, {t[1]:.3f}]"
    for k in [None, 1, 0.5, 0.25, 0.1, 0.05, 0.02] + ([round(kmin * 1.25, 3), round(kmin, 3)] if kmin else []):
        r = lp(k)
        lab = "none (M0)" if k is None else str(k)
        out.append(f"| {lab} | " + (" | ".join(f(r[n]) for n in ("vet", "non", "gap")) if r else "empty | empty | empty") + " |")
    return out

if __name__ == "__main__":
    wi, ia = pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2])
    lines = ["# kappa dial (variance-bounded contextual effects) on the suffrage data, pooled linear model, plug-in", ""]
    for lab, d in (("Wisconsin 1857", prep_wi(wi, "1857")), ("Wisconsin 1865", prep_wi(wi, "1865")),
                   ("Iowa 1857", prep_ia(ia, "57"))):
        lines += run(d, lab) + [""]
    txt = "\n".join(lines); print(txt); open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
