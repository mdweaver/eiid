"""Population gain of the full conditional law over the regression-only LCE set in 2 x C
(C = 3, 4), with the exact subset inequalities (verify_next_steps.py N1) as the full-law rows.
Oracle: the law of T given X is known at replicated composition types (Monte Carlo with many
draws per type), large precincts.  Pooled LCE (M0), one free slope k_c per column c < C.

Column c < C: m0c(x) = a0c + k_c x;  m1c(x) = (a0c + a1c - k_c) + (a2c + k_c) x.
Last column: 1 - sum.  Theta_O: all m_rc(x) in [0,1] at each type.
Theta_X(coord): adds averaged unit cell limits  E[max(0, x-(1-T_c))] <= x m1c(x) <= E[min(x,T_c)].
Theta_X(subset): adds  sum_{c in A} x m1c(x) <= E[min(x, T(A))] for every proper subset A.
Unit rates: beta_r ~ Dirichlet(phi * m_r(x)).  Target: group-1 rate for each column
(population average over types, group-1 weighted); the table reports the mean width ratio.
Run: python3 -I sharp_gain_study_2xC.py
"""
import itertools
import numpy as np
from scipy.optimize import linprog
rng = np.random.default_rng(11)

def subsets(C, coord):
    out = []
    for k in range(1, C):
        for A in itertools.combinations(range(C), k):
            if coord and not (len(A) == 1 or len(A) == C - 1): continue
            out.append(A)
    return out

def run(C, xs, m0, m1, phi, ndraw=40000):
    # m0, m1: functions x -> length-C probability vectors (linear in x)
    X = np.array(xs)
    # identified a's from the truth: for c<C, tau_c(x) = x m1c + (1-x) m0c = a0 + a1 x + a2 x^2
    a = []
    for c in range(C - 1):
        tau = np.array([x * m1(x)[c] + (1 - x) * m0(x)[c] for x in xs])
        a.append(np.linalg.lstsq(np.vstack([np.ones_like(X), X, X ** 2]).T, tau, rcond=None)[0])
    a = np.array(a)
    K = C - 1
    # m1c(x) = base1_c(x) + k_c (x - 1), m0c(x) = a0c + k_c x  (c < C);  last col via sums
    def rows_O():
        Aub, bub = [], []
        for x in xs:
            for c in range(K):
                b1 = a[c, 0] + a[c, 1] + a[c, 2] * x
                e = np.zeros(K); e[c] = 1
                # 0 <= m0c <= 1
                Aub.append(-e * x); bub.append(a[c, 0]); Aub.append(e * x); bub.append(1 - a[c, 0])
                Aub.append(-e * (x - 1)); bub.append(b1); Aub.append(e * (x - 1)); bub.append(1 - b1)
            # last column: m0C = 1 - sum_c m0c in [0,1]; m1C similarly
            s0 = sum(a[c, 0] for c in range(K)); s1 = sum(a[c, 0] + a[c, 1] + a[c, 2] * x for c in range(K))
            Aub.append(np.ones(K) * x); bub.append(1 - s0); Aub.append(-np.ones(K) * x); bub.append(s0)
            Aub.append(np.ones(K) * (x - 1)); bub.append(1 - s1); Aub.append(-np.ones(K) * (x - 1)); bub.append(s1)
        return Aub, bub
    def rows_X(coord):
        Aub, bub = rows_O()
        for x in xs:
            T = Tdraw[x]
            for A in subsets(C, coord):
                rhs = np.mean(np.minimum(x, T[:, list(A)].sum(1)))
                # group-1 mass sum_{c in A} x m1c(x); columns < C linear in k, column C = x - sum_{c<C}
                coef = np.zeros(K); const = 0.0
                for c in A:
                    if c < K:
                        coef[c] += x * (x - 1); const += x * (a[c, 0] + a[c, 1] + a[c, 2] * x)
                    else:
                        coef -= x * (x - 1); const += x * (1 - sum(a[cc, 0] + a[cc, 1] + a[cc, 2] * x for cc in range(K)))
                Aub.append(coef); bub.append(rhs - const)
        return Aub, bub
    Tdraw = {}
    for x in xs:   # one set of draws per type, shared by the coordinate and subset versions
        b1 = rng.dirichlet(phi * m1(x) + 1e-9, ndraw); b0 = rng.dirichlet(phi * m0(x) + 1e-9, ndraw)
        Tdraw[x] = x * b1 + (1 - x) * b0
    w = X / X.sum()
    res = {}
    for name, (Aub, bub) in {"O": rows_O(), "coord": rows_X(True), "subset": rows_X(False)}.items():
        widths = []
        for c in range(K):
            g = np.zeros(K); g[c] = np.sum(w * (X - 1))   # B1c(k) = const + k_c * sum w (x-1)
            lo = linprog(g, A_ub=np.array(Aub), b_ub=bub, bounds=[(None, None)] * K, method="highs")
            hi = linprog(-g, A_ub=np.array(Aub), b_ub=bub, bounds=[(None, None)] * K, method="highs")
            if lo.status != 0 or hi.status != 0: widths.append(np.nan); continue
            widths.append(-hi.fun - lo.fun)
        res[name] = np.array(widths)
    return res

def lin(p0, p1):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    return lambda x: (1 - x) * p0 + x * p1   # valid probability vectors for x in [0,1]

cases = {
    3: [("polarised", lin([.7, .2, .1], [.6, .2, .2]), lin([.1, .3, .6], [.05, .25, .7])),
        ("mild", lin([.4, .35, .25], [.35, .35, .3]), lin([.3, .35, .35], [.25, .35, .4]))],
    4: [("polarised", lin([.6, .2, .1, .1], [.5, .2, .2, .1]), lin([.05, .15, .3, .5], [.05, .1, .3, .55])),
        ("mild", lin([.3, .3, .2, .2], [.25, .3, .25, .2]), lin([.2, .25, .25, .3], [.2, .2, .3, .3]))],
}
designs = {"wide [.05,.95]": np.linspace(.05, .95, 8), "middle [.2,.6]": np.linspace(.2, .6, 7)}
print("mean (over columns) width ratio to Theta_O for the group-1 rates; oracle law, large precincts")
print(f"{'C':>2s} {'truth':10s} {'design':16s} {'phi':>4s} {'coord/O':>8s} {'subset/O':>9s}")
allr = []
for C, lst in cases.items():
    for tn, m1, m0 in lst:
        for dn, xs in designs.items():
            for phi in [3, 10, 50]:
                r = run(C, xs, m0, m1, phi)
                rc = np.nanmean(r["coord"] / r["O"]); rs = np.nanmean(r["subset"] / r["O"])
                allr.append((C, rc, rs))
                print(f"{C:2d} {tn:10s} {dn:16s} {phi:4d} {rc:8.2f} {rs:9.2f}")
allr = np.array(allr)
for C in (3, 4):
    s = allr[allr[:, 0] == C]
    print(f"C={C}: median subset/O {np.median(s[:,2]):.2f} (min {s[:,2].min():.2f}); "
          f"median coord/O {np.median(s[:,1]):.2f}; max extra gain of subsets over coordinates "
          f"{(s[:,1]-s[:,2]).max():.3f}")
