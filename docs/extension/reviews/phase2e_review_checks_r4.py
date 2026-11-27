"""PHASE2e proof review, round 4 (Prop. E35, realised within-precinct parallel trends with tolerance).
Independent of verify_phase2e.py: its own LP in the veteran yes counts x[t,i], and M* by bisection
instead of the M-as-variable LP.  Run: python3 phase2e_review_checks_r4.py"""
import numpy as np
from scipy.optimize import linprog

res = []
def check(name, ok, info=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{info}]" if info else ""))

aV = np.array([10, 50, 150, 190.]); aN = np.array([190, 150, 50, 10.])
yes = np.array([[60, 62, 61, 56], [62, 66, 63, 61], [74, 86, 113, 120.]])   # t = 0,1,2 ; base 1, post 2
wV = aV / aV.sum(); P = 4
AUD = {(1, 2): 20, (2, 2): 22, (1, 3): 4, (2, 3): 5}                        # hypothetical N yes counts

# variables: x1[i], x2[i] (veteran yes at base/post), th[i], S[i]  -> 16
def idx(k, i): return 4 * k + i
def solve(M, Lam=0.0, audit=False, extra=None, tau=None):
    lo, hi = [], []
    for k, t in ((0, 1), (1, 2)):
        for i in range(P):
            l, u = max(0, yes[t, i] - aN[i]), min(aV[i], yes[t, i])
            if audit and (t, i) in AUD:
                l = u = yes[t, i] - AUD[(t, i)]
            lo.append(l); hi.append(u)
    lo += [0] * P + [-Lam] * P; hi += [1] * P + [Lam] * P
    Aub, bub = [], []
    for i in range(P):
        c = np.sqrt(1 / aV[i] + 1 / aN[i]) * M if tau is None else tau[i]
        # gap = th - x1/aV - ((y2-x2)/aN - S - (y1-x1)/aN)
        g = np.zeros(16); g[idx(2, i)] = 1; g[idx(0, i)] = -1 / aV[i] - 1 / aN[i]; g[idx(1, i)] = 1 / aN[i]; g[idx(3, i)] = 1
        g0 = -(yes[2, i] - yes[1, i]) / aN[i]           # constant part of gap
        Aub += [g, -g]; bub += [c - g0, c + g0]
        # untreated N post rate (y2-x2)/aN - S in [0,1]
        h = np.zeros(16); h[idx(1, i)] = -1 / aN[i]; h[idx(3, i)] = -1; h0 = yes[2, i] / aN[i]
        Aub += [h, -h]; bub += [1 - h0, h0]
    att = np.zeros(16)
    for i in range(P):
        att[idx(1, i)] += wV[i] / aV[i]; att[idx(2, i)] -= wV[i]
    if extra is not None:
        Aub.append(extra[0]); bub.append(extra[1])
    out = []
    for s in (1, -1):
        r = linprog(s * att, A_ub=np.array(Aub), b_ub=np.array(bub), bounds=list(zip(lo, hi)), method="highs")
        out.append(s * r.fun if r.status == 0 else np.nan)
    return out, att

def DW_range(audit):
    # D_W = sum wV[(x2-x1)/aV - ((y2-x2)-(y1-x1))/aN]  over the count polytope
    lo, hi = [], []
    for k, t in ((0, 1), (1, 2)):
        for i in range(P):
            l, u = max(0, yes[t, i] - aN[i]), min(aV[i], yes[t, i])
            if audit and (t, i) in AUD: l = u = yes[t, i] - AUD[(t, i)]
            lo.append(l); hi.append(u)
    c = np.zeros(8); c0 = 0
    for i in range(P):
        c[4 + i] += wV[i] * (1 / aV[i] + 1 / aN[i]); c[i] -= wV[i] * (1 / aV[i] + 1 / aN[i])
        c0 -= wV[i] * (yes[2, i] - yes[1, i]) / aN[i]
    r = [s * linprog(s * c, bounds=list(zip(lo, hi)), method="highs").fun + c0 for s in (1, -1)]
    return r

R = {M: solve(M, audit=True)[0] for M in (0.0, 0.1, 0.2)}
mo = solve(0.0)[0]
dw = DW_range(True)
# no-assumption set: theta free, no rows -> tau = 1 (rows slack) and Lam big
na = solve(0, Lam=1, tau=[2] * 4)[0]
print("  margins only M=0:", np.round(mo, 4), " no-assumption:", np.round(na, 4))
print("  audited:", {k: np.round(v, 4) for k, v in R.items()}, " D_W count set:", np.round(dw, 4))
check("R4.1 worked example reproduced by an independent LP (M = 0, 0.1; margins only; D_W set)",
      np.allclose(R[0.0], [0.0525, 0.3522], atol=6e-5) and np.allclose(R[0.1], [0.0310, 0.3740], atol=6e-5)
      and np.allclose(mo, [-0.5675, 0.7322], atol=6e-5) and np.allclose(dw, [-0.0087, 0.3772], atol=6e-5)
      and np.allclose(na, [-0.5675, 0.7325], atol=6e-5))

lo_, hi_ = 0.0, 1.0
for _ in range(60):
    m = (lo_ + hi_) / 2
    (l, _u), _a = solve(m, audit=True)
    lo_, hi_ = (lo_, m) if l <= 0 else (m, hi_)
check("R4.2 breakdown M* by bisection on the lower end = 0.2438 (the note's one-LP value)", abs(hi_ - 0.2438) < 6e-5, f"{hi_:.5f}")

# nesting in M and in Lam on a grid
Ms = np.linspace(0, 1, 21); Ls = (0, 0.02, 0.05)
nest = True
for La in Ls:
    prev = None
    for M in Ms:
        cur = solve(M, Lam=La, audit=True)[0]
        if prev is not None and not (cur[0] <= prev[0] + 1e-9 and cur[1] >= prev[1] - 1e-9): nest = False
        prev = cur
check("R4.3 nested in M for each spillover bound Lam (and in Lam for each M), audited example", nest)

# LCE row D_W^S(n) <= u : a linear row in n (here all precincts supported); it can only shrink the set
_, att = solve(0, audit=False)
row = np.zeros(16)
for i in range(P):
    row[idx(1, i)] += wV[i] * (1 / aV[i] + 1 / aN[i]); row[idx(0, i)] -= wV[i] * (1 / aV[i] + 1 / aN[i])
c0 = -sum(wV[i] * (yes[2, i] - yes[1, i]) / aN[i] for i in range(P))
cut = solve(0.0, extra=(row, 0.30 - c0))[0]; cut2 = solve(0.0, extra=(-row, -(0.10 - c0)))[0]
check("R4.4 LCE row on D_W^S is linear in n; at M=0, Lam=0 the cut ATT set equals the D_W cut (ATT = D_W)",
      cut[1] <= 0.30 + 1e-9 and cut2[0] >= 0.10 - 1e-9, f"D_W<=.30: {np.round(cut,4)}, D_W>=.10: {np.round(cut2,4)}")

# Calibration of M (R4-1): K precincts, persons' untreated changes iid within precinct with a common
# expected change for V and N (so expected WPT holds exactly). Share of draws in which RWPT-tau holds
# in EVERY supported precinct, by M.
rng = np.random.default_rng(4)
def allhold(K, M, reps=2000, aVk=40, aNk=60):
    ok = 0
    for _ in range(reps):
        p0 = rng.uniform(.2, .6, K); d = rng.uniform(-.1, .1, K)
        g = np.empty(K)
        for grp, a in ((0, aVk), (1, aNk)):
            b0 = rng.binomial(a, p0) / a; b1 = rng.binomial(a, np.clip(p0 + d, 0, 1)) / a
            g = (b1 - b0) if grp == 0 else g - (b1 - b0)
        tau = M * np.sqrt(1 / aVk + 1 / aNk)
        ok += np.all(np.abs(g) <= tau)
    return ok / reps
tab = {(K, M): allhold(K, M) for K in (10, 100) for M in (1, 2, 3, 4)}
for K in (10, 100):
    print(f"  K={K}: P(RWPT-tau holds in every precinct) by M=1..4:", [float(tab[(K, M)]) for M in (1, 2, 3, 4)])
check("R4.5 (calibration, R4-1) with pure sampling noise and exact expected WPT, M = 1 fails in some precinct in most draws (K = 10) and in essentially all (K = 100); the M that holds everywhere grows with K",
      tab[(10, 1)] < 0.5 and tab[(100, 1)] < 0.05 and tab[(100, 2)] < tab[(10, 2)] and tab[(100, 3)] > 0.95,
      f"K=10: M=1 {tab[(10,1)]:.3f}, M=2 {tab[(10,2)]:.3f}; K=100: M=1 {tab[(100,1)]:.3f}, M=2 {tab[(100,2)]:.3f}, M=3 {tab[(100,3)]:.3f}")

print(f"\n{sum(res)}/{len(res)} PASS")
