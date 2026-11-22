"""
Independent checks for the proof review of PHASE2e_did_general.md (2e).

Written by the reviewer from the statements in 2e, without reusing verify_phase2e.py's helpers.
Each check prints PASS/FAIL against what the REVIEW expects. A check tagged [CX] is a counterexample to a
statement of 2e as written (PASS = the counterexample was found / behaves as the review says).

Run:  python3 phase2e_review_checks.py  > phase2e_review_checks_output.txt
"""
import itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.stats import norm, chi2

rng = np.random.default_rng(20261122)
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(("PASS  " if ok else "FAIL  ") + name + (("\n      " + detail) if detail else ""), flush=True)


# ------------------------------------------------------------------------------------------
# A generic LP over the T-election product of R x C transportation polytopes, plus extra variables.
# ------------------------------------------------------------------------------------------
class Prod:
    def __init__(self, a, b):
        self.a = np.asarray(a, float)            # p x R (fixed population)
        self.b = np.asarray(b, float)            # T x p x C
        self.T, self.p, self.C = self.b.shape
        self.R = self.a.shape[1]
        self.nn = self.T * self.p * self.R * self.C
        self.nx = self.nn
        self.lo = [0.0] * self.nn
        self.hi = [None] * self.nn
        self.eq, self.ub = [], []
        for t in range(self.T):
            for i in range(self.p):
                for r in range(self.R):
                    v = {self.ix(t, i, r, c): 1.0 for c in range(self.C)}
                    self.eq.append((v, self.a[i, r]))
                for c in range(self.C):
                    v = {self.ix(t, i, r, c): 1.0 for r in range(self.R)}
                    self.eq.append((v, self.b[t, i, c]))

    def ix(self, t, i, r, c):
        return ((t * self.p + i) * self.R + r) * self.C + c

    def var(self, lo=None, hi=None):
        self.lo.append(lo); self.hi.append(hi); self.nx += 1
        return self.nx - 1

    def add_ub(self, v, rhs):  # sum v x <= rhs ; v dict
        self.ub.append((dict(v), rhs))

    def add_ge(self, v, rhs):
        self.ub.append(({k: -x for k, x in v.items()}, -rhs))

    def add_eq(self, v, rhs):
        self.eq.append((dict(v), rhs))

    def _mat(self, rows):
        A = np.zeros((len(rows), self.nx)); bb = np.zeros(len(rows))
        for j, (v, r) in enumerate(rows):
            for k, x in v.items():
                A[j, k] += x
            bb[j] = r
        return A, bb

    def opt(self, obj, sense=1):
        c = np.zeros(self.nx)
        for k, x in obj.items():
            c[k] += sense * x
        Ae, be = self._mat(self.eq)
        kw = {}
        if self.ub:
            Au, bu = self._mat(self.ub); kw = dict(A_ub=Au, b_ub=bu)
        res = linprog(c, A_eq=Ae, b_eq=be, bounds=list(zip(self.lo, self.hi)), method="highs", **kw)
        if res.status != 0:
            return None
        return sense * res.fun, res.x

    def range(self, obj, const=0.0):
        lo = self.opt(obj, 1); hi = self.opt(obj, -1)
        if lo is None or hi is None:
            return None
        return (lo[0] + const, hi[0] + const)

    def copy(self):
        import copy
        return copy.deepcopy(self)

    # linear forms -------------------------------------------------------------------
    def rate(self, r, c, t):     # B_{rct}
        A = self.a[:, r].sum()
        return {self.ix(t, i, r, c): 1.0 / A for i in range(self.p)}

    def rest(self, r, c, t):     # rate of everyone outside row r
        A = self.a.sum() - self.a[:, r].sum()
        return {self.ix(t, i, s, c): 1.0 / A for i in range(self.p) for s in range(self.R) if s != r}


def lin(*terms):
    """lin((coef, dict), ...) -> dict"""
    out = {}
    for co, d in terms:
        for k, x in d.items():
            out[k] = out.get(k, 0.0) + co * x
    return out


def gap(P, c, t, r=0, s=1):
    return lin((1, P.rate(r, c, t)), (-1, P.rate(s, c, t) if s != "rest" else P.rest(r, c, t)))


def did(P, c, u, t, r=0, s=1):
    return lin((1, gap(P, c, t, r, s)), (-1, gap(P, c, u, r, s)))


def val(d, x):
    return sum(co * x[k] for k, co in d.items())


# ------------------------------------------------------------------------------------------
# Causal LP on top (Theorem E16): extra variables theta (treated untreated rate at post) and S.
# b = D(n) + S - (B_r,post(n) - theta).
# ------------------------------------------------------------------------------------------
def causal(P, base, post, Lam=0.0, limits=True, r=0, s=1, c=0):
    th = P.var(0.0 if limits else None, 1.0 if limits else None)
    S = P.var(-Lam, Lam)
    Bs = P.rate(s, c, post) if s != "rest" else P.rest(r, c, post)
    if limits:   # B_s - S in [0,1]
        P.add_ge(lin((1, Bs), (-1, {S: 1})), 0.0)
        P.add_ub(lin((1, Bs), (-1, {S: 1})), 1.0)
    att = lin((1, P.rate(r, c, post)), (-1, {th: 1}))
    bgap = lin((1, did(P, c, base, post, r, s)), (1, {S: 1}), (-1, att))
    return th, S, att, bgap


# ==========================================================================================
# 1. Worked example of 2e Sec. 9.2, recomputed from scratch.
# ==========================================================================================
def worked():
    N = 200
    vets = [10, 50, 150, 190]
    yes = [[60, 62, 74], [62, 66, 86], [61, 63, 113], [56, 61, 120]]
    a = np.array([[v, N - v] for v in vets])
    b = np.array([[[yes[i][t], N - yes[i][t]] for i in range(4)] for t in range(3)])
    return a, b


def cr_rows(P, u, t, c=0, r=0):
    """CR on changes for (r, c, u, t), all units, R = 2: D between NM_D and ER_D (Prop. E5(d)), as LP branches.
    Returns list of branches, each a list of (dict, lo, hi) constraints on D."""
    a = P.a; N_i = a.sum(1); p = N_i / N_i.sum(); X = a[:, r] / N_i
    dT = (P.b[t, :, c] - P.b[u, :, c]) / N_i
    mu = (p * X).sum(); var = (p * (X - mu) ** 2).sum(); cov = (p * (X - mu) * (dT - (p * dT).sum())).sum()
    ER, NM = cov / var, cov / (mu * (1 - mu))
    D = did(P, c, u, t, r, 1)
    return [(D, NM, ER)] if ER >= NM else [(D, ER, NM)], ER, NM


def with_rows(P, rows):
    Q = P.copy()
    for (d, lo, hi) in rows:
        Q.add_ge(d, lo); Q.add_ub(d, hi)
    return Q


def union(ranges):
    rr = [x for x in ranges if x is not None]
    return None if not rr else (min(x[0] for x in rr), max(x[1] for x in rr))


def rm_set(P0, Mbar, base, post, pre_pairs, Lam=0.0, limits=True, rows=()):
    """Relative magnitudes, one post period: |b| <= Mbar * sigma * (g_s - g_{s-1}) over branches."""
    out = []
    for (s0, s1) in pre_pairs:
        for sg in (1, -1):
            P = with_rows(P0, rows)
            th, S, att, bg = causal(P, base, post, Lam, limits)
            A = lin((sg, gap(P, 0, s1)), (-sg, gap(P, 0, s0)))
            P.add_ub(lin((1, bg), (-Mbar, A)), 0.0)
            P.add_ub(lin((-1, bg), (-Mbar, A)), 0.0)
            out.append(P.range(att))
    return union(out)


a, b = worked()
P = Prod(a, b)
Dpost = P.range(did(P, 0, 1, 2)); Dpre = P.range(did(P, 0, 0, 1))
Pn = P.copy(); _, _, att, _ = causal(Pn, 1, 2); noinfo = Pn.range(att)
br, ER, NM = cr_rows(P, 1, 2)
Pcr = with_rows(P, br)
Dcr = Pcr.range(did(Pcr, 0, 1, 2))
# PT, no spillover
Q = Pcr.copy(); th, S, att, bg = causal(Q, 1, 2); Q.add_eq(bg, 0.0); pt = Q.range(att)
Q = Pcr.copy(); th, S, att, bg = causal(Q, 1, 2); Q.add_ub(bg, .05); Q.add_ge(bg, -.05); dial = Q.range(att)
Q = Pcr.copy(); th, S, att, bg = causal(Q, 1, 2, Lam=.02); Q.add_ub(bg, .05); Q.add_ge(bg, -.05); dialS = Q.range(att)
rm0 = rm_set(P, 1.0, 1, 2, [(0, 1)])
brpre, _, _ = cr_rows(P, 0, 1)
Dpre_cr = with_rows(P, brpre).range(did(P, 0, 0, 1))
rm_crpre_only = rm_set(P, 1.0, 1, 2, [(0, 1)], rows=brpre)
rm_cr_both = rm_set(P, 1.0, 1, 2, [(0, 1)], rows=brpre + br)


def mtr_nd_ps(P0, rows=(), Lam=0.0):
    Q = with_rows(P0, rows); th, S, att, bg = causal(Q, 1, 2, Lam)
    Q.add_ub(lin((1, {th: 1}), (-1, Q.rate(0, 0, 2))), 0)                  # MTR
    Q.add_ge(lin((1, {th: 1}), (-1, Q.rate(0, 0, 1))), 0)                  # ND
    Q.add_ge(lin((1, {th: 1}), (-1, Q.rate(1, 0, 2)), (1, {S: 1})), 0)     # PS (untreated comparison rate)
    return Q.range(att)


mnp = mtr_nd_ps(P); mnp_cr = mtr_nd_ps(P, br)
r3 = lambda x: None if x is None else tuple(round(v, 3) for v in x)
printed = {"DiD": (-0.407, 0.792), "effect none": (-0.568, 0.732), "CR DiD": (0.143, 0.270), "PT": (0.143, 0.270),
           "dial .05": (0.093, 0.320), "dial+spill": (0.073, 0.340), "RM margins": (-0.568, 0.732),
           "RM + CR (pre and post)": (0.135, 0.278), "MTR/ND/PS": (0.0, 0.482), "MTR/ND/PS + CR": (0.0, 0.311)}
got = {"DiD": r3(Dpost), "effect none": r3(noinfo), "CR DiD": r3(Dcr), "PT": r3(pt), "dial .05": r3(dial),
       "dial+spill": r3(dialS), "RM margins": r3(rm0), "RM + CR (pre and post)": r3(rm_cr_both),
       "MTR/ND/PS": r3(mnp), "MTR/ND/PS + CR": r3(mnp_cr)}
okW = all(np.allclose(got[k], printed[k], atol=1.5e-3) for k in printed)
check("W  Sec. 9.2 worked example: every printed interval reproduced by an independent LP", okW,
      "; ".join(f"{k} {got[k]}" for k in printed) +
      f"\n      pre DiD margins {r3(Dpre)} (2e prints -0.60..0.60); pre DiD under CR on the pre change {r3(Dpre_cr)} (2e: [0.004,0.008])"
      f"\n      RM with CR on the PRE change only: {r3(rm_crpre_only)} (the printed [0.135,0.278] needs CR on both changes)")

# ==========================================================================================
# 2. N21 (NS 3.2) recomputed.
# ==========================================================================================
a21 = np.array([[20, 80], [60, 40]])
Y21 = [(30, 45), (32, 48), (30, 60)]
b21 = np.array([[[Y21[t][i], 100 - Y21[t][i]] for i in range(2)] for t in range(3)])
P21 = Prod(a21, b21)
dpre, dpost = P21.range(did(P21, 0, 0, 1)), P21.range(did(P21, 0, 1, 2))
rmj = rm_set(P21, 1.0, 1, 2, [(0, 1)], limits=False)
rml = rm_set(P21, 1.0, 1, 2, [(0, 1)], limits=True)
Pq = P21.copy(); _, _, att, _ = causal(Pq, 1, 2); none21 = Pq.range(att)
Asup = max(abs(dpre[0]), abs(dpre[1]))
ew = (dpost[0] - Asup, dpost[1] + Asup)
ok = (np.allclose(np.round(dpre, 2), (-1.23, 1.27)) and np.allclose(np.round(dpost, 2), (-1.08, 1.42))
      and np.allclose(np.round(rmj, 2), (-2.35, 2.65)) and np.allclose(np.round(ew, 2), (-2.35, 2.69))
      and np.allclose(rml, none21, atol=1e-9) and np.allclose(none21, (-0.75, 1.0)))
check("N21 numbers of 2e Sec. 5.5 (and NS) reproduced; RM with rate limits = no-assumption set [-0.75, 1.00]", ok,
      f"pre {r3(dpre)}, post {r3(dpost)}, RM joint {r3(rmj)}, endpoint-wise {r3(ew)}, with limits {r3(rml)}, none {r3(none21)}")

# ==========================================================================================
# 3. [CX] Cor. E17(b): the stated sufficient condition for equality is incomplete (needs also
#    theta_PT(n+) <= 1 + M + Lam and theta_PT(n-) >= -(M + Lam)). A count-layer instance with one row.
#    Variables (R = 2, C = 2): y = treated-yes counts. Abstract version: F a segment from
#    n+ = (B_r0,B_s0,B_r1,B_s1) = (.9,0,1,.2) to n- = (.5,.5,.1,.5).
# ==========================================================================================
def e17_abstract():
    lam = np.linspace(0, 1, 100001)
    np_, nm_ = np.array([.9, 0, 1, .2]), np.array([.5, .5, .1, .5])
    pts = lam[:, None] * np_ + (1 - lam[:, None]) * nm_
    Br0, Bs0, Br1, Bs1 = pts.T
    D = (Br1 - Br0) - (Bs1 - Bs0); thPT = Br0 + Bs1 - Bs0
    feas = (thPT >= 0) & (thPT <= 1)
    return D.max(), D[feas].max(), thPT[-1], thPT[0]


Dmax, attmax, th_plus, th_minus = e17_abstract()
check("E17 [CX] Cor. E17(b) (M = Lam = 0): stated conditions hold at n+ and n-, yet max ATT < max D",
      th_plus >= 0 and th_minus <= 1 and attmax < Dmax - 1e-3,
      f"theta_PT(n+) = {th_plus:.2f} (>= 0 as required, but > 1), theta_PT(n-) = {th_minus:.2f} (<= 1); "
      f"max D = {Dmax:.3f}, max ATT under PT = {attmax:.3f}")


# realised with counts: one unit, a = (100, 100), 100 yes at both elections (each rate free in [0,1]),
# and one linear count row coupling the elections, 3 y_0 - 2 y_1 = 140 (y_t = treated-yes count).
a17 = np.array([[100, 100]]); b17 = np.array([[[100, 100]], [[100, 100]]])
P17 = Prod(a17, b17)
P17.add_eq({P17.ix(0, 0, 0, 0): 3.0, P17.ix(1, 0, 0, 0): -2.0}, 140.0)
D17 = did(P17, 0, 0, 1)
hi17, lo17 = P17.opt(D17, -1), P17.opt(D17, 1)
thPT17 = lambda x: val(P17.rate(0, 0, 0), x) + val(P17.rate(1, 0, 1), x) - val(P17.rate(1, 0, 0), x)
Q = P17.copy(); th, S, att, bg = causal(Q, 0, 1); Q.add_eq(bg, 0.0); pt17 = Q.range(att)
cond17 = thPT17(hi17[1]) >= 0 and thPT17(lo17[1]) <= 1      # Lam = M = 0: E17(b)'s conditions
check("E17 [CX] same failure with counts: one 2x2 unit, one count row coupling the elections (3 y_0 - 2 y_1 = 140)",
      cond17 and pt17 is not None and pt17[1] < hi17[0] - 1e-6,
      f"max D = {hi17[0]:.3f} at theta_PT(n+) = {thPT17(hi17[1]):.3f}; theta_PT(n-) = {thPT17(lo17[1]):.3f}; "
      f"sharp ATT set under PT {r3(pt17)}")


# margins only: n+ minimises theta_PT, so the missing condition is implied (no counterexample)
def e17_margins(trials=200):
    bad = 0
    for k in range(trials):
        p, R = 3, rng.choice([2, 3])
        N = rng.integers(10, 40, p)
        a = np.array([rng.multinomial(n, rng.dirichlet(np.ones(R))) for n in N])
        b = np.array([[rng.multinomial(n, rng.dirichlet(np.ones(2))) for n in N] for t in range(2)])
        if (a.sum(0) == 0).any():
            continue
        comp = 1 if R == 2 or rng.random() < .5 else "rest"
        P = Prod(a, b)
        Dd = did(P, 0, 0, 1, 0, comp)
        Dr = P.range(Dd)
        for M in (0.0, 0.05):
            for Lam in (0.0, 0.03):
                Q = P.copy(); th, S, att, bg = causal(Q, 0, 1, Lam, s=comp)
                Q.add_ub(bg, M); Q.add_ge(bg, -M)
                rr = Q.range(att)
                hi = P.opt(Dd, -1)[1]; lo = P.opt(Dd, 1)[1]
                Bs = P.rate(1, 0, 1) if comp != "rest" else P.rest(0, 0, 1)
                thPT = lambda x: val(P.rate(0, 0, 0), x) + val(Bs, x) - val(P.rate(1, 0, 0) if comp != "rest" else P.rest(0, 0, 0), x)
                cond = (thPT(hi) >= M + Lam and val(Bs, hi) >= Lam and thPT(lo) <= 1 - M - Lam and val(Bs, lo) <= 1 - Lam)
                if cond and (rr is None or abs(rr[1] - (Dr[1] + M + Lam)) > 1e-6 or abs(rr[0] - (Dr[0] - M - Lam)) > 1e-6):
                    bad += 1
    return bad


check("E17 [V] margins only (R = 2, named pair or rest in R = 3): stated condition => equality (no failure)",
      e17_margins() == 0)

# ==========================================================================================
# 4. [CX] Prop. E6(b): emptiness criterion omits the case Dbar^+ < 0 <= cap. Pure units.
# ==========================================================================================
# bin L: unit A all group 1 (20), unit B all group 2 (80); bin H: unit C all group 1 (80), unit D all group 2 (20)
a6 = np.array([[20, 0], [0, 80], [80, 0], [0, 20]])
yes_u = [10, 40, 40, 4]; yes_t = [8, 40, 72, 16]
b6 = np.array([[[y, n - y] for y, n in zip(yy, a6.sum(1))] for yy in (yes_u, yes_t)])
P6 = Prod(a6, b6)
L_, H_ = [0, 1], [2, 3]


def binrate(P, r, t, U):
    A = P.a[U, r].sum()
    return {P.ix(t, i, r, 0): 1.0 / A for i in U}


def bin_change_gap(P, U):
    return lin((1, binrate(P, 0, 1, U)), (-1, binrate(P, 0, 0, U)), (-1, binrate(P, 1, 1, U)), (1, binrate(P, 1, 0, U)))


Dbar = lin((.5, bin_change_gap(P6, L_)), (.5, bin_change_gap(P6, H_)))
xb = lambda U: a6[U, 0].sum() / a6[U].sum()
lam6 = xb(H_) - xb(L_)
dtau = lambda U: (b6[1, U, 0].sum() - b6[0, U, 0].sum()) / a6[U].sum()
cap = (dtau(H_) - dtau(L_)) / lam6
Drange = P6.range(Dbar)
Q = P6.copy(); Q.add_ge(Dbar, 0.0); Q.add_ub(Dbar, cap)
check("E6 [CX] Prop. E6(b): Delta tau_H >= Delta tau_L and cap >= Dbar^-, yet the set is empty (Dbar^+ < 0)",
      dtau(H_) >= dtau(L_) and cap >= Drange[0] and Q.range(Dbar) is None and Drange[1] < 0,
      f"Dbar range {r3(Drange)}, cap = {cap:.3f}, lambda = {lam6:.2f}; LP with both beliefs infeasible")

# ==========================================================================================
# 5. Prop. E3(c): corner infeasibility does not need the other elections to have ranges of positive length.
# ==========================================================================================
# T = 3, elections -1, 0, 1 = codes 0, 1, 2; elections -1 and 1 from pure units only (gap degenerate)
a3 = np.array([[30, 0], [0, 30], [15, 15]])
b3 = np.array([[[10, 20], [20, 10], [15, 15]],      # election -1: mixed unit forced? make it degenerate via b=0
               [[10, 20], [20, 10], [8, 22]],
               [[12, 18], [21, 9], [15, 15]]])
b3[0, 2] = [0, 30]; b3[2, 2] = [30, 0]            # mixed unit all-no at -1, all-yes at 1: rates identified
P3 = Prod(a3, b3)
G = [P3.range(gap(P3, 0, t)) for t in range(3)]
Dp_, Dq_ = did(P3, 0, 1, 0), did(P3, 0, 1, 2)       # D_{-1} = g_{-1} - g_0 ; D_1 = g_1 - g_0
Rp, Rq = P3.range(Dp_), P3.range(Dq_)
Q = P3.copy(); Q.add_ge(Dq_, Rq[1] - 1e-9); Q.add_ub(Dp_, Rp[0] + 1e-9)
check("E3 [W] Prop. E3(c): corner (D_1 max, D_-1 min) infeasible although only the base gap has positive length",
      abs(G[0][1] - G[0][0]) < 1e-9 and abs(G[2][1] - G[2][0]) < 1e-9 and G[1][1] - G[1][0] > 0 and Q.range(Dq_) is None,
      f"gap widths {[round(g[1] - g[0], 3) for g in G]}")


# ==========================================================================================
# 6. Brute force over integer two-election 2x2 tables: E5(a),(b),(d) with E4 identities from 2c Def. 1,
#    and E26 (lattice; CR integer emptiness).
# ==========================================================================================
def assoc(a, n1, b, u, t):
    """realised one-vs-rest associations on changes, row 0, from 2c Def. 1 (R = 2), and D."""
    N = a.sum(1); p = N / N.sum(); X = a[:, 0] / N
    beta = np.zeros((2, len(N), 2))         # t, i, r
    for k, tt in enumerate((u, t)):
        beta[k, :, 0] = n1[k] / a[:, 0]
        beta[k, :, 1] = (b[tt, :, 0] - n1[k]) / a[:, 1]
    dB = beta[1] - beta[0]
    mu = (p * X).sum()
    dT = X * dB[:, 0] + (1 - X) * dB[:, 1]
    deltaB = (p * X * (1 - X) * (dB[:, 0] - dB[:, 1])).sum()
    w1 = p * X / mu; w0 = p * (1 - X) / (1 - mu)
    cov = lambda w, f, g: (w * f * g).sum() - (w * f).sum() * (w * g).sum()
    deltaW = mu * cov(w1, dB[:, 0], X) + (1 - mu) * cov(w0, dB[:, 1], X)
    A1, A0 = a[:, 0].sum(), a[:, 1].sum()
    D = (n1[1].sum() - n1[0].sum()) / A1 - ((b[t, :, 0] - n1[1]).sum() - (b[u, :, 0] - n1[0]).sum()) / A0
    var = (p * (X - mu) ** 2).sum(); cv = (p * (X - mu) * (dT - (p * dT).sum())).sum()
    return deltaB, deltaW, D, cv / var, cv / (mu * (1 - mu)), var


ok6, nd, cr_int_empty, lattice_ok = True, 0, 0, True
for rep in range(40):
    p = 3
    N = rng.integers(3, 7, p)
    a = np.stack([rng.integers(1, N), N - rng.integers(1, N)], 1); a[:, 1] = N - a[:, 0]
    X_ = a[:, 0] / N; p_ = N / N.sum()
    if (a[:, 1] <= 0).any() or (p_ * (X_ - (p_ * X_).sum()) ** 2).sum() < 1e-9:
        continue                                   # 2e assumes Var_p(X_r) > 0
    b = np.zeros((2, p, 2), int)
    for t in range(2):
        b[t, :, 0] = rng.integers(0, N + 1); b[t, :, 1] = N - b[t, :, 0]
    rngs = [[range(max(0, a[i, 0] + b[t, i, 0] - N[i]), min(a[i, 0], b[t, i, 0]) + 1) for i in range(p)] for t in range(2)]
    Ds, crD, wD = [], [], []
    for y0 in itertools.product(*rngs[0]):
        for y1 in itertools.product(*rngs[1]):
            n1 = np.array([y0, y1], float)
            dB, dW, D, ER, NM, var = assoc(a, n1, b, 0, 1)
            ok6 &= abs((ER - D) - dW / var) < 1e-10           # E4(c) first identity
            mu = a[:, 0].sum() / a.sum()
            ok6 &= abs((NM - D) + dB / (mu * (1 - mu))) < 1e-10  # E4(c) second identity
            Ds.append(D)
            if dW * dB >= -1e-12:
                crD.append(D)
            if dW >= -1e-12:
                wD.append(D)
    nd += 1
    Ds = np.unique(np.round(Ds, 12))
    A1, A0 = a[:, 0].sum(), a[:, 1].sum(); kap = 1 / A1 + 1 / A0
    steps = np.diff(Ds) / kap
    lattice_ok &= np.allclose(steps, np.round(steps)) and np.allclose(steps, 1)   # every lattice point attained
    # continuous ends from E5(a)
    Zlo = sum(r.start for r in rngs[1]) - sum(r.stop - 1 for r in rngs[0])
    Zhi = sum(r.stop - 1 for r in rngs[1]) - sum(r.start for r in rngs[0])
    Dz = lambda Z: kap * Z - (b[1, :, 0].sum() - b[0, :, 0].sum()) / A0
    lattice_ok &= abs(Ds[0] - Dz(Zlo)) < 1e-10 and abs(Ds[-1] - Dz(Zhi)) < 1e-10
    # CR integer set = lattice points between NM and ER
    _, _, _, ER, NM, _ = assoc(a, np.array([list(r.start for r in rngs[0]), list(r.start for r in rngs[1])], float), b, 0, 1)
    lo_, hi_ = min(ER, NM), max(ER, NM)
    pred = [d for d in Ds if lo_ - 1e-10 <= d <= hi_ + 1e-10]
    ok6 &= np.allclose(sorted(set(np.round(crD, 9))), np.round(pred, 9))
    if len(crD) == 0:
        cr_int_empty += 1
check("E4/E5/E26 [V] brute force, 40 integer two-election 2x2 designs: E4(c) identities from 2c Def. 1; lattice of step kappa "
      "with DD ends; integer CR set = lattice points between NM and ER", ok6 and lattice_ok,
      f"{nd} designs; integer CR empty in {cr_int_empty}")


# ==========================================================================================
# 7. Prop. E18 with L = 2 pre-periods (2e checks L = 1 only): union of 2L = 4 LPs = grid brute force,
#    with and without rate limits.  One mixed unit + pure units: one free count per election.
# ==========================================================================================
a7 = np.array([[40, 0], [0, 60], [25, 25]])
yes7 = [[20, 30, 20], [22, 33, 15], [21, 30, 30], [30, 32, 35]]   # elections -2,-1,0,1 (codes 0..3)
b7 = np.array([[[y, n - y] for y, n in zip(yy, a7.sum(1))] for yy in yes7])
P7 = Prod(a7, b7)
rm7 = rm_set(P7, 1.0, 2, 3, [(0, 1), (1, 2)], limits=False)
rm7l = rm_set(P7, 1.0, 2, 3, [(0, 1), (1, 2)], limits=True)
# grid: free counts y_t in the mixed unit
lo7 = [max(0, 25 + b7[t, 2, 0] - 50) for t in range(4)]; hi7 = [min(25, b7[t, 2, 0]) for t in range(4)]
grids = np.meshgrid(*[np.linspace(lo7[t], hi7[t], 41) for t in range(4)], indexing="ij")
A1, A0 = 65, 85
Br = [(b7[t, 0, 0] + grids[t]) / A1 for t in range(4)]
Bs = [(b7[t, 1, 0] + b7[t, 2, 0] - grids[t]) / A0 for t in range(4)]
g = [Br[t] - Bs[t] for t in range(4)]
Arm = np.maximum(np.abs(g[1] - g[0]), np.abs(g[2] - g[1]))
D1 = g[3] - g[2]
grid_nl = ((D1 - Arm).min(), (D1 + Arm).max())
thPT = Br[2] + Bs[3] - Bs[2]
tlo = np.maximum(0, thPT - Arm); thi = np.minimum(1, thPT + Arm)
feas = tlo <= thi
grid_l = ((Br[3] - thi)[feas].min(), (Br[3] - tlo)[feas].max())
check("E18 [V] L = 2: union of 4 branch LPs = grid brute force (consecutive-difference RM), with and without rate limits",
      np.allclose(rm7, grid_nl, atol=2e-3) and np.allclose(rm7l, grid_l, atol=2e-3),
      f"LP {r3(rm7)} / grid {r3(grid_nl)}; with limits LP {r3(rm7l)} / grid {r3(grid_l)}")

# ==========================================================================================
# 8. Prop. E21(c) with a large spillover allowance (B^L_s < Lam): closed form still = LP (max{.,.} absorbs).
# ==========================================================================================
ok8, n8 = True, 0
for rep in range(60):
    p, R = 4, 3
    N = rng.integers(10, 40, p)
    a = np.array([rng.multinomial(n, rng.dirichlet(np.ones(R))) for n in N])
    b = np.array([[rng.multinomial(n, rng.dirichlet(np.ones(3))) for n in N] for t in range(2)])
    if (a.sum(0) == 0).any():
        continue
    Lam = 0.3
    for comp in (1, "rest"):
        Q = Prod(a, b); th, S, att, bg = causal(Q, 0, 1, Lam, s=comp)
        Bs1 = Q.rate(1, 0, 1) if comp != "rest" else Q.rest(0, 0, 1)
        Q.add_ub(lin((1, {th: 1}), (-1, Q.rate(0, 0, 1))), 0)
        Q.add_ge(lin((1, {th: 1}), (-1, Q.rate(0, 0, 0))), 0)
        Q.add_ge(lin((1, {th: 1}), (-1, Bs1), (1, {S: 1})), 0)
        rr = Q.range(att)
        P0 = Prod(a, b)
        BU = P0.range(P0.rate(0, 0, 1))[1]; BL0 = P0.range(P0.rate(0, 0, 0))[0]
        if comp == 1:
            Bl = P0.range(P0.rate(1, 0, 1))[0]
        else:
            Pq = P0.copy(); Pq.add_ge(Pq.rate(0, 0, 1), BU - 1e-9); Bl = Pq.range(Pq.rest(0, 0, 1))[0]
        U = BU - max(BL0, Bl - Lam)
        n8 += 1
        ok8 &= (rr is None and U < 0) or (rr is not None and abs(rr[1] - U) < 1e-6 and abs(rr[0]) < 1e-7)
check("E21 [V] Prop. E21(c), R = 3, Lam = 0.3 (spillover can exceed B^L_s): closed form = LP", ok8, f"{n8} cases")

# ==========================================================================================
# 9. Theorem E27: TU of the change-row matrix by Ghouila-Houri (all row subsets), and what lies outside it.
#    Block rows and the other columns are unit rows/columns after the substitution, so TU of the reduced
#    matrix is TU of the change rows on column c alone (review remark).
# ==========================================================================================
def gh_tu(M):
    """Ghouila-Houri over all row subsets (exact for small M)."""
    m = M.shape[0]
    for k in range(1, m + 1):
        for S in itertools.combinations(range(m), k):
            sub = M[list(S)]
            found = False
            for signs in itertools.product((1, -1), repeat=k - 1):
                s = np.array((1,) + signs)
                cs = s @ sub
                if np.all(np.abs(cs) <= 1):
                    found = True; break
            if not found:
                return False, S
    return True, None


def change_matrix(p, T, pairs_sets):
    cols = {(i, t): j for j, (i, t) in enumerate(itertools.product(range(p), range(T)))}
    rows = []
    for (u, t), fam in pairs_sets:
        for U in fam:
            v = np.zeros(len(cols))
            for i in U:
                v[cols[(i, t)]] += 1; v[cols[(i, u)]] -= 1
            rows.append(v)
    return np.array(rows)


fam_a = [[0, 1, 2], [0, 1], [0]]
fam_b = [[0, 1, 2], [1, 2], [2]]
fam_c = [[0, 1, 2], [0, 2], [1]]
cases = {
    "one pair": [((0, 1), fam_a)],
    "path of 2 pairs": [((0, 1), fam_a), ((1, 2), fam_b)],
    "two pairs sharing a base": [((1, 0), fam_a), ((1, 2), fam_b)],
    "star of 3 pairs (degree 3)": [((0, 1), fam_a), ((0, 2), fam_b), ((0, 3), fam_c)],
    "cycle of 3 pairs": [((0, 1), fam_a), ((1, 2), fam_b), ((0, 2), fam_c)],
}
tu = {}
for name, ps in cases.items():
    T = 1 + max(max(u, t) for (u, t), _ in ps)
    tu[name] = gh_tu(change_matrix(3, T, ps))[0]
check("E27 [V] GH over all row subsets: change-row matrices of the theorem's classes are TU",
      tu["one pair"] and tu["path of 2 pairs"] and tu["two pairs sharing a base"],
      "; ".join(f"{k}: {'TU' if v else 'not TU'}" for k, v in tu.items()))

# ==========================================================================================
# 10. Prop. E8(b): sharp FD set inside FD4 (Jensen), 2x2 discrete design; and degenerate case equal.
# ==========================================================================================
Lf = lambda t, x: np.maximum(0, (t - 1 + x) / x)
Uf = lambda t, x: np.minimum(1, t / x)
ok10 = True
for rep in range(200):
    x = rng.uniform(.05, .95)
    Tu = rng.uniform(0, 1, 5); Tt = np.clip(Tu + rng.normal(0, .2, 5), 0, 1)
    sh = (Lf(Tt, x).mean() - Uf(Tu, x).mean(), Uf(Tt, x).mean() - Lf(Tu, x).mean())
    ou = (Lf(Tt.mean(), x) - Uf(Tu.mean(), x), Uf(Tt.mean(), x) - Lf(Tu.mean(), x))
    ok10 &= sh[0] >= ou[0] - 1e-12 and sh[1] <= ou[1] + 1e-12
check("E8 [V] Prop. E8(b),(d): per-composition sharp interval for Dm_v lies inside FD4's (L convex, U concave)", ok10)

# ==========================================================================================
# 11. Lemma E13 Monte Carlo: dependent within-unit laws across elections and groups; Var(E) <= S_D^2,
#     tail P(|E| > z S_D) <= alpha, and the extreme law attains Var = S_D^2.
# ==========================================================================================
p = 300
a11 = np.stack([rng.integers(5, 50, p), rng.integers(5, 50, p)], 1)
wr, ws = a11[:, 0] / a11[:, 0].sum(), a11[:, 1] / a11[:, 1].sum()
SD = np.sqrt(((wr + ws) ** 2).sum())
reps = 4000
# law: common shock drives all four rates with opposite signs for the two groups
m = rng.uniform(.2, .8, (p, 4))   # means of (b_r0, b_r1, b_s0, b_s1)
Es, Eext = [], []
for k in range(reps):
    z = rng.random(p) < .5
    ext = np.where(z[:, None], [0, 1, 1, 0], [1, 0, 0, 1])      # extreme law
    Eext.append((wr * (ext[:, 1] - ext[:, 0])).sum() - (ws * (ext[:, 3] - ext[:, 2])).sum())
    u = rng.random(p)
    br = np.clip(m + 0.3 * np.where(u[:, None] < .5, [-1, 1, 1, -1], [1, -1, -1, 1]), 0, 1)
    e = br - np.clip(m + 0.3 * np.array([-1, 1, 1, -1]), 0, 1) * .5 - np.clip(m + 0.3 * np.array([1, -1, -1, 1]), 0, 1) * .5
    Es.append((wr * (e[:, 1] - e[:, 0])).sum() - (ws * (e[:, 3] - e[:, 2])).sum())
Es, Eext = np.array(Es), np.array(Eext)
z = norm.ppf(.975)
check("E13 [V] Lemma E13: Var(E)/S_D^2 <= 1 under dependent laws, = 1 for the extreme law; tail <= alpha",
      Es.var() / SD ** 2 <= 1 and abs(Eext.var() / SD ** 2 - 1) < .06 and np.mean(np.abs(Es) > z * SD) <= .05,
      f"ratio {Es.var() / SD ** 2:.3f}; extreme {Eext.var() / SD ** 2:.3f}; tail {np.mean(np.abs(Es) > z * SD):.4f}")

# ==========================================================================================
# 12. Props E24/E25 Monte Carlo with CORRELATED gap estimates (2e uses independent ones), truth at the
#     RM upper end, and at the kink D_pre = 0.
# ==========================================================================================
def mc_e25(gtrue, Sig, Mbar=1.0, reps=20000, upper=False):
    L = np.linalg.cholesky(Sig)
    gh = gtrue + rng.normal(size=(reps, 3)) @ L.T
    Dp, Dq = gtrue[1] - gtrue[0], gtrue[2] - gtrue[1]
    truth = Dq + (Mbar * abs(Dp) if upper else -Mbar * abs(Dp))
    kq = np.sqrt(chi2.ppf(.95, 3))
    sd = lambda c: np.sqrt(c @ Sig @ c)
    lows = [np.array([Mbar, -1 - Mbar, 1.0]), np.array([-Mbar, -1 + Mbar, 1.0])]  # D_post - Mbar * (+-D_pre)
    ups = lows
    jlo = np.min([gh @ c - kq * sd(c) for c in lows], 0)
    jhi = np.max([gh @ c + kq * sd(c) for c in ups], 0)
    # wait: upper end is D_post + Mbar |D_pre| = max over the two signs of D_post + Mbar*sigma*D_pre
    ups = [np.array([-Mbar, -1 + Mbar, 1.0]), np.array([Mbar, -1 - Mbar, 1.0])]
    jhi = np.max([gh @ c + kq * sd(c) for c in ups], 0)
    Dph, Dqh = gh[:, 1] - gh[:, 0], gh[:, 2] - gh[:, 1]
    cpost = np.array([0, -1.0, 1])
    nlo = Dqh - Mbar * np.abs(Dph) - z * sd(cpost); nhi = Dqh + Mbar * np.abs(Dph) + z * sd(cpost)
    return np.mean((jlo <= truth) & (truth <= jhi)), np.mean((nlo <= truth) & (truth <= nhi))


Sig = np.array([[1, .8, .6], [.8, 1, .8], [.6, .8, 1.0]])
cj1, cn1 = mc_e25(np.array([0, 3., 4.]), Sig)
cj2, cn2 = mc_e25(np.array([0, 3., 4.]), Sig, upper=True)
cj3, cn3 = mc_e25(np.array([0, 0., 1.]), Sig)
check("E25 [V] joint projection region covers with correlated estimates (lower end, upper end, kink); plug-in can under-cover",
      min(cj1, cj2, cj3) >= .945 and min(cn1, cn2) < .93,
      f"joint {cj1:.3f}/{cj2:.3f}/{cj3:.3f}; plug-in {cn1:.3f}/{cn2:.3f}/{cn3:.3f}")

# union of branch intervals (Prop. E24), truth exactly on the boundary between branches (theta = 0)
reps = 20000
thh = rng.normal(size=reps)
lo_p, hi_p = np.maximum(0, thh - z), thh + z
lo_m, hi_m = thh - z, np.minimum(0, thh + z)
cov0 = np.mean(((thh + z >= 0) & (lo_p <= 0) & (0 <= hi_p)) | ((thh - z <= 0) & (lo_m <= 0) & (0 <= hi_m)))
check("E24 [V] union of branch intervals at the branch boundary covers >= .95", cov0 >= .945, f"coverage {cov0:.3f}")

print(f"\n{sum(RESULTS)}/{len(RESULTS)} PASS")
