"""Checks for PHASE2c_monotone.md (monotone ecological inference in R x C).

Run:  python3 -I verify_phase2c.py > verify_phase2c_output.txt
"""
import itertools
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(2026_11_01)
results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))


def wcov(a, b, w):
    w = w / w.sum()
    return np.sum(w * a * b) - np.sum(w * a) * np.sum(w * b)


# ---------------------------------------------------------------------------
# helpers: unit-level LP over beta (R rows, C columns) for one population
# ---------------------------------------------------------------------------
class Pop:
    """Finite population: p units, weights pw, X (p,R), T (p,C)."""

    def __init__(self, N, X, beta):
        self.N, self.X, self.beta = N, X, beta
        self.p, self.R, self.C = beta.shape
        self.pw = N / N.sum()
        self.T = np.einsum("ir,irc->ic", X, beta)
        self.mu = self.pw @ X

    def idx(self, i, r, c):
        return (i * self.R + r) * self.C + c

    def nvar(self):
        return self.p * self.R * self.C

    def equalities(self):
        A, b = [], []
        for i in range(self.p):
            for r in range(self.R):          # rows sum to 1
                a = np.zeros(self.nvar())
                for c in range(self.C):
                    a[self.idx(i, r, c)] = 1
                A.append(a); b.append(1.0)
            for c in range(self.C - 1):      # accounting identity
                a = np.zeros(self.nvar())
                for r in range(self.R):
                    a[self.idx(i, r, c)] = self.X[i, r]
                A.append(a); b.append(self.T[i, c])
        return np.array(A), np.array(b)

    # linear forms in beta (vectors of length nvar) --------------------------
    def B(self, r, c):
        a = np.zeros(self.nvar())
        for i in range(self.p):
            a[self.idx(i, r, c)] = self.pw[i] * self.X[i, r] / self.mu[r]
        return a

    def dB_pair(self, r, s, c):        # E_p[X_r X_s (beta_r - beta_s)]
        a = np.zeros(self.nvar())
        for i in range(self.p):
            f = self.pw[i] * self.X[i, r] * self.X[i, s]
            a[self.idx(i, r, c)] += f
            a[self.idx(i, s, c)] -= f
        return a

    def dB_rest(self, r, c):
        return sum(self.dB_pair(r, s, c) for s in range(self.R) if s != r)

    def dW_rest(self, r, c):
        """W1: mu_r Cov_{w^r}(beta_r, X_r) + (1-mu_r) Cov_{w^notr}(beta_notr, X_r), linear in beta."""
        pw, X, mu = self.pw, self.X, self.mu
        a = np.zeros(self.nvar())
        wr = pw * X[:, r]                        # un-normalised, sums to mu_r
        mxr = np.sum(wr * X[:, r]) / mu[r]
        for i in range(self.p):
            a[self.idx(i, r, c)] += wr[i] * (X[i, r] - mxr)
        wn = pw * (1 - X[:, r])                  # sums to 1-mu_r
        mxn = np.sum(wn * X[:, r]) / (1 - mu[r])
        for i in range(self.p):                  # (1-X_r) beta_notr = sum_{s!=r} X_s beta_s
            for s in range(self.R):
                if s != r:
                    a[self.idx(i, s, c)] += pw[i] * X[i, s] * (X[i, r] - mxn)
        return a

    def value(self, a, beta=None):
        beta = self.beta if beta is None else beta
        return a @ beta.reshape(-1)

    def lp_range(self, target, rows_ge=(), rows_le=()):
        Aeq, beq = self.equalities()
        Aub, bub = [], []
        for a, k in rows_ge:                     # a.beta >= k
            Aub.append(-a); bub.append(-k)
        for a, k in rows_le:                     # a.beta <= k
            Aub.append(a); bub.append(k)
        kw = dict(A_eq=Aeq, b_eq=beq, bounds=[(0, 1)] * self.nvar(), method="highs")
        if Aub:
            kw.update(A_ub=np.array(Aub), b_ub=np.array(bub))
        lo = linprog(target, **kw)
        hi = linprog(-target, **kw)
        if lo.status != 0 or hi.status != 0:
            return None
        return lo.fun, -hi.fun


def random_pop(p=25, R=2, C=2, conc=None):
    N = rng.integers(50, 2000, p).astype(float)
    X = rng.dirichlet(conc if conc is not None else [2.0] * R, p)
    beta = rng.dirichlet([1.0] * C, (p, R))
    return Pop(N, X, beta)


# ---------------------------------------------------------------------------
# C1. Lemma M1: R x C identities (R = 3, C = 3), one-vs-rest collapse.
# ---------------------------------------------------------------------------
worst = 0.0
for _ in range(100):
    P = random_pop(p=30, R=3, C=3)
    for c in range(3):
        t = P.T[:, c]
        M = (P.pw[:, None, None] * P.X[:, :, None] * P.X[:, None, :]).sum(0)
        g = np.linalg.solve(M, (P.pw[:, None] * P.X * t[:, None]).sum(0))
        B = np.array([P.value(P.B(r, c)) for r in range(3)])
        Delta = sum(P.mu[s] * np.array([wcov(P.X[:, r], P.beta[:, s, c], P.pw * P.X[:, s])
                                        for r in range(3)]) for s in range(3))
        worst = max(worst, np.max(np.abs(g - B - np.linalg.solve(M, Delta))), abs(Delta.sum()))
        for r in range(3):
            NM = np.sum(P.pw * P.X[:, r] * t) / P.mu[r]
            worst = max(worst, abs(NM - B[r] + P.value(P.dB_rest(r, c)) / P.mu[r]))
            # collapse to r vs rest and check EGGHM identities with the linear forms
            x = P.X[:, r]
            mu = P.mu[r]
            Bnot = sum(P.mu[s] * B[s] for s in range(3) if s != r) / (1 - mu)
            D = B[r] - Bnot
            varxn = wcov(x, x, P.pw)
            ER = wcov(t, x, P.pw) / varxn
            NMd = NM - np.sum(P.pw * (1 - x) * t) / (1 - mu)
            dW = P.value(P.dW_rest(r, c))
            dB = P.value(P.dB_rest(r, c))
            worst = max(worst, abs(ER - D - dW / varxn), abs(NMd - D + dB / (mu * (1 - mu))))
check("C1 Lemma M1 (a)-(c): Goodman/NM bias formulas, 1'Delta=0, one-vs-rest collapse (R=C=3)",
      worst < 1e-10, f"max err {worst:.1e}")


# ---------------------------------------------------------------------------
# C2. Prop M3 (2 rows): closed forms, emptiness iff ER < MOB-, NM always
#     feasible, CR never empty and single-branch, sign(D) = sign(ER) under CR.
# ---------------------------------------------------------------------------
worst, n_empty, n_emp_ok, cr_ok, sign_ok, n_cr_true = 0.0, 0, True, True, True, 0
for _ in range(300):
    P = random_pop(p=20)
    c = 0
    tgt = P.B(0, c) - P.B(1, c)
    D = P.value(tgt)
    mob = P.lp_range(tgt)
    t = P.T[:, c]; x = P.X[:, 0]; mu = P.mu[0]
    ER = wcov(t, x, P.pw) / wcov(x, x, P.pw)
    NM = np.sum(P.pw * x * t) / mu - np.sum(P.pw * (1 - x) * t) / (1 - mu)
    gam = wcov(x, x, P.pw) / (mu * (1 - mu))
    worst = max(worst, abs(NM - gam * ER))
    dW, dB = P.dW_rest(0, c), P.dB_rest(0, c)
    w_ge = P.lp_range(tgt, rows_ge=[(dW, 0)])
    if ER < mob[0] - 1e-12:
        n_empty += 1
        n_emp_ok &= w_ge is None
    else:
        worst = max(worst, abs(w_ge[0] - mob[0]), abs(w_ge[1] - min(mob[1], ER)))
    b_ge = P.lp_range(tgt, rows_ge=[(dB, 0)])
    b_le = P.lp_range(tgt, rows_le=[(dB, 0)])
    worst = max(worst, abs(b_ge[0] - NM), abs(b_ge[1] - mob[1]), abs(b_le[0] - mob[0]), abs(b_le[1] - NM))
    plus = P.lp_range(tgt, rows_ge=[(dW, 0), (dB, 0)])
    minus = P.lp_range(tgt, rows_le=[(dW, 0), (dB, 0)])
    branches = [b for b in (plus, minus) if b is not None]
    cr_ok &= len(branches) >= 1 and (len(branches) == 1 or abs(ER) < 1e-12)
    cf = (NM, min(ER, mob[1])) if ER >= 0 else (max(ER, mob[0]), NM)
    worst = max(worst, abs(branches[0][0] - cf[0]), abs(branches[0][1] - cf[1]))
    vW, vB = P.value(dW), P.value(dB)
    if vW * vB >= 0:
        n_cr_true += 1
        sign_ok &= np.sign(D) == np.sign(ER) or abs(D) < 1e-12
check("C2 Prop M3: 2x2 closed forms = LP; NM = gamma*ER; dB-sign sets never empty", worst < 1e-8,
      f"max err {worst:.1e}")
check("C2 Prop M3(b): dW>=0 set empty exactly when ER < MOB-", n_emp_ok, f"{n_empty} empty cases")
check("C2 Prop M3(d): CR never empty (continuous, no U_i), exactly one branch (ER != 0, gamma < 1)", cr_ok)
check("C2 Prop M3(d): under realised CR, sign(D) = sign(ER)", sign_ok, f"{n_cr_true} populations satisfy CR")


# ---------------------------------------------------------------------------
# C3. Prop M4: joint restrictions (several columns) are weakly tighter than
#     imposing only the target column's restriction; find a strict example.
# ---------------------------------------------------------------------------
gains, weak_ok = [], True
for _ in range(150):
    P = random_pop(p=12, R=2, C=3)
    tgt = P.B(0, 0)
    W0, W1 = P.dW_rest(0, 0), P.dW_rest(0, 1)
    single = P.lp_range(tgt, rows_ge=[(W0, 0)])
    joint = P.lp_range(tgt, rows_ge=[(W0, 0), (W1, 0)])
    if single is None or joint is None:
        continue
    weak_ok &= joint[0] >= single[0] - 1e-9 and joint[1] <= single[1] + 1e-9
    gains.append((single[1] - single[0]) - (joint[1] - joint[0]))
gains = np.array(gains)
check("C3 Prop M4: joint <= single (nested), and strictly narrower in some cases", weak_ok and gains.max() > 1e-4,
      f"{(gains > 1e-6).sum()}/{gains.size} strict; max width gain {gains.max():.3f}")


# ---------------------------------------------------------------------------
# C4. Prop M5 (K-rows): under LCE, superpopulation associations are affine in k;
#     mean-level ordering implies dB*(r,s) >= 0.
# ---------------------------------------------------------------------------
def skew_basis(R):
    return [(r, s) for r in range(R) for s in range(r + 1, R)]


def H_from(Q, k, R, C):
    """Q: (C,R,R) symmetric with sum_c Q_c = J; k: ((C-1)*npairs,) -> H (C,R,R)."""
    pairs = skew_basis(R)
    K = np.zeros((C, R, R))
    for c in range(C - 1):
        for j, (r, s) in enumerate(pairs):
            v = k[c * len(pairs) + j]
            K[c, r, s] += v; K[c, s, r] -= v
    K[C - 1] = -K[:C - 1].sum(0)
    return Q + K


def m_of(H, X):      # m[i, r, c] = e_r' H_c X_i
    return np.einsum("crs,is->irc", H, X)


worst, order_ok = 0.0, True
for _ in range(100):
    R, C, p = 3, 3, 40
    N = rng.integers(50, 2000, p).astype(float); pw = N / N.sum()
    X = rng.dirichlet([2, 2, 2], p)
    H0 = rng.dirichlet([1] * C, (R, R)).transpose(2, 0, 1)        # H0[c, r, s], sum_c = J
    Q = 0.5 * (H0 + H0.transpose(0, 2, 1))
    d0 = (C - 1) * len(skew_basis(R))
    k1, k2 = rng.normal(0, .1, d0), rng.normal(0, .1, d0)
    def assoc(k):
        m = m_of(H_from(Q, k, R, C), X)
        out = []
        for c in range(C):
            for r, s in itertools.permutations(range(R), 2):
                out.append(np.sum(pw * X[:, r] * X[:, s] * (m[:, r, c] - m[:, s, c])))
            mu = pw @ X
            out += list(sum(mu[s] * np.array([wcov(X[:, r], m[:, s, c], pw * X[:, s]) for r in range(R)])
                            for s in range(R)))
        return np.array(out)
    a = 0.37
    worst = max(worst, np.max(np.abs(assoc(a * k1 + (1 - a) * k2) - a * assoc(k1) - (1 - a) * assoc(k2))))
    # ordering (rebuilt after the Batch B review): H0[:, r, s] is row r's transition vector at vertex e_s.
    # At each vertex where row 0 has the smaller column-0 rate, swap the WHOLE transition vectors of rows 0
    # and 1. Every row stays a probability vector at every vertex, so M_H(x) is a valid transition matrix on
    # the simplex, and M_H(x)[0,0] >= M_H(x)[1,0] there by linearity. (The old code swapped column-0 entries
    # only, which breaks the row sums.)
    Hord = H0.copy()
    for s_ in range(R):
        if H0[0, 0, s_] < H0[0, 1, s_]:
            Hord[:, 0, s_], Hord[:, 1, s_] = H0[:, 1, s_].copy(), H0[:, 0, s_].copy()
    m = m_of(Hord, X)
    valid = np.all(m >= -1e-14) and np.allclose(m.sum(2), 1) and np.all(m[:, 0, 0] >= m[:, 1, 0] - 1e-14)
    order_ok &= valid and np.sum(pw * X[:, 0] * X[:, 1] * (m[:, 0, 0] - m[:, 1, 0])) >= -1e-14
check("C4 Prop M5: superpopulation associations affine in k (R=C=3)", worst < 1e-12, f"max err {worst:.1e}")
check("C4 Prop M5(b): mean-level ordering (valid transition matrices, rebuilt) => dB*(r,s) >= 0", order_ok)


# ---------------------------------------------------------------------------
# C5. N2 (R = 2): outer LP on support points.
#  (a) truth inside under a monotone nonlinear DGP;
#  (b) M0 + signed slopes (outer) is contained in N2 (outer) when tau is quadratic;
#  (c) R >= 3: monotone in a coordinate does not sign Cov(m, X_s).
# ---------------------------------------------------------------------------
def n2_bounds(xs, tau, wts, target_rc, mono, C):
    """Variables m[j, r, c] (J points, 2 rows, C cols). mono: dict (r,c) -> +1/-1."""
    J = len(xs)
    nv = J * 2 * C
    id_ = lambda j, r, c: (j * 2 + r) * C + c
    Aeq, beq, Aub, bub = [], [], [], []
    for j in range(J):
        for r in range(2):
            a = np.zeros(nv); a[[id_(j, r, c) for c in range(C)]] = 1; Aeq.append(a); beq.append(1)
        for c in range(C - 1):
            a = np.zeros(nv); a[id_(j, 0, c)] = xs[j]; a[id_(j, 1, c)] = 1 - xs[j]
            Aeq.append(a); beq.append(tau[j, c])
    order = np.argsort(xs)
    for (r, c), sgn in mono.items():
        for j0, j1 in zip(order[:-1], order[1:]):
            a = np.zeros(nv); a[id_(j1, r, c)] = -sgn; a[id_(j0, r, c)] = sgn
            Aub.append(a); bub.append(0)
    r, c = target_rc
    obj = np.zeros(nv)
    for j in range(J):
        obj[id_(j, r, c)] = wts[j]
    kw = dict(A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, 1)] * nv, method="highs")
    if Aub:
        kw.update(A_ub=np.array(Aub), b_ub=np.array(bub))
    lo, hi = linprog(obj, **kw), linprog(-obj, **kw)
    return (lo.fun, -hi.fun) if lo.status == 0 and hi.status == 0 else None


inside_ok, nest_ok, n_cases = True, True, 0
for _ in range(60):
    C, J = 3, 15
    xs = np.sort(rng.uniform(0.05, 0.95, J))
    N = rng.integers(50, 2000, J).astype(float)
    # (a) monotone nonlinear m: row v col 0 increasing logistic, row n col 0 increasing step-ish
    m = np.zeros((J, 2, C))
    a1 = 0.15 + 0.6 / (1 + np.exp(-10 * (xs - rng.uniform(.3, .7))))
    a0 = 0.05 + 0.3 * xs ** 3
    for j in range(J):
        rest1 = rng.dirichlet([1, 1]); rest0 = rng.dirichlet([1, 1])
        m[j, 0] = [a1[j], *(1 - a1[j]) * rest1]
        m[j, 1] = [a0[j], *(1 - a0[j]) * rest0]
    tau = np.einsum("j,jc->jc", xs, m[:, 0]) + np.einsum("j,jc->jc", 1 - xs, m[:, 1])
    wts = N * xs / np.sum(N * xs)
    truth = np.sum(wts * m[:, 0, 0])
    b = n2_bounds(xs, tau, wts, (0, 0), {(0, 0): 1, (1, 0): 1}, C)
    inside_ok &= b is not None and b[0] - 1e-9 <= truth <= b[1] + 1e-9
    # (b) linear H with nonnegative slopes in column 0 for both rows; compare outer sets
    Hs = rng.dirichlet([1] * C, (2, 2)).transpose(2, 0, 1)        # H[c, r, s]; s: (v-share, n-share)
    # slope of row r in x (= share of v): H[c, r, 0] - H[c, r, 1]; force >= 0 in col 0 by sorting
    for r in range(2):
        hi_, lo_ = max(Hs[0, r, 0], Hs[0, r, 1]), min(Hs[0, r, 0], Hs[0, r, 1])
        diff_hi = hi_ - Hs[0, r, 0]
        Hs[0, r, 0], Hs[0, r, 1] = hi_, lo_
        # keep columns summing to one: move the change to column 1
        Hs[1, r, 0] -= diff_hi; Hs[1, r, 1] += diff_hi
    if (Hs < 0).any():
        continue
    Xm = np.column_stack([xs, 1 - xs])
    mlin = m_of(Hs, Xm)
    taul = np.einsum("jr,jrc->jc", Xm, mlin)
    Q = 0.5 * (Hs + Hs.transpose(0, 2, 1))
    # M0 outer set in kappa_c (c < C), vertices = observed range; plus slope signs in column 0
    lo_x, hi_x = xs.min(), xs.max()
    verts = [np.array([lo_x, 1 - lo_x]), np.array([hi_x, 1 - hi_x])]
    def Hk(kap):
        K = np.zeros((C, 2, 2))
        for c in range(C - 1):
            K[c, 0, 1], K[c, 1, 0] = kap[c], -kap[c]
        K[C - 1] = -K[:C - 1].sum(0)
        return Q + K
    Aub, bub = [], []
    E = np.eye(C - 1)
    def Hlin(c):       # H_c as affine function of kappa: returns (const, coef[(C-1)])
        const = Q[c]
        coef = np.zeros((C - 1, 2, 2))
        for cc in range(C - 1):
            K = np.zeros((2, 2))
            if c < C - 1 and cc == c:
                K[0, 1], K[1, 0] = 1, -1
            if c == C - 1:
                K[0, 1], K[1, 0] = -1, 1
            coef[cc] = K
        return const, coef
    for c in range(C):
        const, coef = Hlin(c)
        for v in verts:
            for r in range(2):
                # (const + sum kap coef) v >= 0
                Aub.append(-np.array([coef[cc][r] @ v for cc in range(C - 1)])); bub.append(const[r] @ v)
    for r in range(2):     # slope of row r in column 0: H[0,r,0]-H[0,r,1] >= 0
        const, coef = Hlin(0)
        Aub.append(-np.array([coef[cc][r, 0] - coef[cc][r, 1] for cc in range(C - 1)]))
        bub.append(const[r, 0] - const[r, 1])
    # target: B~_{v,0} = sum_j wts_j m_{j,v,0}(kappa) affine
    const, coef = Hlin(0)
    g0 = np.sum(wts * (Xm @ const[0]))
    gk = np.array([np.sum(wts * (Xm @ coef[cc][0])) for cc in range(C - 1)])
    kw = dict(A_ub=np.array(Aub), b_ub=np.array(bub), bounds=[(None, None)] * (C - 1), method="highs")
    l1, h1 = linprog(gk, **kw), linprog(-gk, **kw)
    if l1.status != 0 or h1.status != 0:
        continue
    m0b = (g0 + l1.fun, g0 - h1.fun)
    n2b = n2_bounds(xs, taul, wts, (0, 0), {(0, 0): 1, (1, 0): 1}, C)
    n_cases += 1
    nest_ok &= n2b is not None and n2b[0] <= m0b[0] + 1e-9 and m0b[1] <= n2b[1] + 1e-9
check("C5(a) N2 outer bound contains the truth under a monotone nonlinear DGP", inside_ok)
check("C5(b) M0 + signed slopes (outer) is inside N2 (outer) when tau is quadratic", nest_ok and n_cases > 20,
      f"{n_cases} cases")
X3 = rng.dirichlet([2, 2, 2], 400)
w = np.ones(400) / 400
m3 = X3[:, 2]                                   # constant (weakly nondecreasing) in x_0
check("C5(c) R=3: m nondecreasing in x_0 but Cov(m, X_0) < 0 (no Chebyshev nesting)",
      wcov(m3, X3[:, 0], w) < 0, f"Cov = {wcov(m3, X3[:, 0], w):.4f}")
# C5(c') (added after review, Batch B): the W1 association uses the tilted weights w^r ~ X_0 and
# w^{not r} ~ 1 - X_0. E[X_3 | X_0] = (1 - X_0)/2 is preserved by both tilts, so both covariances are negative.
cw1, cw0 = wcov(m3, X3[:, 0], X3[:, 0]), wcov(m3, X3[:, 0], 1 - X3[:, 0])
check("C5(c') R=3: the same m has negative covariance with X_0 under BOTH W1 weights", cw1 < 0 and cw0 < 0,
      f"Cov_w1 = {cw1:.4f}, Cov_w(not1) = {cw0:.4f}")


# ---------------------------------------------------------------------------
# C6. Local version (Prop N2.3): tau'(x) = D(x) + dW(x); local CR <=> D(x)
#     between 0 and tau'(x).
# ---------------------------------------------------------------------------
ok = True
for _ in range(500):
    x = rng.uniform(0.01, 0.99)
    m1, m0, d1, d0 = rng.uniform(0, 1, 2).tolist() + rng.normal(0, 1, 2).tolist()
    h = 1e-6
    tau = lambda z: z * (m1 + d1 * (z - x)) + (1 - z) * (m0 + d0 * (z - x))
    tp = (tau(x + h) - tau(x - h)) / (2 * h)
    D, dW = m1 - m0, x * d1 + (1 - x) * d0
    ok &= abs(tp - (D + dW)) < 1e-6
    cr = D * dW >= 0
    between = (0 <= D <= tp) or (tp <= D <= 0)
    ok &= cr == between
check("C6 Prop N2.3 (structural tau_m, differentiable m): local identity and local CR <=> D(x) between 0 and tau_m'(x)", ok)


# ===========================================================================
# Added after the Batch B review (Nov 2026). Each check fails on the old statement and passes on the new one
# (or, where marked, strengthens a check without anything failing).
# ===========================================================================
from fractions import Fraction as Fr
from numpy.polynomial import Polynomial as Poly


def int_gaps(sizes, ones, yes):
    """All realised gaps D = B_1 - B_2 over integer tables: unit i has sizes[i] people, ones[i] in group 1 and
    yes[i] positive outcomes. Returns the set of exact Fractions."""
    N1, N0 = sum(ones), sum(n - o for n, o in zip(sizes, ones))
    legal = [range(max(0, y - (n - o)), min(o, y) + 1) for n, o, y in zip(sizes, ones, yes)]
    return {Fr(sum(z), N1) - Fr(sum(yes) - sum(z), N0) for z in itertools.product(*legal)}


def two_by_two(sizes, ones, yes):
    N = np.array(sizes, float); X = np.array(ones) / N; T = np.array(yes) / N
    beta = np.zeros((len(N), 2, 2))
    beta[:, 0, 0] = beta[:, 1, 0] = T; beta[:, :, 1] = 1 - beta[:, :, 0]   # any valid table: margins only matter
    return Pop(N, np.column_stack([X, 1 - X]), beta)


def cr_branches(P):
    tgt = P.B(0, 0) - P.B(1, 0)
    dW, dB = P.dW_rest(0, 0), P.dB_rest(0, 0)
    return tgt, P.lp_range(tgt, rows_ge=[(dW, 0), (dB, 0)]), P.lp_range(tgt, rows_le=[(dW, 0), (dB, 0)])


# C7. M2 is sharp for CONTINUOUS tables only (P-B1, P-E11). Old claim: every value in phi(F) is attained by a
# unit table. Both referees' examples: the continuous CR set is {0}, no integer table satisfies CR; and under
# a positive between sign the continuous lower end is 0 but the smallest integer gap is 2/15.
ok7 = True; det = []
for name, sizes, ones, yes in (("PI: 3 units of 4", [4, 4, 4], [1, 2, 3], [0, 1, 0]),
                               ("EI: 2 precincts of 4", [4, 4], [1, 2], [1, 1])):
    P = two_by_two(sizes, ones, yes)
    tgt, plus, minus = cr_branches(P)
    cont = [b for b in (plus, minus) if b is not None]
    cont_has_zero = any(lo - 1e-9 <= 0 <= hi + 1e-9 for lo, hi in cont)
    ig = int_gaps(sizes, ones, yes)
    t = P.T[:, 0]; x = P.X[:, 0]
    ER = wcov(t, x, P.pw) / wcov(x, x, P.pw)
    gam = wcov(x, x, P.pw) / (P.mu[0] * (1 - P.mu[0]))
    lo_cr, hi_cr = sorted([gam * ER, ER])
    int_cr = [g for g in ig if lo_cr - 1e-12 <= float(g) <= hi_cr + 1e-12]
    ok7 &= cont_has_zero and not int_cr
    det.append(f"{name}: integer gaps {sorted(str(g) for g in ig)}, continuous CR contains 0: {cont_has_zero}, "
               f"integer CR {'empty' if not int_cr else int_cr}")
    if name.startswith("EI"):
        dB = P.dB_rest(0, 0)
        b_ge = P.lp_range(tgt, rows_ge=[(dB, 0)])
        NMd = np.sum(P.pw * x * t) / P.mu[0] - np.sum(P.pw * (1 - x) * t) / (1 - P.mu[0])
        int_lo = min(g for g in ig if float(g) >= NMd - 1e-12)
        ok7 &= abs(b_ge[0]) < 1e-9 and int_lo == Fr(2, 15)
        det.append(f"positive between sign: continuous lower end {b_ge[0]:.3f}, integer lower end {int_lo}")
check("C7 Remark M2.1: continuous N1 accepts CR (and lower end 0) while no integer table does "
      "(old 'attained by unit tables' fails for integer counts)", ok7, "; ".join(det))

# C7a. Remark M2.1, lattice statement (2x2, margins only, sign rows on D): the integer gap set is the
# continuous DD interval intersected with the lattice -Y/N0 + j (1/N1 + 1/N0); endpoints differ by less than
# one step 1/(N mu (1-mu)). Brute force on tiny populations.
ok7a, n7a = True, 0
for _ in range(300):
    p_ = int(rng.integers(2, 5)); sizes = rng.integers(1, 6, p_).tolist()
    ones = [int(rng.integers(0, n + 1)) for n in sizes]; yes = [int(rng.integers(0, n + 1)) for n in sizes]
    N1, N0 = sum(ones), sum(sizes) - sum(ones)
    if N1 == 0 or N0 == 0:
        continue
    ig = sorted(int_gaps(sizes, ones, yes))
    step = Fr(1, N1) + Fr(1, N0)
    Ntot = sum(sizes); mu_ = Fr(N1, Ntot)
    P = two_by_two(sizes, ones, yes)
    lo, hi = P.lp_range(P.B(0, 0) - P.B(1, 0))
    full = all(b - a == step for a, b in zip(ig[:-1], ig[1:]))
    ok7a &= full and abs(float(ig[0]) - lo) < 1e-9 and abs(float(ig[-1]) - hi) < 1e-9 \
        and step == 1 / (Ntot * mu_ * (1 - mu_))
    n7a += 1
check("C7a Remark M2.1 (2x2, margins only): integer gap set = DD interval on the lattice of step "
      "1/(N mu(1-mu)), endpoints attained", ok7a, f"{n7a} tiny populations")

# C7b. M2(d): with two slack dials (within signs on two columns, 2x3) the set of (eps1, eps2) that preserves a
# conclusion (here: the width of the B_{1,1} interval stays <= c) is a frontier, not a
# product of intervals, so the old 'the breakdown value is well defined' needs a scalar path.
found = False
for _ in range(60):
    P = random_pop(p=10, R=2, C=3)
    tgt = P.B(0, 0)
    W0, W1 = P.dW_rest(0, 0), P.dW_rest(0, 1)
    if P.lp_range(tgt, rows_ge=[(W0, 0), (W1, 0)]) is None:
        continue
    grid = np.linspace(0, 0.05, 7)
    U = np.array([[(lambda r: np.nan if r is None else r[1] - r[0])(P.lp_range(tgt, rows_ge=[(W0, -e1), (W1, -e2)]))
                   for e2 in grid] for e1 in grid])
    if np.isnan(U).any():
        continue
    if not (U[6, 0] > U[0, 0] + 1e-6 and U[0, 6] > U[0, 0] + 1e-6):   # both dials must matter
        continue
    c = U[3, 3]
    S = U <= c + 1e-12
    prod = np.outer(S.any(1), S.any(0))
    if (prod & ~S).any():
        found = True
        break
check("C7b M2(d): the region of (eps1, eps2) keeping 'width <= c' is not a product of intervals "
      "(breakdown frontier)", found)

# C8. M3(d) at gamma = 1 (P-B3, P-E13): two pure units (0,0), (1,1): both CR branches are non-empty,
# contradicting 'exactly one'. Under the corrected statement (gamma < 1 required) both give {D_ER}.
P = Pop(np.array([1.0, 1.0]), np.array([[0.0, 1.0], [1.0, 0.0]]),
        np.array([[[0.0, 1.0], [0.0, 1.0]], [[1.0, 0.0], [1.0, 0.0]]]))
tgt, plus, minus = cr_branches(P)
x = P.X[:, 0]; t = P.T[:, 0]
ER = wcov(t, x, P.pw) / wcov(x, x, P.pw); gam = wcov(x, x, P.pw) / (P.mu[0] * (1 - P.mu[0]))
check("C8 M3(d): gamma = 1 (pure units): both CR branches non-empty and equal to {D_ER} (old 'exactly one' fails)",
      plus is not None and minus is not None and abs(gam - 1) < 1e-12 and
      max(abs(plus[0] - ER), abs(plus[1] - ER), abs(minus[0] - ER), abs(minus[1] - ER)) < 1e-9,
      f"gamma = {gam:.3f}, ER = {ER:.3f}, branches {plus}, {minus}")

# C8b. M3(c) needs 'no unit-level restrictions' (P-B4, P-E12): a unit ordering beta_i1 >= beta_i2 + 0.1 removes
# the equal-rates witness and makes {delta_B <= 0} empty.
P = random_pop(p=15)
P.T[:, 0] = np.clip(P.T[:, 0], 0.2, 0.8); P.T[:, 1] = 1 - P.T[:, 0]
tgt = P.B(0, 0) - P.B(1, 0)
dB = P.dB_rest(0, 0)
Aeq, beq = P.equalities()
Aub = [P.dB_rest(0, 0)]; bub = [0.0]                       # delta_B <= 0
for i in range(P.p):
    a = np.zeros(P.nvar()); a[P.idx(i, 1, 0)] = 1; a[P.idx(i, 0, 0)] = -1
    Aub.append(a); bub.append(-0.1)                        # beta_i2 - beta_i1 <= -0.1
r_ord = linprog(tgt, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=Aeq, b_eq=beq, bounds=[(0, 1)] * P.nvar(), method="highs")
r_noord = P.lp_range(tgt, rows_le=[(dB, 0)])
Aub2 = np.array(Aub[1:]); r_ordonly = linprog(tgt, A_ub=Aub2, b_ub=np.array(bub[1:]), A_eq=Aeq, b_eq=beq,
                                             bounds=[(0, 1)] * P.nvar(), method="highs")
check("C8b M3(c): with a unit ordering the delta_B <= 0 set is empty although each part alone is feasible "
      "(old 'never refuted' fails once U_i are present)", r_ord.status == 2 and r_noord is not None and r_ordonly.status == 0)

# C9. M5(c) under M1 (P-B7): six equally weighted units, districts {.1,.2,.3} and {.7,.8,.9}, both rows' mean
# .8 + .1x in district 1 and .2 + .1x in district 2 (common positive slopes, valid rates). Pooled delta*_W < 0,
# contradicting the old 'both contextual slopes >= 0 => delta*_W >= 0' for M1. The corrected M5(c) adds the
# between-district term: delta*_W = sum_s mu_s [ E_{w^s} Cov(X, m_s | D) + Cov_{w^s}(E[X|D], E[m_s|D]) ].
x6 = np.array([.1, .2, .3, .7, .8, .9]); d6 = np.array([0, 0, 0, 1, 1, 1]); p6 = np.ones(6) / 6
m6 = np.where(d6 == 0, .8, .2) + .1 * x6
mu6 = p6 @ x6
def dW2(x, p, m1, m0):
    mu = p @ x
    return mu * wcov(x, m1, p * x) + (1 - mu) * wcov(x, m0, p * (1 - x))
pooled = dW2(x6, p6, m6, m6)
def decomp(x, p, d, m1, m0):
    tot = 0.0
    for wts, m in ((p * x, m1), (p * (1 - x), m0)):
        share = wts.sum()
        wq = wts / share
        within = sum(wq[d == g].sum() * wcov(x[d == g], m[d == g], wq[d == g]) for g in (0, 1))
        ex = np.array([np.sum(wq[d == g] * x[d == g]) / wq[d == g].sum() for g in d])
        em = np.array([np.sum(wq[d == g] * m[d == g]) / wq[d == g].sum() for g in d])
        tot += share * (within + wcov(ex, em, wq))
    return tot
dec = decomp(x6, p6, d6, m6, m6)
per_district = [dW2(x6[d6 == g], p6[d6 == g] / p6[d6 == g].sum(), m6[d6 == g], m6[d6 == g]) for g in (0, 1)]
check("C9 M5(c): common positive slopes under M1, pooled delta*_W < 0 (old M5(c) fails); the corrected "
      "decomposition with the between-district term is exact; within each district delta*_W > 0",
      pooled < 0 and abs(pooled + 0.0492711) < 1e-6 and abs(dec - pooled) < 1e-12 and min(per_district) > 0,
      f"pooled {pooled:.7f}; decomposition {dec:.7f}; per district {per_district[0]:.5f}, {per_district[1]:.5f}")

# C10. N2.2(b) needs the complete finite support (P-B9). True support {1/4, 1/2, 3/4}, row means (.4,.4,.4)
# and (8/15,.7,.8): tau = (.5,.55,.5). On a grid keeping only the extremes, (.5,.5) at both extremes passes the
# LP, but no monotone extension reaches tau(1/2) = .55: the old 'every point extends; nothing lost' fails.
xs3 = np.array([.25, .5, .75]); m1_3 = np.array([.4, .4, .4]); m0_3 = np.array([8 / 15, .7, .8])
tau3 = xs3 * m1_3 + (1 - xs3) * m0_3
grid_point_ok = all(abs(x_ * .5 + (1 - x_) * .5 - tau3[j]) < 1e-12 for j, x_ in ((0, .25), (2, .75)))
# extension: m1(1/2), m0(1/2) with .5 <= m(1/2) <= .5 (both rows nondecreasing) and the accounting at 1/2
ext = linprog(np.zeros(2), A_eq=[[.5, .5]], b_eq=[tau3[1]], bounds=[(.5, .5), (.5, .5)], method="highs")
check("C10 N2.2(b): a grid point on a strict subset of the support has no monotone extension "
      "(old 'restricting to support points loses nothing' fails on a grid)",
      np.allclose(tau3, [.5, .55, .5]) and grid_point_ok and ext.status == 2)

# C10b. N2.2(a): conditioning on X alone is not the sharp version when size is informative (P-B10). X = 1/2,
# two equiprobable sizes; size 1: T in {0,1} equiprobably; size 2: T = 1/2. The model E[beta|X,N] = m(X) forces
# m(1/2) = (.5,.5) (both rates equal T at size 1). E[S(X,T)|X] and S(x, tau) both admit (.75, .25).
def in_S(b, x_, t_):
    return 0 <= b[0] <= 1 and 0 <= b[1] <= 1 and abs(x_ * b[0] + (1 - x_) * b[1] - t_) < 1e-12
cand = (.75, .25)
# E[S|X=1/2] = 1/2 E[S|size 1] + 1/2 S(1/2, 1/2), with E[S|size 1] = {(.5,.5)}
x_cond = in_S((2 * cand[0] - .5, 2 * cand[1] - .5), .5, .5)
s_tau = in_S(cand, .5, .5)
v_cond = cand == (.5, .5)                                    # size 1 alone forces (.5, .5)
check("C10b N2.2(a): the X-conditional set and S(x, tau) admit (.75,.25), which the V=(X,N) model excludes",
      x_cond and s_tau and not v_cond)

# C11. N2.3 as an observable restriction (P-B11): on support {1/4,1/2,3/4} add h = (x-1/4)(x-1/2)(x-3/4) to both
# rows of m1 = .7, m0 = .3. The law on the support is unchanged, but tau'(1/2) moves by -1/16. Model 1 has local
# within effect 0 and gap .4, so local CR holds with D = tau'(1/2); the other compatible derivative gives the
# threshold .3375 < .4, so an LP using it excludes the true gap.
z = Poly([0, 1]); h = (z - .25) * (z - .5) * (z - .75)
grid01 = np.linspace(0, 1, 1001)
valid2 = np.all((.7 + h(grid01) >= 0) & (.7 + h(grid01) <= 1) & (.3 + h(grid01) >= 0) & (.3 + h(grid01) <= 1))
tau1 = .3 + .4 * z; tau2 = z * (.7 + h) + (1 - z) * (.3 + h)
same_law = np.allclose(tau1(xs3), tau2(xs3)) and np.allclose(h(xs3), 0)
d1, d2 = tau1.deriv()(.5), tau2.deriv()(.5)
check("C11 N2.3: two compatible models, same law on the support, local thresholds .4 and .3375; the second "
      "excludes the first model's gap .4 (old N2.3 as an observable restriction fails)",
      valid2 and same_law and abs(d2 - d1 + 1 / 16) < 1e-12 and d2 < .4 <= d1 + 1e-12,
      f"tau'(1/2) = {d1:.4f} vs {d2:.4f}")

# C12. Section 6 binning (P-B13): equal bin weights at x = .25, .75, m1 = x, m0 = 0. E[T|bin] = E[X^2] = .3125,
# but the ordinary bin-average rates at E[X] = .5 give .25: averaging Lemma N2.1 with [CE-cvx] fails. The
# person-weighted bin rate E_b[X m1]/E_b[X] = .625 reproduces .3125 (Lemma N2.4).
xb = np.array([.25, .75]); m1b = xb; m0b = np.zeros(2)
ET = np.mean(xb * m1b + (1 - xb) * m0b); xbar = xb.mean()
old_bin = xbar * m1b.mean() + (1 - xbar) * m0b.mean()
m1bar = np.mean(xb * m1b) / xbar; m0bar = np.mean((1 - xb) * m0b) / (1 - xbar)
new_bin = xbar * m1bar + (1 - xbar) * m0bar
# Lemma N2.4 on random monotone designs: bin rates reproduce E_b T, the row targets, and the ordering
ok12 = True
for _ in range(200):
    J = 30; xs_ = np.sort(rng.uniform(0, 1, J)); pw_ = rng.uniform(1, 5, J); pw_ /= pw_.sum()
    a1 = np.sort(rng.uniform(0, 1, J)); a0 = np.sort(rng.uniform(0, 1, J))      # nondecreasing rows
    cuts = np.sort(rng.choice(np.arange(1, J), 4, replace=False)); bins = np.split(np.arange(J), cuts)
    r1, r0, xbars = [], [], []
    for b_ in bins:
        w_ = pw_[b_] / pw_[b_].sum(); xb_ = w_ @ xs_[b_]
        r1.append(w_ @ (xs_[b_] * a1[b_]) / xb_); r0.append(w_ @ ((1 - xs_[b_]) * a0[b_]) / (1 - xb_)); xbars.append(xb_)
        ok12 &= abs(xb_ * r1[-1] + (1 - xb_) * r0[-1] - w_ @ (xs_[b_] * a1[b_] + (1 - xs_[b_]) * a0[b_])) < 1e-12
    Pb = np.array([pw_[b_].sum() for b_ in bins])
    ok12 &= abs(np.sum(Pb * np.array(xbars) * np.array(r1)) - pw_ @ (xs_ * a1)) < 1e-12
    ok12 &= np.all(np.diff(r1) >= -1e-12) and np.all(np.diff(r0) >= -1e-12)
check("C12 Lemma N2.4: ordinary bin averages break accounting (.25 vs .3125; old binning step fails); "
      "person-weighted bin rates reproduce it (.625) and keep row targets and monotone order",
      abs(ET - .3125) < 1e-12 and abs(old_bin - .25) < 1e-12 and abs(m1bar - .625) < 1e-12
      and abs(new_bin - ET) < 1e-12 and ok12)

# C12b. Section 6, degree selection (P-B14): LOO choice between intercept-only and intercept+slope in a balanced
# two-level design (true slope 0). Given that the slope model is chosen, the nominal 95% slope interval covers
# about 1 - .05/P(|Z| > sqrt 2) = 0.682.
from scipy.stats import norm as _norm
nb, reps_b = 200, 40000
xx = np.repeat([-1.0, 1.0], nb // 2)
E_ = rng.standard_normal((reps_b, nb))
b1 = E_ @ xx / nb                                          # slope estimate, var 1/n (sigma known = 1)
rss1 = np.sum((E_ - E_.mean(1, keepdims=True)) ** 2, 1)
rss2 = np.sum((E_ - E_.mean(1, keepdims=True) - b1[:, None] * xx) ** 2, 1)
sel = rss2 / (1 - 2 / nb) ** 2 < rss1 / (1 - 1 / nb) ** 2
cover = np.abs(b1) * np.sqrt(nb) <= 1.959964
cond = cover[sel].mean()
theory = 1 - .05 / (2 * _norm.sf(np.sqrt(2)))
check("C12b Section 6: fixed-degree inference after LOO degree selection under-covers (conditional coverage "
      "~0.68, not 0.95)", abs(cond - theory) < 0.02 and cond < 0.8,
      f"conditional coverage {cond:.3f} (theory {theory:.3f}); selected in {sel.mean():.3f}")

# C13. N2.3(d), smooth witness (P-B12; strengthens 1b V6, nothing fails): X ~ U[1/5, 2/5], m1 = x - x^2/2,
# m0 = 1 - x^2/2. Local within effect x m1' + (1-x) m0' = 0 everywhere; aggregate delta*_W = -1/47250
# (finite design with equal mass at the endpoints: -1/5250).
m1s = z - z * z / 2; m0s = 1 - z * z / 2
loc0 = np.allclose((z * m1s.deriv() + (1 - z) * m0s.deriv()).coef, 0)
def Eu(f):
    F_ = f.integ(); return 5 * (F_(.4) - F_(.2))
def wterm(wf, m):
    return Eu(wf * z * m) - Eu(wf * z) * Eu(wf * m) / Eu(wf)
dW_cont = wterm(z, m1s) + wterm(1 - z, m0s)
xe = np.array([.2, .4]); pe = np.array([.5, .5])
dW_fin = dW2(xe, pe, m1s(xe), m0s(xe))
check("C13 N2.3(d): smooth rows with zero local within effect everywhere, aggregate delta*_W < 0",
      loc0 and abs(dW_cont + 1 / 47250) < 1e-14 and abs(dW_fin + 1 / 5250) < 1e-14,
      f"continuous {dW_cont:.4e} (-1/47250), finite {dW_fin:.4e} (-1/5250)")

# C14. Definition 1, W2 is ONE pooled sum (P-E1): two units of equal size; every original group's rate rises by
# 2 points from unit 1 to unit 2, yet W1 = -0.02633 and W2 = +0.00110. The old label 'each group separately'
# would make W2 a per-group direction; it is not, and W1 and W2 can disagree in sign.
X14 = np.array([[.2, .7, .1], [.6, .05, .35]]); m14 = np.array([[.4, .8, .1], [.42, .82, .12]]); p14 = np.array([.5, .5])
beta14 = np.stack([m14, 1 - m14], axis=2)
P14 = Pop(np.array([1.0, 1.0]), X14, beta14)
W1_14 = P14.value(P14.dW_rest(0, 0))
W2_14 = sum(P14.mu[s] * wcov(X14[:, 0], m14[:, s], p14 * X14[:, s]) for s in range(3))
check("C14 Def. 1: every group's rate rises, W1 < 0 < W2 (W2 is one pooled sum, not a per-group direction)",
      abs(W1_14 + 0.02633) < 5e-6 and abs(W2_14 - 0.00110) < 5e-6, f"W1 {W1_14:.5f}, W2 {W2_14:.5f}")

print(f"\n{sum(results)}/{len(results)} PASS")
