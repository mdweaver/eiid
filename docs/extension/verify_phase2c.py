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
check("C2 Prop M3(d): CR never empty, exactly one branch (ER != 0)", cr_ok)
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
    # ordering: build H with row r dominating row s elementwise in column c=0
    Hord = H0.copy()
    Hord[0, 0, :] = np.maximum(H0[0, 0, :], H0[0, 1, :])
    Hord[0, 1, :] = np.minimum(H0[0, 0, :], H0[0, 1, :])
    m = m_of(Hord, X)
    order_ok &= np.sum(pw * X[:, 0] * X[:, 1] * (m[:, 0, 0] - m[:, 1, 0])) >= -1e-14
check("C4 Prop M5: superpopulation associations affine in k (R=C=3)", worst < 1e-12, f"max err {worst:.1e}")
check("C4 Prop M5(a): mean-level ordering => dB*(r,s) >= 0", order_ok)


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
check("C6 Prop N2.3: local identity and local CR <=> D(x) between 0 and tau'(x)", ok)

print(f"\n{sum(results)}/{len(results)} PASS")
