"""Numerical checks of the facts used in PHASE1_monotone_did.md (Phase-1 proposal).

Run:  python3 -I verify_phase1_mei_did.py  > verify_phase1_mei_did_output.txt

Each check prints PASS/FAIL. Checks are on finite populations of units
(neighbourhoods) with sizes N_i, group-1 share X_i and unit rates beta_i, i.e.
Elzayn, Goldin, Guage, Ho & Morton's setting with p_n = N_i / sum N.
"""
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(20261029)
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))


def wmean(a, w):
    return np.sum(w * a) / np.sum(w)


def wcov(a, b, w):
    return wmean(a * b, w) - wmean(a, w) * wmean(b, w)


def mei_objects(N, X, b1, b0):
    """Population objects of the 2-group problem (group 1 share X)."""
    p = N / N.sum()
    T = X * b1 + (1 - X) * b0
    EX = np.sum(p * X)
    VarX = EX * (1 - EX)                      # individual-level variance of X
    VarXN = wcov(X, X, p)                     # between-unit variance
    Y1 = np.sum(p * X * b1) / EX
    Y0 = np.sum(p * (1 - X) * b0) / (1 - EX)
    D = Y1 - Y0
    dB = np.sum(p * X * (1 - X) * (b1 - b0))
    dW = EX * wcov(b1, X, p * X) + (1 - EX) * wcov(b0, X, p * (1 - X))
    ER = wcov(T, X, p) / VarXN
    ER1 = wmean(T, p) + ER * (1 - EX)
    ER0 = wmean(T, p) - ER * EX
    NM1 = np.sum(p * X * T) / EX
    NM0 = np.sum(p * (1 - X) * T) / (1 - EX)
    NM = NM1 - NM0
    gamma = VarXN / VarX
    return dict(p=p, T=T, EX=EX, VarX=VarX, VarXN=VarXN, Y1=Y1, Y0=Y0, D=D, dB=dB,
                dW=dW, ER=ER, ER1=ER1, ER0=ER0, NM=NM, NM1=NM1, NM0=NM0, gamma=gamma)


def random_pop(p=40, lo=0.0, hi=1.0):
    N = rng.integers(50, 2000, p).astype(float)
    X = rng.beta(2, 3, p) * 0.9 + 0.05
    b1 = rng.uniform(lo, hi, p)
    b0 = rng.uniform(lo, hi, p)
    return N, X, b1, b0


# ---------------------------------------------------------------------------
# V1. Elzayn et al. identities (any bounded outcome)
# ---------------------------------------------------------------------------
worst = 0.0
for _ in range(300):
    N, X, b1, b0 = random_pop(lo=-1.0, hi=1.0)       # outcome in [-1, 1] (e.g. a change)
    o = mei_objects(N, X, b1, b0)
    errs = [o["ER"] - o["D"] - o["dW"] / o["VarXN"],
            o["NM"] - o["D"] + o["dB"] / o["VarX"],
            o["D"] * (1 - o["gamma"]) * o["VarX"] - (o["dB"] + o["dW"]),
            o["NM"] - o["gamma"] * o["ER"],
            o["ER1"] - o["Y1"] - o["dW"] * (1 - o["EX"]) / o["VarXN"],
            o["ER0"] - o["Y0"] + o["dW"] * o["EX"] / o["VarXN"],
            o["NM1"] - o["Y1"] + o["dB"] / o["EX"],
            o["NM0"] - o["Y0"] - o["dB"] / (1 - o["EX"])]
    worst = max(worst, max(abs(e) for e in errs))
check("V1 ER/NM bias identities, decomposition, NM = gamma*ER (300 pops, Y in [-1,1])",
      worst < 1e-10, f"max abs error {worst:.1e}")


# ---------------------------------------------------------------------------
# V2. Closed-form monotone bounds for D = LP over unit-level rates (sharpness
#     in the finite-population sense), for dW>=0, dW<=0, dB>=0, dB<=0, CR.
# ---------------------------------------------------------------------------
def lp_bounds_D(N, X, T, extra):
    """min/max D over beta1_i in its DD interval, beta0 implied, plus linear rows.
    extra: list of (coef_vector_on_beta1, const, sense) meaning coef@b1 + const sense 0."""
    p = N / N.sum()
    L = np.maximum(0, (T - (1 - X)) / X)
    U = np.minimum(1, T / X)
    EX = np.sum(p * X)
    # D = sum p X b1 /EX - sum p (T - X b1) /(1-EX)
    cD = p * X * (1 / EX + 1 / (1 - EX))
    kD = -np.sum(p * T) / (1 - EX)
    A, b = [], []
    for coef, const, sense in extra:
        if sense == ">=":
            A.append(-coef); b.append(const)
        else:
            A.append(coef); b.append(-const)
    A = np.array(A) if A else None
    b = np.array(b) if b else None
    out = []
    for s in (1, -1):
        r = linprog(s * cD, A_ub=A, b_ub=b, bounds=list(zip(L, U)), method="highs")
        if r.status != 0:
            return None
        out.append(s * r.fun + kD)
    return out


def linear_forms(N, X, T):
    """dW and dB as affine functions of the unit vector beta1 (beta0 implied)."""
    p = N / N.sum()
    EX = np.sum(p * X)
    VarX = EX * (1 - EX)
    VarXN = wcov(X, X, p)
    cD = p * X * (1 / EX + 1 / (1 - EX))
    kD = -np.sum(p * T) / (1 - EX)
    ER = wcov(T, X, p) / VarXN
    NM = np.sum(p * X * T) / EX - np.sum(p * (1 - X) * T) / (1 - EX)
    # dW = (ER - D) VarXN ; dB = (D - NM) VarX
    dW = (-cD * VarXN, (ER - kD) * VarXN)
    dB = (cD * VarX, (kD - NM) * VarX)
    return dW, dB, ER, NM


worst, n_cases, n_empty = 0.0, 0, 0
for _ in range(200):
    N, X, b1, b0 = random_pop()
    o = mei_objects(N, X, b1, b0)
    T = o["T"]
    p, EX, VarX = o["p"], o["EX"], o["VarX"]
    L = np.maximum(0, (T - (1 - X)) / X); U = np.minimum(1, T / X)
    Lo = np.maximum(0, (T - X) / (1 - X)); Uo = np.minimum(1, T / (1 - X))
    mobm = np.sum(p * X * L) / EX - np.sum(p * (1 - X) * Uo) / (1 - EX)
    mobp = np.sum(p * X * U) / EX - np.sum(p * (1 - X) * Lo) / (1 - EX)
    dW, dB, ER, NM = linear_forms(N, X, T)
    cases = {
        "dW>=0": ([(dW[0], dW[1], ">=")], (mobm, min(mobp, ER))),
        "dW<=0": ([(dW[0], dW[1], "<=")], (max(mobm, ER), mobp)),
        "dB>=0": ([(dB[0], dB[1], ">=")], (NM, mobp)),
        "dB<=0": ([(dB[0], dB[1], "<=")], (mobm, NM)),
        "both>=0": ([(dW[0], dW[1], ">="), (dB[0], dB[1], ">=")], (NM, min(mobp, ER))),
    }
    for nm, (rows, cf) in cases.items():
        lp = lp_bounds_D(N, X, T, rows)
        if cf[0] > cf[1] + 1e-12:
            n_empty += 1
            ok = lp is None
            worst = max(worst, 0 if ok else 1)
        else:
            n_cases += 1
            worst = max(worst, abs(lp[0] - cf[0]), abs(lp[1] - cf[1]))
    # contextual reinforcement: union of (both >=0) and (both <=0)
    a = lp_bounds_D(N, X, T, [(dW[0], dW[1], ">="), (dB[0], dB[1], ">=")])
    c = lp_bounds_D(N, X, T, [(dW[0], dW[1], "<="), (dB[0], dB[1], "<=")])
    cf = (NM, min(mobp, ER)) if ER >= NM else (max(mobm, ER), NM)
    lpu = a if a is not None else c
    if a is not None and c is not None:
        lpu = (min(a[0], c[0]), max(a[1], c[1]))
    n_cases += 1
    worst = max(worst, abs(lpu[0] - cf[0]), abs(lpu[1] - cf[1]))
check("V2 Elzayn closed forms (Prop.5, 7, 8 of arXiv v1) = unit-level LP",
      worst < 1e-8, f"{n_cases} non-empty cases, {n_empty} empty cases agree; max err {worst:.1e}")


# ---------------------------------------------------------------------------
# V3. Under 2x2 linear contextual effects (beta = m(X) exactly, so finite and
#     super-population objects coincide), dW and dB are affine in w1:
#     dW = a*b1 + b*w1 with b1 = w1 + d1,
#     a = E[X] Var_{w1}(X),  b = E[1-X] Var_{w0}(X);
#     dB = E_N[X(1-X)(c1 - w1 + d1 X)].
# ---------------------------------------------------------------------------
worst = 0.0
for _ in range(200):
    N = rng.integers(50, 2000, 60).astype(float)
    X = rng.uniform(0.05, 0.95, 60)
    w0, w1, b0, b1 = rng.uniform(-0.3, 0.3, 4) + np.array([0.5, 0, 0.5, 0])
    m1, m0 = b0 + b1 * X, w0 + w1 * X
    o = mei_objects(N, X, m1, m0)
    p, EX = o["p"], o["EX"]
    c1, d1 = w1 + b0 - w0, b1 - w1
    a = EX * wcov(X, X, p * X)
    bb = (1 - EX) * wcov(X, X, p * (1 - X))
    dW_cf = a * (w1 + d1) + bb * w1
    dB_cf = np.sum(p * X * (1 - X) * (c1 - w1 + d1 * X))
    worst = max(worst, abs(o["dW"] - dW_cf), abs(o["dB"] - dB_cf))
check("V3 under LCE, dW and dB are affine in the unidentified w1 (closed forms)",
      worst < 1e-12, f"max err {worst:.1e}")


# ---------------------------------------------------------------------------
# V4. R x C generalisation of the two bias formulas (one column c, R = 3):
#     Goodman (no intercept) g - B_c = M^{-1} sum_s mu_s Cov_{w_s}(X, beta_sc),
#     NM_rc - B_rc = sum_{s != r} E_p[X_r X_s (beta_sc - beta_rc)] / mu_r.
# ---------------------------------------------------------------------------
worst = 0.0
for _ in range(200):
    p_ = 50
    N = rng.integers(50, 2000, p_).astype(float)
    X = rng.dirichlet([2, 2, 1.5], p_)
    beta = rng.uniform(0, 1, (p_, 3))                       # beta[i, s] = beta_{i, s, c}
    w = N / N.sum()
    T = np.sum(X * beta, axis=1)
    mu = w @ X
    B = (w[:, None] * X * beta).sum(0) / mu
    M = (w[:, None, None] * X[:, :, None] * X[:, None, :]).sum(0)
    g = np.linalg.solve(M, (w[:, None] * X * T[:, None]).sum(0))
    Delta = sum(mu[s] * np.array([wcov(X[:, r], beta[:, s], w * X[:, s]) for r in range(3)])
                for s in range(3))
    worst = max(worst, np.max(np.abs(g - B - np.linalg.solve(M, Delta))), abs(Delta.sum()))
    for r in range(3):
        NMr = np.sum(w * X[:, r] * T) / mu[r]
        rhs = sum(np.sum(w * X[:, r] * X[:, s] * (beta[:, s] - beta[:, r]))
                  for s in range(3) if s != r) / mu[r]
        worst = max(worst, abs(NMr - B[r] - rhs))
check("V4 R x C bias formulas for Goodman regression and the neighbourhood model; 1'Delta = 0",
      worst < 1e-10, f"max err {worst:.1e}")


# ---------------------------------------------------------------------------
# V5. Monotone conditional means for both groups (Manski-Pepper-type) imply
#     dW >= 0 (Chebyshev's association inequality, any weights).
# ---------------------------------------------------------------------------
ok = True
for _ in range(500):
    N, X, _, _ = random_pop()
    order = np.argsort(X)
    m1 = np.empty_like(X); m0 = np.empty_like(X)
    m1[order] = np.sort(rng.uniform(0, 1, X.size))
    m0[order] = np.sort(rng.uniform(0, 1, X.size))
    ok &= mei_objects(N, X, m1, m0)["dW"] >= -1e-14
check("V5 nondecreasing m1 and m0 in X  =>  dW >= 0 (500 random cases)", bool(ok))


# ---------------------------------------------------------------------------
# V6. The local sign condition x m1'(x) + (1-x) m0'(x) >= 0 for all x does NOT
#     imply the aggregate dW >= 0: search for a counterexample by LP over
#     piecewise-linear m1, m0 on a grid (slopes are the variables).
# ---------------------------------------------------------------------------
K = 9
xs = np.linspace(0.05, 0.95, K)                    # unit compositions (one unit per point)
Nn = np.array([5, 1, 1, 1, 1, 1, 1, 1, 5], float)  # heavy mass at the extremes
p = Nn / Nn.sum()
EX = np.sum(p * xs)
# m_g(x_k) = m_g(x_0) + sum_{j<k} s_gj h ; variables: s1 (K-1), s0 (K-1), m1_0, m0_0
h = xs[1] - xs[0]
nv = 2 * (K - 1) + 2
def mvec(g, k):
    v = np.zeros(nv)
    v[(g * (K - 1)):(g * (K - 1) + k)] = h
    v[2 * (K - 1) + g] = 1.0
    return v
# dW is linear in m-values: EX*Cov_{w1}(m1,X) + (1-EX)*Cov_{w0}(m0,X)
w1 = p * xs / np.sum(p * xs); w0 = p * (1 - xs) / np.sum(p * (1 - xs))
obj = np.zeros(nv)
for k in range(K):
    obj += EX * w1[k] * (xs[k] - np.sum(w1 * xs)) * mvec(0, k)
    obj += (1 - EX) * w0[k] * (xs[k] - np.sum(w0 * xs)) * mvec(1, k)
A, b = [], []
for j in range(K - 1):               # local sign at both ends of each segment
    for x in (xs[j], xs[j + 1]):
        row = np.zeros(nv); row[j] = -x; row[(K - 1) + j] = -(1 - x)
        A.append(row); b.append(0.0)
for g in (0, 1):                     # 0 <= m <= 1 at grid points
    for k in range(K):
        A.append(mvec(g, k)); b.append(1.0)
        A.append(-mvec(g, k)); b.append(0.0)
r = linprog(obj, A_ub=np.array(A), b_ub=np.array(b), bounds=[(None, None)] * nv, method="highs")
check("V6 local within-sign everywhere does not imply aggregate dW >= 0 (counterexample exists)",
      r.status == 0 and r.fun < -1e-6, f"min aggregate dW = {r.fun:.4f}")


# ---------------------------------------------------------------------------
# V7. DiD: with the same units, N and X in both periods,
#  (a) the model-free DiD bound is the difference of the two periods' DD intervals
#      (product structure), checked against a joint LP;
#  (b) first-difference LCE with ARBITRARY unit effects (levels not linear)
#      plus Delta w1 = 0 point-identifies the DiD from the regression of
#      Delta T on X, X^2;
#  (c) monotone EI on changes: ER on Delta T minus DiD = dW(Delta)/Var(X_N).
# ---------------------------------------------------------------------------
worst_a = 0.0
for _ in range(100):
    N, X, b1a, b0a = random_pop(p=25)
    _, _, b1b, b0b = random_pop(p=25)
    Ta = X * b1a + (1 - X) * b0a
    Tb = X * b1b + (1 - X) * b0b
    pw = N / N.sum(); EX = np.sum(pw * X)
    def dd(T):
        L = np.maximum(0, (T - (1 - X)) / X); U = np.minimum(1, T / X)
        c = pw * X * (1 / EX + 1 / (1 - EX)); k = -np.sum(pw * T) / (1 - EX)
        return c, k, L, U
    ca, ka, La, Ua = dd(Ta); cb, kb, Lb, Ub = dd(Tb)
    obj = np.concatenate([cb, -ca])
    lo = linprog(obj, bounds=list(zip(Lb, Ub)) + list(zip(La, Ua)), method="highs").fun + kb - ka
    hi = -linprog(-obj, bounds=list(zip(Lb, Ub)) + list(zip(La, Ua)), method="highs").fun + kb - ka
    da = (ca @ La + ka, ca @ Ua + ka); db = (cb @ Lb + kb, cb @ Ub + kb)
    worst_a = max(worst_a, abs(lo - (db[0] - da[1])), abs(hi - (db[1] - da[0])))
check("V7a model-free DiD bound = difference of the periods' DD intervals", worst_a < 1e-9,
      f"max err {worst_a:.1e}")

worst_b = 0.0
for _ in range(200):
    pn = 80
    N = rng.integers(50, 2000, pn).astype(float)
    X = rng.uniform(0.05, 0.6, pn)
    # arbitrary unit effects zeta (not linear in X) shared by both periods
    zeta1 = 0.15 * np.sin(9 * X) + 0.05 * rng.standard_normal(pn)
    zeta0 = 0.15 * np.cos(7 * X) + 0.05 * rng.standard_normal(pn)
    w0, w1, b0, b1 = 0.4, 0.1, 0.3, 0.2
    dw0, db0, db1 = rng.uniform(-0.1, 0.1, 3)
    dw1 = 0.0                                                  # the restriction
    pre1, pre0 = 0.3 + zeta1 + b0 + b1 * X, 0.3 + zeta0 + w0 + w1 * X
    post1 = pre1 + db0 + db1 * X
    post0 = pre0 + dw0 + dw1 * X
    dT = X * (post1 - pre1) + (1 - X) * (post0 - pre0)
    pw = N / N.sum(); EX = np.sum(pw * X)
    DiD = (np.sum(pw * X * (post1 - pre1)) / EX) - (np.sum(pw * (1 - X) * (post0 - pre0)) / (1 - EX))
    # regression of dT on 1, X, X^2 (weights irrelevant: exact fit)
    Z = np.column_stack([np.ones(pn), X, X ** 2])
    a0, c1, d1 = np.linalg.lstsq(Z, dT, rcond=None)[0]
    # with dw1 = 0: dw0 = a0, db0 = a0 + c1, db1 = d1
    est = (np.sum(pw * X * (a0 + c1 + d1 * X)) / EX) - a0
    worst_b = max(worst_b, abs(est - DiD))
check("V7b first-difference LCE + (Delta w1 = 0) point-identifies the DiD despite arbitrary unit effects",
      worst_b < 1e-10, f"max err {worst_b:.1e}")

worst_c = 0.0
for _ in range(200):
    N, X, b1a, b0a = random_pop()
    _, _, b1b, b0b = random_pop()
    o = mei_objects(N, X, b1b - b1a, b0b - b0a)      # outcome = within-person change
    worst_c = max(worst_c, abs(o["ER"] - o["D"] - o["dW"] / o["VarXN"]))
check("V7c ER on Delta T is biased for the DiD by dW(changes)/Var(X_N)", worst_c < 1e-10,
      f"max err {worst_c:.1e}")

# ---------------------------------------------------------------------------
# V8. M1 across periods (common slopes, row-specific baseline shifts zeta_r):
#     Delta w1 = Delta b1 = 0 (R1 and R2) but Delta kappa = (zeta_v - zeta_n)/2.
# ---------------------------------------------------------------------------
worst = 0.0
for _ in range(200):
    H = rng.uniform(0, 1, (2, 2))            # rows/cols ordered (v, n); H[r, s]
    zeta = rng.uniform(-0.2, 0.2, 2)
    Hp = H + zeta[:, None]                   # H + zeta 1'
    dH = Hp - H
    dw1 = dH[1, 0] - dH[1, 1]                # non-veteran slope change
    db1 = dH[0, 0] - dH[0, 1]                # veteran slope change
    dkap = (dH[0, 1] - dH[1, 0]) / 2
    worst = max(worst, abs(dw1), abs(db1), abs(dkap - (zeta[0] - zeta[1]) / 2))
check("V8 M1 across periods = R1 and R2; Delta kappa = (zeta_v - zeta_n)/2, not 0", worst < 1e-14,
      f"max err {worst:.1e}")

print(f"\n{sum(results)}/{len(results)} PASS")
