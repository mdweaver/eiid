"""Checks for PHASE2d_did.md (difference-in-differences with partial identification).

Run:  python3 -I verify_phase2d.py > verify_phase2d_output.txt
Rows: v (veteran, index 0), n (non-veteran, index 1). Columns: yes (0), no (1), abstain (2).
"""
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(1857_1865)
results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))


def wcov(a, b, w):
    w = w / w.sum()
    return np.sum(w * a * b) - np.sum(w * a) * np.sum(w * b)


def moments(pw, x):
    return dict(mvv=np.sum(pw * x * x), mvn=np.sum(pw * x * (1 - x)), mnn=np.sum(pw * (1 - x) ** 2),
                mv=np.sum(pw * x), mn=np.sum(pw * (1 - x)))


def fitQ(pw, x, t):
    """Weighted LS of t on x^2, 2x(1-x), (1-x)^2 -> (Qvv, Qvn, Qnn) (homogeneous form, exact for quadratics)."""
    Z = np.column_stack([x ** 2, 2 * x * (1 - x), (1 - x) ** 2])
    W = np.sqrt(pw)
    return np.linalg.lstsq(Z * W[:, None], t * W, rcond=None)[0]


def did_of_dkappa(Q, mom, dk):
    """Model-implied DiD for one column as an affine function of Delta kappa (2a Lemma 2 applied to Delta H)."""
    qvv, qvn, qnn = Q
    dBv = (qvv * mom["mvv"] + (qvn + dk) * mom["mvn"]) / mom["mv"]
    dBn = (qnn * mom["mnn"] + (qvn - dk) * mom["mvn"]) / mom["mn"]
    return dBv - dBn


# ---------------------------------------------------------------------------
# D1. First-difference LCE with arbitrary (non-linear) unit levels.
#   For column yes: R1 (Delta w1 = 0), R2 (Delta b1 = 0), R4 (a Delta b1 + b Delta w1 = 0)
#   each point-identify the DiD; DiD_R4 lies between DiD_R1 and DiD_R2;
#   the sensitivity of the DiD to Delta w1 is -s, s = mu_vn (1/mu_v + 1/mu_n).
# ---------------------------------------------------------------------------
def fd_population(rule, p=120):
    N = rng.integers(50, 2000, p).astype(float); pw = N / N.sum()
    x = rng.uniform(0.05, 0.6, p)
    # arbitrary levels in the simplex (not linear in x), interior so that changes stay feasible
    lev = np.stack([np.stack([0.15 + 0.1 * np.sin(5 * x + k) ** 2 for k in range(3)], 1) for _ in range(2)], 1)
    lev = lev / lev.sum(2, keepdims=True)              # lev[i, r, c]
    # Delta H for column yes in homogeneous form; entries H[r, s] (s = composition coordinate v, n)
    dHvn, dHnn = rng.uniform(-0.05, 0.05, 2)
    db1, dw1 = rng.uniform(-0.05, 0.05, 2)
    mom = moments(pw, x)
    a = mom["mv"] * wcov(x, x, pw * x)
    b = mom["mn"] * wcov(x, x, pw * (1 - x))
    if rule == "R1":
        dw1 = 0.0
    elif rule == "R2":
        db1 = 0.0
    elif rule == "R4":
        dw1 = -a * db1 / b
    dH = np.zeros((3, 2, 2))
    dH[0] = [[dHvn + db1, dHvn], [dHnn + dw1, dHnn]]
    dH[1] = rng.uniform(-0.03, 0.03, (2, 2))
    dH[2] = -dH[0] - dH[1]
    X = np.column_stack([x, 1 - x])
    dbeta = np.einsum("crs,is->irc", dH, X)
    b0 = lev
    b1 = lev + dbeta
    assert (b1 > 0).all()
    T0 = np.einsum("ir,irc->ic", X, b0); T1 = np.einsum("ir,irc->ic", X, b1)
    wv, wn = pw * x / mom["mv"], pw * (1 - x) / mom["mn"]
    truth = np.sum(wv * dbeta[:, 0, 0]) - np.sum(wn * dbeta[:, 1, 0])
    return dict(pw=pw, x=x, dT=T1[:, 0] - T0[:, 0], mom=mom, truth=truth, a=a, b=b, T0=T0, T1=T1,
                b0=b0, b1=b1, dH=dH)


def estimates(P):
    Q = fitQ(P["pw"], P["x"], P["dT"])
    qvv, qvn, qnn = Q
    dk = {"R1": qvn - qnn, "R2": qvv - qvn}
    dk["R4"] = (P["a"] * dk["R2"] + P["b"] * dk["R1"]) / (P["a"] + P["b"])
    return Q, {k: did_of_dkappa(Q, P["mom"], v) for k, v in dk.items()}, dk


worst, between_ok, er_ok, s_ok = 0.0, True, True, True
for rule in ("R1", "R2", "R4"):
    for _ in range(100):
        P = fd_population(rule)
        Q, est, dk = estimates(P)
        worst = max(worst, abs(est[rule] - P["truth"]))
        lo, hi = sorted([est["R1"], est["R2"]])
        between_ok &= lo - 1e-12 <= est["R4"] <= hi + 1e-12
        if rule == "R4":     # Goodman (with intercept) on the change equals the DiD
            ER = wcov(P["dT"], P["x"], P["pw"]) / wcov(P["x"], P["x"], P["pw"])
            er_ok &= abs(ER - P["truth"]) < 1e-10
        s = P["mom"]["mvn"] * (1 / P["mom"]["mv"] + 1 / P["mom"]["mn"])
        # DiD(Delta w1) = DiD(Delta kappa_R1 - Delta w1): slope -s in Delta w1
        d1 = did_of_dkappa(Q, P["mom"], dk["R1"] - 0.01) - did_of_dkappa(Q, P["mom"], dk["R1"])
        s_ok &= abs(d1 + 0.01 * s) < 1e-12
check("D1 Prop FD2: R1, R2 and R4 each point-identify the DiD under FD-LCE with non-linear levels",
      worst < 1e-10, f"max err {worst:.1e}")
check("D1 Prop FD2(c): DiD_R4 lies between DiD_R1 and DiD_R2", between_ok)
check("D1 Prop FD2(d): under R4, Goodman on the change equals the DiD", er_ok)
check("D1 Prop FD3: dDiD/dDelta_w1 = -s, s = mu_vn (1/mu_v + 1/mu_n)", s_ok)


# ---------------------------------------------------------------------------
# D2. FD outer set via an LP with auxiliary level variables at the support points:
#   (a) the truth is feasible and the DiD bound contains the true DiD;
#   (b) when levels are linear in each period (per-period LCE), the product of
#       the per-period outer sets gives a DiD bound inside the FD outer bound.
# ---------------------------------------------------------------------------
def fd_outer_bounds(x, tau0, tau1, Q3, mom, C=3):
    """Variables: dk (C-1), m0[j, r, c]. Q3[c] = (qvv, qvn, qnn) of Delta T_c, c = 0..C-1."""
    J = len(x)
    nk = C - 1
    nv = nk + J * 2 * C
    im = lambda j, r, c: nk + (j * 2 + r) * C + c
    Aeq, beq, Aub, bub = [], [], [], []
    def dm_coef(j, r, c):
        """M_DeltaH(x_j)[r, c] = const + coef . dk."""
        qvv, qvn, qnn = Q3[c]
        xv, xn = x[j], 1 - x[j]
        coef = np.zeros(nk)
        if r == 0:
            const = qvv * xv + qvn * xn
            sgn = 1.0
        else:
            const = qvn * xv + qnn * xn
            sgn = -1.0
        k_coef = sgn * (xn if r == 0 else xv)          # K_vn multiplies x_n for row v; K_nv = -K_vn multiplies x_v
        if c < C - 1:
            coef[c] = k_coef
        else:
            coef[:] = -k_coef
        return const, coef
    for j in range(J):
        for r in range(2):
            a = np.zeros(nv); a[[im(j, r, c) for c in range(C)]] = 1; Aeq.append(a); beq.append(1)
        for c in range(C - 1):
            a = np.zeros(nv); a[im(j, 0, c)] = x[j]; a[im(j, 1, c)] = 1 - x[j]; Aeq.append(a); beq.append(tau0[j, c])
        for r in range(2):
            for c in range(C):           # 0 <= m0 + dm <= 1
                const, coef = dm_coef(j, r, c)
                a = np.zeros(nv); a[im(j, r, c)] = -1; a[:nk] = -coef; Aub.append(a); bub.append(const)
                a = np.zeros(nv); a[im(j, r, c)] = 1; a[:nk] = coef; Aub.append(a); bub.append(1 - const)
    # target: DiD for column 0 = did_of_dkappa(Q3[0], mom, dk0) = const + s * dk0
    s = mom["mvn"] * (1 / mom["mv"] + 1 / mom["mn"])
    const = did_of_dkappa(Q3[0], mom, 0.0)
    obj = np.zeros(nv); obj[0] = s
    kw = dict(A_eq=np.array(Aeq), b_eq=np.array(beq), A_ub=np.array(Aub), b_ub=np.array(bub),
              bounds=[(None, None)] * nk + [(0, 1)] * (J * 2 * C), method="highs")
    lo, hi = linprog(obj, **kw), linprog(-obj, **kw)
    return (const + lo.fun, const - hi.fun) if lo.status == 0 and hi.status == 0 else None


inside_ok, cnt = True, 0
for _ in range(40):
    P = fd_population("none", p=40)
    Q3 = [fitQ(P["pw"], P["x"], P["T1"][:, c] - P["T0"][:, c]) for c in range(3)]
    b = fd_outer_bounds(P["x"], P["T0"], P["T1"], Q3, P["mom"])
    cnt += b is not None
    inside_ok &= b is not None and b[0] - 1e-9 <= P["truth"] <= b[1] + 1e-9
check("D2(a) FD outer LP contains the true DiD (levels non-linear, exact tau)", inside_ok and cnt == 40)


def lce_period_interval(pw, x, T, mom, C=3):
    """Per-period LCE outer set (2a Prop 6, vertices = observed range) -> interval of kappa_0 (column 0)
    and the map to B_v(col 0), B_n(col 0)."""
    Q3 = [fitQ(pw, x, T[:, c]) for c in range(C)]
    verts = [x.min(), x.max()]
    Aub, bub = [], []
    for c in range(C):
        qvv, qvn, qnn = Q3[c]
        for xv in verts:
            xn = 1 - xv
            for r in range(2):
                const = (qvv * xv + qvn * xn) if r == 0 else (qvn * xv + qnn * xn)
                kc = (xn if r == 0 else -xv)
                coef = np.zeros(C - 1)
                if c < C - 1:
                    coef[c] = kc
                else:
                    coef[:] = -kc
                Aub.append(-coef); bub.append(const)
    kw = dict(A_ub=np.array(Aub), b_ub=np.array(bub), bounds=[(None, None)] * (C - 1), method="highs")
    e = np.zeros(C - 1); e[0] = 1
    lo, hi = linprog(e, **kw), linprog(-e, **kw)
    return Q3, (lo.fun, -hi.fun)


nest_ok, n_nest = True, 0
for _ in range(60):
    p = 40
    N = rng.integers(50, 2000, p).astype(float); pw = N / N.sum()
    x = rng.uniform(0.05, 0.6, p); X = np.column_stack([x, 1 - x])
    H0 = rng.dirichlet([2, 2, 2], (2, 2)).transpose(2, 0, 1)
    H1 = np.clip(H0 + rng.normal(0, 0.04, H0.shape), 0.01, None); H1 = H1 / H1.sum(0, keepdims=True)
    b0 = np.einsum("crs,is->irc", H0, X); b1 = np.einsum("crs,is->irc", H1, X)
    T0 = np.einsum("ir,irc->ic", X, b0); T1 = np.einsum("ir,irc->ic", X, b1)
    mom = moments(pw, x)
    Qa, ka = lce_period_interval(pw, x, T0, mom)
    Qb, kb = lce_period_interval(pw, x, T1, mom)
    # product bound for the DiD of column 0: DiD = [B_v1 - B_n1] - [B_v0 - B_n0], each affine in its kappa
    s = mom["mvn"] * (1 / mom["mv"] + 1 / mom["mn"])
    gap = lambda Q, k: did_of_dkappa(Q, mom, k)    # level gap B_v - B_n with kappa k
    prod = (gap(Qb[0], kb[0]) - gap(Qa[0], ka[1]), gap(Qb[0], kb[1]) - gap(Qa[0], ka[0]))
    Q3 = [fitQ(pw, x, T1[:, c] - T0[:, c]) for c in range(3)]
    fd = fd_outer_bounds(x, T0, T1, Q3, mom)
    if fd is None:
        continue
    n_nest += 1
    nest_ok &= fd[0] <= prod[0] + 1e-8 and prod[1] <= fd[1] + 1e-8
check("D2(b) per-period LCE (product of outer sets) is inside the FD outer bound", nest_ok and n_nest > 40,
      f"{n_nest} cases")


# ---------------------------------------------------------------------------
# D3. Monotone EI on changes (fixed population, exact): CR on changes is never
#     empty, has one branch, and sign(DiD) = sign(ER_Delta) when realised CR holds.
# ---------------------------------------------------------------------------
def dd_lp_did(pw, x, T0, T1, extra_ge=()):
    """Variables beta_v,i,0 and beta_v,i,1 (2x2 per period, column yes); beta_n implied."""
    p = len(x)
    L = lambda T: np.maximum(0, (T - (1 - x)) / x)
    U = lambda T: np.minimum(1, T / x)
    mv, mn = np.sum(pw * x), np.sum(pw * (1 - x))
    cv = pw * x * (1 / mv + 1 / mn)                 # B_v - B_n = cv . beta_v - sum pw T / mn
    obj = np.concatenate([-cv, cv])                  # DiD = (gap1 - gap0)
    const = -np.sum(pw * T1) / mn + np.sum(pw * T0) / mn
    Aub, bub = [], []
    for a, k in extra_ge:
        Aub.append(-a); bub.append(-k)
    kw = dict(bounds=list(zip(L(T0), U(T0))) + list(zip(L(T1), U(T1))), method="highs")
    if Aub:
        kw.update(A_ub=np.array(Aub), b_ub=np.array(bub))
    lo, hi = linprog(obj, **kw), linprog(-obj, **kw)
    return (const + lo.fun, const - hi.fun) if lo.status == 0 and hi.status == 0 else None


ok_one, ok_sign, n_cr = True, True, 0
for _ in range(200):
    p = 25
    N = rng.integers(50, 2000, p).astype(float); pw = N / N.sum()
    x = rng.uniform(0.05, 0.9, p)
    bv0, bn0, bv1, bn1 = rng.uniform(0, 1, (4, p))
    T0, T1 = x * bv0 + (1 - x) * bn0, x * bv1 + (1 - x) * bn1
    mv, mn = np.sum(pw * x), np.sum(pw * (1 - x))
    dT = T1 - T0
    ER = wcov(dT, x, pw) / wcov(x, x, pw)
    NM = np.sum(pw * x * dT) / mv - np.sum(pw * (1 - x) * dT) / mn
    # DiD <= ER_D  <=> dW_D >= 0 ;  DiD >= NM_D <=> dB_D >= 0
    cv = pw * x * (1 / mv + 1 / mn); const = -np.sum(pw * dT) / mn
    a = np.concatenate([-cv, cv])                    # DiD = a.beta + const
    plus = dd_lp_did(pw, x, T0, T1, [(-a, const - ER), (a, NM - const)])
    minus = dd_lp_did(pw, x, T0, T1, [(a, ER - const), (-a, const - NM)])
    branches = [b for b in (plus, minus) if b is not None]
    ok_one &= len(branches) == 1 or abs(ER) < 1e-12
    truth = np.sum(pw * x * (bv1 - bv0)) / mv - np.sum(pw * (1 - x) * (bn1 - bn0)) / mn
    dW = (ER - truth) * wcov(x, x, pw); dB = (truth - NM) * mv * mn
    if dW * dB >= 0:
        n_cr += 1
        ok_sign &= np.sign(truth) == np.sign(ER)
check("D3 Prop D4: CR on changes never empty, one branch", ok_one)
check("D3 Prop D4: under realised CR on changes, sign(DiD) = sign(ER on Delta T)", ok_sign, f"{n_cr} cases")


# ---------------------------------------------------------------------------
# D4. Changing population counts (levels LCE per period, R1 across periods):
#     DiD = I' + (s1 - s0) kappa_0 exactly; truth in the interval.
# ---------------------------------------------------------------------------
ok_formula, ok_inside = True, True
for _ in range(100):
    p = 40
    x = rng.uniform(0.05, 0.6, p); X = np.column_stack([x, 1 - x])
    N0 = rng.integers(50, 2000, p).astype(float)
    N1 = N0 * rng.uniform(0.7, 1.6, p)                      # counts change, unevenly
    pw0, pw1 = N0 / N0.sum(), N1 / N1.sum()
    H0 = rng.dirichlet([2, 2, 2], (2, 2)).transpose(2, 0, 1)
    dH = rng.normal(0, 0.03, (3, 2, 2)); dH[0, 1, 0] = dH[0, 1, 1]      # R1 in column yes
    dH[2] = -dH[0] - dH[1]
    H1 = H0 + dH
    if (H1 < 0).any():
        continue
    m0, m1 = np.einsum("crs,is->irc", H0, X), np.einsum("crs,is->irc", H1, X)
    T0, T1 = np.einsum("ir,irc->ic", X, m0), np.einsum("ir,irc->ic", X, m1)
    mo0, mo1 = moments(pw0, x), moments(pw1, x)
    truth = (np.sum(pw1 * x * m1[:, 0, 0]) / mo1["mv"] - np.sum(pw1 * (1 - x) * m1[:, 1, 0]) / mo1["mn"]) - \
            (np.sum(pw0 * x * m0[:, 0, 0]) / mo0["mv"] - np.sum(pw0 * (1 - x) * m0[:, 1, 0]) / mo0["mn"])
    Qa, ka = lce_period_interval(pw0, x, T0, mo0)
    Qb, kb = lce_period_interval(pw1, x, T1, mo1)
    # R1: kappa_1 = kappa_0 + (dQvn - dQnn) with dQ = Q1 - Q0 (column yes)
    shift = (Qb[0][1] - Qa[0][1]) - (Qb[0][2] - Qa[0][2])
    f = lambda k0: did_of_dkappa(Qb[0], mo1, k0 + shift) - did_of_dkappa(Qa[0], mo0, k0)
    k0_true = 0.5 * (H0[0, 0, 1] - H0[0, 1, 0])
    ok_formula &= abs(f(k0_true) - truth) < 1e-10
    s0 = mo0["mvn"] * (1 / mo0["mv"] + 1 / mo0["mn"]); s1 = mo1["mvn"] * (1 / mo1["mv"] + 1 / mo1["mn"])
    ok_formula &= abs((f(0.1) - f(0.0)) - 0.1 * (s1 - s0)) < 1e-10
    lo_k, hi_k = max(ka[0], kb[0] - shift), min(ka[1], kb[1] - shift)
    vals = sorted([f(lo_k), f(hi_k)])
    ok_inside &= lo_k <= hi_k + 1e-12 and vals[0] - 1e-10 <= truth <= vals[1] + 1e-10
check("D4 Prop D6: with changing counts, DiD = I' + (s1 - s0) kappa_0 under R1", ok_formula)
check("D4 Prop D6: truth inside the resulting interval", ok_inside)


# ---------------------------------------------------------------------------
# D5. Aging-in decomposition (post-war base for both periods).
# ---------------------------------------------------------------------------
ok = True
for _ in range(300):
    p = 30
    O = rng.integers(20, 500, (p, 2)).astype(float)      # eligible in 1857, rows v, n
    A = rng.integers(0, 200, (p, 2)).astype(float)       # aged in
    A[:, 0] = A[:, 0] * 1.5                              # veterans younger
    y0 = rng.uniform(0, 0.15, (p, 2)) * O                 # 1857 yes among O
    y1O = rng.uniform(0, 0.8, (p, 2)) * O                 # 1865 yes among O
    y1A = rng.uniform(0, 0.8, (p, 2)) * A
    om = O.sum(0) / (O + A).sum(0)
    BO0, BO1 = y0.sum(0) / O.sum(0), y1O.sum(0) / O.sum(0)
    BA1 = y1A.sum(0) / A.sum(0)
    Borig0 = y0.sum(0) / (O + A).sum(0)
    B1 = (y1O + y1A).sum(0) / (O + A).sum(0)
    did_orig = (B1[0] - Borig0[0]) - (B1[1] - Borig0[1])
    did_O = (BO1[0] - BO0[0]) - (BO1[1] - BO0[1])
    rhs = (om[0] * (BO1[0] - BO0[0]) + (1 - om[0]) * BA1[0]) - (om[1] * (BO1[1] - BO0[1]) + (1 - om[1]) * BA1[1])
    ok &= abs(did_orig - rhs) < 1e-12
    # under B^A_1 = B^O_1: did_orig - did_O = (1-om_v) BO0_v - (1-om_n) BO0_n
    B1_eq = BO1
    did_orig_eq = (B1_eq[0] - Borig0[0]) - (B1_eq[1] - Borig0[1])
    ok &= abs((did_orig_eq - did_O) - ((1 - om[0]) * BO0[0] - (1 - om[1]) * BO0[1])) < 1e-12
check("D5 Prop C2: aging-in decomposition of the post-war-base DiD; bias formula when aged-in vote like stayers", ok)


# ---------------------------------------------------------------------------
# D6. Composition polytope (aged-in split interval-valued, stayers' DiD):
#     model-free bound by LP for fixed V_O, scanned over V_O, contains the truth.
# ---------------------------------------------------------------------------
def composition_bounds(O, A, V, Y0, Y1, grid):
    """Per unit: O (eligible 1857), A (aged in), V (veterans, all present). Unknown split Ov in
    [max(0, V-A), min(V, O)]. Cells: y0v, y0n (1857 yes among stayers), y1vO, y1nO, y1vA, y1nA."""
    p = len(O)
    k = 7                                        # Ov, y0v, y0n, y1vO, y1nO, y1vA, y1nA
    nv = p * k
    I = lambda i, j: i * k + j
    Aeq, beq, Aub, bub, bnds = [], [], [], [], []
    for i in range(p):
        bnds += [(max(0, V[i] - A[i]), min(V[i], O[i]))] + [(0, None)] * 6
        a = np.zeros(nv); a[[I(i, 1), I(i, 2)]] = 1; Aeq.append(a); beq.append(Y0[i])
        a = np.zeros(nv); a[[I(i, 3), I(i, 4), I(i, 5), I(i, 6)]] = 1; Aeq.append(a); beq.append(Y1[i])
        rows = [(1, +1, 0), (2, -1, O[i]), (3, +1, 0), (4, -1, O[i]), (5, -1, V[i]), (6, +1, A[i] - V[i])]
        # y0v <= Ov ; y0n <= O - Ov ; y1vO <= Ov ; y1nO <= O - Ov ; y1vA <= V - Ov ; y1nA <= A - (V - Ov)
        for j, sgn, cst in rows:
            a = np.zeros(nv); a[I(i, j)] = 1
            if j in (1, 3):
                a[I(i, 0)] = -1; Aub.append(a); bub.append(0)
            elif j in (2, 4):
                a[I(i, 0)] = 1; Aub.append(a); bub.append(O[i])
            elif j == 5:
                a[I(i, 0)] = 1; Aub.append(a); bub.append(V[i])
            else:
                a[I(i, 0)] = -1; Aub.append(a); bub.append(A[i] - V[i])
    lo_all, hi_all = np.inf, -np.inf
    Otot = O.sum()
    for VO in grid:
        a = np.zeros(nv); a[[I(i, 0) for i in range(p)]] = 1
        Aeq2, beq2 = Aeq + [a], beq + [VO]
        obj = np.zeros(nv)
        for i in range(p):
            obj[I(i, 3)] += 1 / VO; obj[I(i, 1)] -= 1 / VO
            obj[I(i, 4)] -= 1 / (Otot - VO); obj[I(i, 2)] += 1 / (Otot - VO)
        kw = dict(A_eq=np.array(Aeq2), b_eq=np.array(beq2), A_ub=np.array(Aub), b_ub=np.array(bub),
                  bounds=bnds, method="highs")
        lo, hi = linprog(obj, **kw), linprog(-obj, **kw)
        if lo.status == 0:
            lo_all = min(lo_all, lo.fun); hi_all = max(hi_all, -hi.fun)
    return lo_all, hi_all


ok = True
for _ in range(15):
    p = 12
    O = rng.integers(50, 300, p).astype(float); A = rng.integers(20, 120, p).astype(float)
    Ov = np.floor(rng.uniform(0.1, 0.3, p) * O); Av = np.floor(rng.uniform(0.3, 0.6, p) * A)
    V = Ov + Av
    y0v = np.floor(rng.uniform(0, .1, p) * Ov); y0n = np.floor(rng.uniform(0, .1, p) * (O - Ov))
    y1vO = np.floor(rng.uniform(.2, .6, p) * Ov); y1nO = np.floor(rng.uniform(.05, .3, p) * (O - Ov))
    y1vA = np.floor(rng.uniform(.2, .6, p) * Av); y1nA = np.floor(rng.uniform(.05, .3, p) * (A - Av))
    truth = (y1vO.sum() - y0v.sum()) / Ov.sum() - (y1nO.sum() - y0n.sum()) / (O - Ov).sum()
    lo_V = np.sum(np.maximum(0, V - A)); hi_V = np.sum(np.minimum(V, O))
    grid = np.unique(np.concatenate([np.linspace(lo_V, hi_V, 25), [Ov.sum()]]))
    grid = grid[(grid > 0) & (grid < O.sum())]
    lo, hi = composition_bounds(O, A, V, y0v + y0n, y1vO + y1nO + y1vA + y1nA, grid)
    ok &= lo - 1e-9 <= truth <= hi + 1e-9
check("D6 Prop C1: composition polytope (interval aged-in split) bound contains the stayers' DiD", ok)


# ---------------------------------------------------------------------------
# D7. Causal layer: ATT bounds under MTR, no-decline and positive selection.
#     Epigraph LP = closed form B_v1^U - max(B_v0^L, B_n1 at B_v1^U) (2 rows, product).
# ---------------------------------------------------------------------------
ok = True
for _ in range(200):
    p = 20
    N = rng.integers(50, 2000, p).astype(float); pw = N / N.sum()
    x = rng.uniform(0.05, 0.9, p)
    T0 = x * rng.uniform(0, .3, p) + (1 - x) * rng.uniform(0, .3, p)
    T1 = x * rng.uniform(.2, 1, p) + (1 - x) * rng.uniform(0, .6, p)
    mv, mn = np.sum(pw * x), np.sum(pw * (1 - x))
    L = lambda T: np.maximum(0, (T - (1 - x)) / x); U = lambda T: np.minimum(1, T / x)
    cv = pw * x / mv
    # variables: beta_v0 (p), beta_v1 (p), t
    Bv0 = np.concatenate([cv, np.zeros(p), [0]])
    Bv1 = np.concatenate([np.zeros(p), cv, [0]])
    # B_n1 = (sum pw T1 - mv B_v1)/mn
    Bn1_c, Bn1 = np.sum(pw * T1) / mn, -np.concatenate([np.zeros(p), pw * x / mn, [0]])
    tvec = np.zeros(2 * p + 1); tvec[-1] = 1
    Aub = [tvec - (Bv1 - Bv0), tvec - (Bv1 - Bn1)]
    bub = [0, -Bn1_c]
    bnds = list(zip(L(T0), U(T0))) + list(zip(L(T1), U(T1))) + [(None, None)]
    r = linprog(-tvec, A_ub=np.array(Aub), b_ub=np.array(bub), bounds=bnds, method="highs")
    bv1U = cv @ U(T1); bv0L = cv @ L(T0); bn1_at = Bn1_c - (pw * x / mn) @ U(T1)
    cf = bv1U - max(bv0L, bn1_at)
    ok &= r.status == 0 and abs(-r.fun - cf) < 1e-9
check("D7 Prop CL: ATT upper bound (epigraph LP) = closed form in the 2-row product case", ok)


# ---------------------------------------------------------------------------
# D8. Party information.
#   P2: veteran x party x suffrage with zero (Democrat, yes) gives the same
#       veteran-yes range as Duncan-Davis (coupling).
#   P3: selection on party only (pre-war) gives vet-yes = v_R * Y / R, a range
#       inside Duncan-Davis that contains the truth.
# ---------------------------------------------------------------------------
ok2, ok3, strict3 = True, True, 0
for _ in range(300):
    Nn = rng.integers(100, 1000)
    R_, D_ = rng.integers(10, Nn // 2), rng.integers(10, Nn // 3)
    Ab = Nn - R_ - D_                                         # did not vote for governor
    V = rng.integers(5, Nn // 2)
    # suffrage: yes only from R (and Abstain-governor gives no yes here), no from anyone
    piR = rng.uniform(0, 0.3)
    vR = rng.integers(max(0, V - (Nn - R_)), min(V, R_) + 1)
    Y = int(round(piR * R_))
    # P2: LP over 3-way counts n[g, party, s], g in {v, n}, party in {R, D, Ab}, s in {yes, other}
    idx = lambda g, q, s: (g * 3 + q) * 2 + s
    Aeq, beq = [], []
    for g, tot in ((0, V), (1, Nn - V)):
        a = np.zeros(12); a[[idx(g, q, s) for q in range(3) for s in range(2)]] = 1; Aeq.append(a); beq.append(tot)
    for q, tot in ((0, R_), (1, D_), (2, Ab)):
        a = np.zeros(12); a[[idx(g, q, s) for g in range(2) for s in range(2)]] = 1; Aeq.append(a); beq.append(tot)
    a = np.zeros(12); a[[idx(g, q, 0) for g in range(2) for q in range(3)]] = 1; Aeq.append(a); beq.append(Y)
    bnds = [(0, None)] * 12
    for g in range(2):
        bnds[idx(g, 1, 0)] = (0, 0)                            # Democrats never yes
    tgt = np.zeros(12); tgt[[idx(0, q, 0) for q in range(3)]] = 1
    lo = linprog(tgt, A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=bnds, method="highs")
    hi = linprog(-tgt, A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=bnds, method="highs")
    dd = (max(0, V + Y - Nn), min(V, Y))
    ok2 &= abs(lo.fun - dd[0]) < 1e-7 and abs(-hi.fun - dd[1]) < 1e-7
    # P3 (yes only from R, same rate for v and n within R): vet-yes = v_R * Y / R
    lo3 = max(0, V + R_ - Nn) * Y / R_; hi3 = min(V, R_) * Y / R_
    truth = vR * Y / R_
    ok3 &= dd[0] - 1e-9 <= lo3 <= truth + 1e-9 and truth - 1e-9 <= hi3 <= dd[1] + 1e-9
    strict3 += (hi3 - lo3) < (dd[1] - dd[0]) - 1e-9
check("D8 Prop P2: party margin + (Democrat, yes) = 0 leaves the veteran-yes range at Duncan-Davis", ok2)
check("D8 Prop P3: selection on party only gives a range inside Duncan-Davis containing the truth", ok3,
      f"strictly narrower in {strict3}/300")

print(f"\n{sum(results)}/{len(results)} PASS")
