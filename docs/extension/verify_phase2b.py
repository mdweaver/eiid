# Numerical checks of the Phase 2b results (PHASE2b_districts.md).
# Units nested in districts; models M0 (pooled), M2 (district-specific),
# M1(delta) (bounded deviation), M1 (shared contextual effects + district
# intercepts), L5 (constancy within district).
# Each check prints PASS/FAIL. Run: python3 -I verify_phase2b.py
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

rng = np.random.default_rng(20261029)
ok = lambda name, cond: print(("PASS " if cond else "FAIL ") + name)
TOL = 1e-7


# ---------------- parametrisation (as in Phase 2a, Prop. 7) ----------------
def pairs(R):
    return [(r, s) for r in range(R) for s in range(r + 1, R)]


def K_of(k, R, C):
    """free coordinates k (length (C-1)R(R-1)/2) -> list of C skew matrices summing to 0"""
    P, Ks, i = pairs(R), [], 0
    for c in range(C - 1):
        K = np.zeros((R, R))
        for r, s in P:
            K[r, s], K[s, r] = k[i], -k[i]; i += 1
        Ks.append(K)
    Ks.append(-sum(Ks))
    return Ks


def k_of(Ks):
    """inverse of K_of (reads columns 1..C-1)"""
    R = Ks[0].shape[0]
    return np.array([Ks[c][r, s] for c in range(len(Ks) - 1) for r, s in pairs(R)])


def basis(R, C):
    d0 = (C - 1) * len(pairs(R))
    return [K_of(np.eye(d0)[j], R, C) for j in range(d0)], d0


def sym(A): return (A + A.T) / 2
def skew(A): return (A - A.T) / 2


# ---------------- data-generating processes ----------------
def hull_vertices(X):
    R = X.shape[1]
    if R == 2:
        return [X[X[:, 0].argmin()], X[X[:, 0].argmax()]]
    return list(X[ConvexHull(X[:, :R - 1]).vertices])


def quad_fit(X, T):
    """OLS of each column of T on the homogeneous quadratic monomials of X -> Q_c"""
    R, C = X.shape[1], T.shape[1]
    idx = [(r, s) for r in range(R) for s in range(r, R)]
    Z = np.column_stack([X[:, r] * X[:, s] for r, s in idx])
    coef = np.linalg.lstsq(Z, T, rcond=None)[0]
    Qs = []
    for c in range(C):
        Q = np.zeros((R, R))
        for j, (r, s) in enumerate(idx):
            if r == s: Q[r, r] = coef[j, c]
            else: Q[r, s] = Q[s, r] = coef[j, c] / 2
        Qs.append(Q)
    return Qs


def make_districts(H, zeta, alphas, polls=250, N_rng=(300, 1500)):
    """H: list of C (RxR) common coefficient matrices; zeta[d]: list of C vectors (len R).
    Polls in district d: Dirichlet(alphas[d]). Noise-free T = M(X)'X (identification only)."""
    R, C = H[0].shape[0], len(H)
    out = []
    for d, a in enumerate(alphas):
        X = rng.dirichlet(a, polls)
        Hd = [H[c] + np.outer(zeta[d][c], np.ones(R)) for c in range(C)]
        M = np.stack([np.stack([X @ Hd[c][r] for c in range(C)], 1) for r in range(R)], 1)  # p x R x C
        assert M.min() >= -1e-12, "invalid district rates"
        T = np.einsum("ir,irc->ic", X, M)
        N = rng.integers(*N_rng, polls).astype(float)
        out.append(dict(X=X, T=T, N=N, M=M, Hd=Hd))
    return out


def moments(dist):
    R = dist[0]["X"].shape[1]
    mu_r = sum((d["N"][:, None] * d["X"]).sum(0) for d in dist)
    mu_d = [np.einsum("i,ir,is->rs", d["N"], d["X"], d["X"]) for d in dist]
    return mu_r, mu_d


def truth(dist, r, c):
    mu_r, _ = moments(dist)
    return sum(np.sum(d["N"] * d["X"][:, r] * d["M"][:, r, c]) for d in dist) / mu_r[r]


def dd_bound(dist, r, c):
    mu_r, _ = moments(dist)
    lo = sum(np.sum(d["N"] * np.maximum(0, d["X"][:, r] + d["T"][:, c] - 1)) for d in dist)
    hi = sum(np.sum(d["N"] * np.minimum(d["X"][:, r], d["T"][:, c])) for d in dist)
    return lo / mu_r[r], hi / mu_r[r]


# ---------------- the district LP ----------------
def district_lp(Qd, Vd, mu_r, mu_d, r, c, mode, delta=0.0, zeta_c=None):
    """Bounds on B*_rc. Variables: kbar (d0) then k^(d) for each district (d0 each).
    mode: 'M2' (no link), 'M1' (k^(d) = kbar + c_d), 'M1d' (|k^(d) - kbar - c_d| <= delta)."""
    R, C = Qd[0][0].shape[0], len(Qd[0])
    Kb, d0 = basis(R, C)
    D = len(Qd)
    nv = d0 * (D + 1)
    A, b, Aeq, beq = [], [], [], []
    const = 0.0
    obj = np.zeros(nv)
    for d in range(D):
        off = d0 * (d + 1)
        for v in Vd[d]:
            for cc in range(C):
                for rr in range(R):
                    row = np.zeros(nv)
                    row[off:off + d0] = [-(Kb[j][cc] @ v)[rr] for j in range(d0)]
                    A.append(row); b.append((Qd[d][cc] @ v)[rr])
        const += (Qd[d][c][r] @ mu_d[d][r]) / mu_r[r]
        obj[off:off + d0] = [(Kb[j][c][r] @ mu_d[d][r]) / mu_r[r] for j in range(d0)]
        if mode in ("M1", "M1d"):
            cd = k_of([skew(np.outer(zeta_c[d][cc], np.ones(R))) for cc in range(C)])
            for j in range(d0):
                row = np.zeros(nv); row[off + j] = 1; row[j] = -1
                if mode == "M1":
                    Aeq.append(row); beq.append(cd[j])
                else:
                    A.append(row); b.append(cd[j] + delta)
                    A.append(-row); b.append(-cd[j] + delta)
    kw = dict(A_ub=np.array(A), b_ub=np.array(b), bounds=[(None, None)] * nv, method="highs")
    if Aeq:
        kw.update(A_eq=np.array(Aeq), b_eq=np.array(beq))
    lo, hi = linprog(obj, **kw), linprog(-obj, **kw)
    if lo.status != 0 or hi.status != 0:
        return None
    return const + lo.fun, const - hi.fun


def zeta_from_Q(Qd):
    """Lemma D1: zeta_{c,d} = diag(Q^(d)_c - Q^(1)_c)"""
    return [[np.diag(Qd[d][c] - Qd[0][c]).copy() for c in range(len(Qd[0]))] for d in range(len(Qd))]


def goodman_within(dist, r, c):
    """L5: regress T_c on X without intercept within each district; aggregate."""
    mu_r, _ = moments(dist)
    tot = 0.0
    for d in dist:
        alpha = np.linalg.lstsq(d["X"], d["T"], rcond=None)[0]   # R x C
        tot += np.sum(d["N"] * d["X"][:, r]) * alpha[r, c]
    return tot / mu_r[r]


# ================= checks =================
# common H (3x3), rows/cols (P, O, N); columns ordered (abstain, P, O) as in Phase 2a
def H_common():
    def M(x):
        xP, xO, xN = x
        return np.array([[0.08 - 0.03 * xP, 0.80 + 0.10 * xP, 0.12 - 0.07 * xP],
                         [0.10 + 0.02 * xO, 0.10 - 0.05 * xO, 0.80 + 0.03 * xO],
                         [0.75 - 0.10 * xN, 0.15 + 0.05 * xP, 0.10 + 0.10 * xN - 0.05 * xP]])
    return [np.array([[M(np.eye(3)[s])[r, c] for s in range(3)] for r in range(3)]) for c in range(3)]


def random_zeta(D, R, C, scale):
    Z = [[np.zeros(R) for _ in range(C)]]
    for _ in range(D - 1):
        z = rng.normal(0, scale, (C, R)); z -= z.mean(0)          # sum over columns = 0
        Z.append([z[c] for c in range(C)])
    return Z


H = H_common()
R, C, D = 3, 3, 6
alphas = [[6, 2, 2], [2, 6, 2], [2, 2, 6], [4, 4, 2], [3, 3, 3], [5, 2, 3]]
while True:
    zeta = random_zeta(D, R, C, 0.04)
    try:
        dist = make_districts(H, zeta, alphas); break
    except AssertionError:
        continue

Qd = [quad_fit(d["X"], d["T"]) for d in dist]
Vd = [hull_vertices(d["X"]) for d in dist]
mu_r, mu_d = moments(dist)

# Lemma D1: zeta identified from district quadratic fits; overidentifying restrictions hold
zhat = zeta_from_Q(Qd)
err_z = max(np.abs(zhat[d][c] - zeta[d][c]).max() for d in range(D) for c in range(C))
ok(f"Lemma D1: zeta recovered from diag(Q^(d)-Q^(1)) (max error {err_z:.1e})", err_z < 1e-8)
over = max(abs((Qd[d][c] - Qd[0][c])[r_, s_] - (zhat[d][c][r_] + zhat[d][c][s_]) / 2)
           for d in range(D) for c in range(C) for r_, s_ in pairs(R))
ok(f"Lemma D1: overidentifying restrictions hold under M1 (max violation {over:.1e})", over < 1e-8)
Ktrue = [skew(H[c]) for c in range(C)]
kd_true = [k_of([skew(dist[d]["Hd"][c]) for c in range(C)]) for d in range(D)]
cd = [k_of([skew(np.outer(zhat[d][c], np.ones(R))) for c in range(C)]) for d in range(D)]
link = max(np.abs(kd_true[d] - k_of(Ktrue) - cd[d]).max() for d in range(D))
ok(f"Lemma D1: K^(d) = K + skew(zeta_d 1') (max error {link:.1e})", link < 1e-8)

# ladder table for four targets
targets = [(0, 1, "loyalty P->P"), (0, 2, "defection P->O"), (0, 0, "demobilisation P->abstain"),
           (2, 1, "mobilisation N->P")]
rows = []
nest_ok = cover_ok = mono_ok = True
for r, c, lab in targets:
    tr = truth(dist, r, c)
    dd = dd_bound(dist, r, c)
    m2 = district_lp(Qd, Vd, mu_r, mu_d, r, c, "M2")
    md = [district_lp(Qd, Vd, mu_r, mu_d, r, c, "M1d", dl, zhat) for dl in (0.10, 0.05, 0.02)]
    m1 = district_lp(Qd, Vd, mu_r, mu_d, r, c, "M1", zeta_c=zhat)
    # pooled M0 (misspecified when zeta varies with composition)
    Xall = np.vstack([d["X"] for d in dist]); Tall = np.vstack([d["T"] for d in dist])
    Qp = quad_fit(Xall, Tall)
    m0 = district_lp([Qp], [hull_vertices(Xall)], mu_r, [sum(mu_d)], r, c, "M2")
    l5 = goodman_within(dist, r, c)
    rows.append((lab, tr, dd, m0, m2, md, m1, l5))
    chain = [m2] + md + [m1]
    for a, b in zip(chain, chain[1:]):
        nest_ok &= b[0] >= a[0] - TOL and b[1] <= a[1] + TOL
    for a, b in zip(md, md[1:]):
        mono_ok &= b[1] - b[0] <= a[1] - a[0] + TOL
    cover_ok &= all(x[0] - TOL <= tr <= x[1] + TOL for x in chain) and dd[0] - TOL <= tr <= dd[1] + TOL
    nest_ok &= m2[0] >= dd[0] - TOL and m2[1] <= dd[1] + TOL

fmt = lambda x: "empty" if x is None else f"[{x[0]:.3f}, {x[1]:.3f}]"
print("     3x3 ladder (6 districts, noise-free):")
print("     target                     truth  DD              M0 pooled       M2              M1(.10)         M1(.05)         M1(.02)         M1              L5 point")
for lab, tr, dd, m0, m2, md, m1, l5 in rows:
    print(f"     {lab:26s} {tr:.3f}  {fmt(dd)}  {fmt(m0):14s}  {fmt(m2)}  " + "  ".join(fmt(x) for x in md) + f"  {fmt(m1)}  {l5:.3f}")
ok("Theorem 2 / ladder: DD ⊇ M2 ⊇ M1(.10) ⊇ M1(.05) ⊇ M1(.02) ⊇ M1 for all targets", nest_ok)
ok("Theorem 2: truth inside DD, M2, every M1(delta) and M1", cover_ok)
ok("Prop D6: M1(delta) width non-increasing as delta falls", mono_ok)

# Prop D6: large delta reproduces M2
r, c = 0, 1
m2 = district_lp(Qd, Vd, mu_r, mu_d, r, c, "M2")
mbig = district_lp(Qd, Vd, mu_r, mu_d, r, c, "M1d", 10.0, zhat)
ok("Prop D6: M1(delta) = M2 for delta large", abs(m2[0] - mbig[0]) < 1e-7 and abs(m2[1] - mbig[1]) < 1e-7)

# Prop D5: M2 bound = sum of district-by-district bounds
s_lo = s_hi = 0.0
for d in range(D):
    b_ = district_lp([Qd[d]], [Vd[d]], mu_r, [mu_d[d]], r, c, "M2")
    s_lo += b_[0]; s_hi += b_[1]
ok("Prop D5: M2 bound equals the (weighted) sum of district bounds",
   abs(s_lo - m2[0]) < 1e-7 and abs(s_hi - m2[1]) < 1e-7)

# Prop D4(a): with zeta = 0, M1 intersection = pooled outer set on conv(union of A_d)
dist0 = make_districts(H, [[np.zeros(R) for _ in range(C)] for _ in range(D)], alphas)
Qd0 = [quad_fit(d["X"], d["T"]) for d in dist0]; Vd0 = [hull_vertices(d["X"]) for d in dist0]
mu_r0, mu_d0 = moments(dist0)
X0 = np.vstack([d["X"] for d in dist0]); T0 = np.vstack([d["T"] for d in dist0])
eq = True
for r, c, _ in targets:
    a = district_lp(Qd0, Vd0, mu_r0, mu_d0, r, c, "M1", zeta_c=zeta_from_Q(Qd0))
    b_ = district_lp([quad_fit(X0, T0)], [hull_vertices(X0)], mu_r0, [sum(mu_d0)], r, c, "M2")
    eq &= abs(a[0] - b_[0]) < 1e-6 and abs(a[1] - b_[1]) < 1e-6
ok("Prop D4(a): if zeta = 0, riding dummies give exactly the pooled bounds", eq)

# Prop D7 (L5): with zero contextual effects, L5 point = truth and lies in the M1 bound
Hflat = [np.outer(np.array([H[c][r_, r_] for r_ in range(R)]), np.ones(R)) for c in range(C)]
# Hflat_c[r, s] = H_c[r, r] for all s: rates constant in composition; rows sum to 1 by construction
assert np.allclose(sum(Hflat), np.ones((R, R)))
while True:
    zf = random_zeta(D, R, C, 0.04)
    try:
        distf = make_districts(Hflat, zf, alphas); break
    except AssertionError:
        continue
Qf = [quad_fit(d["X"], d["T"]) for d in distf]; Vf = [hull_vertices(d["X"]) for d in distf]
mu_rf, mu_df = moments(distf)
good = True
for r, c, _ in targets:
    tr = truth(distf, r, c); l5 = goodman_within(distf, r, c)
    m1 = district_lp(Qf, Vf, mu_rf, mu_df, r, c, "M1", zeta_c=zeta_from_Q(Qf))
    good &= abs(l5 - tr) < 1e-8 and m1[0] - TOL <= l5 <= m1[1] + TOL
ok("Prop D7: under within-district constancy, L5 point = truth and lies inside the M1 bound", good)

# Prop D3: R = 2 closed form under M1 vs LP (random 2xC designs)
def m1_r2_closed(Qd, Vd, zeta_c, c):
    lo, hi = -np.inf, np.inf
    for d in range(len(Qd)):
        Q = Qd[d][c]; Q0 = Qd[0][c]
        for v in Vd[d]:
            t = v[0]
            z1, z2 = zeta_c[d][c]
            # K^(d) = K + skew(zeta 1'): kappa^(d) = kappa + (z1 - z2)/2
            sh = (z1 - z2) / 2
            lo = max(lo, -(Q[0, 0] * t + Q[0, 1] * (1 - t)) / (1 - t) - sh)
            hi = min(hi, (Q[0, 1] * t + Q[1, 1] * (1 - t)) / t - sh)
    return lo, hi


def m1_r2_lp_kappa(Qd, Vd, zeta_c, c):
    """LP range of the common kappa_c (reference district) under M1, R = 2, all columns coupled."""
    C_ = len(Qd[0]); Kb, d0 = basis(2, C_)
    A, b = [], []
    for d in range(len(Qd)):
        shift = [skew(np.outer(zeta_c[d][cc], np.ones(2))) for cc in range(C_)]
        for v in Vd[d]:
            for cc in range(C_):
                for rr in range(2):
                    A.append([-(Kb[j][cc] @ v)[rr] for j in range(d0)])
                    # (Q^(d) + K + skew(zeta_d 1')) v >= 0  <=>  -(K v) <= (Q^(d) + skew(zeta_d 1')) v
                    b.append(((Qd[d][cc] + shift[cc]) @ v)[rr])
    obj = np.zeros(d0); obj[c] = 1
    lo = linprog(obj, A_ub=A, b_ub=b, bounds=[(None, None)] * d0, method="highs")
    hi = linprog(-obj, A_ub=A, b_ub=b, bounds=[(None, None)] * d0, method="highs")
    return lo.fun, -hi.fun


bad = n = 0
for trial in range(60):
    C2 = int(rng.integers(2, 5))
    G = rng.gamma(1.0, size=(C2, 2, 2)); H2 = [G[c] / G.sum(0) for c in range(C2)]
    D2 = int(rng.integers(2, 6))
    al = [rng.uniform(0.5, 5, 2) for _ in range(D2)]
    for _ in range(50):
        z = random_zeta(D2, 2, C2, 0.03)
        try:
            ds = make_districts(H2, z, al, polls=60); break
        except AssertionError:
            continue
    else:
        continue
    Q2 = [quad_fit(d["X"], d["T"]) for d in ds]; V2 = [hull_vertices(d["X"]) for d in ds]
    zh = zeta_from_Q(Q2)
    # closed form: per-column box (intersection over districts and their vertices), then Lemma 4(b) projection
    box = [m1_r2_closed(Q2, V2, zh, c) for c in range(C2)]
    lo_v = np.array([b_[0] for b_ in box]); hi_v = np.array([b_[1] for b_ in box])
    for c in range(C2 - 1):
        pl = max(lo_v[c], -(hi_v.sum() - hi_v[c])); pu = min(hi_v[c], -(lo_v.sum() - lo_v[c]))
        L, U = m1_r2_lp_kappa(Q2, V2, zh, c)
        n += 1
        bad += not (abs(L - pl) < 1e-7 and abs(U - pu) < 1e-7)
ok(f"Prop D3: closed-form M1 interval for kappa_c (R = 2) matches the LP ({n} column-instances)", bad == 0)

# Prop D0: model-free cell bounds are attained (feasible tables at the Frechet extremes)
bad = 0
for _ in range(300):
    Rr, Cc = int(rng.integers(2, 5)), int(rng.integers(2, 5))
    x = rng.dirichlet(np.ones(Rr)); t = rng.dirichlet(np.ones(Cc))
    r_, c_ = int(rng.integers(Rr)), int(rng.integers(Cc))
    for target in (max(0, x[r_] + t[c_] - 1), min(x[r_], t[c_])):
        # feasibility of {n >= 0, row sums x, col sums t, n[r_,c_] = target}
        Aeq = []; beq = []
        for i in range(Rr):
            row = np.zeros(Rr * Cc); row[i * Cc:(i + 1) * Cc] = 1; Aeq.append(row); beq.append(x[i])
        for j in range(Cc):
            row = np.zeros(Rr * Cc); row[j::Cc] = 1; Aeq.append(row); beq.append(t[j])
        row = np.zeros(Rr * Cc); row[r_ * Cc + c_] = 1; Aeq.append(row); beq.append(target)
        res = linprog(np.zeros(Rr * Cc), A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, None)] * (Rr * Cc), method="highs")
        bad += res.status != 0
ok("Prop D0: per-cell Frechet extremes are attained by feasible tables (300 random cases)", bad == 0)

# 2x2 worked example (Phase 1 section 7.3) through the ladder
def ex2x2():
    g = 0.10
    spec = [((0.5, 0.9), 0.89, 0.20), ((0.1, 0.5), 0.45, 0.02)]
    ds = []
    for (a, b_), aP, aO in spec:
        t_ = np.linspace(a, b_, 2001); X = np.c_[t_, 1 - t_]
        loy, cro = aP + g * t_, aO + g * t_
        M = np.stack([np.c_[loy, 1 - loy], np.c_[cro, 1 - cro]], 1)
        T = np.einsum("ir,irc->ic", X, M)
        ds.append(dict(X=X, T=T, N=np.ones(len(t_)), M=M))
    return ds


ds = ex2x2()
Q2 = [quad_fit(d["X"], d["T"]) for d in ds]; V2 = [hull_vertices(d["X"]) for d in ds]
mu2, mud2 = moments(ds)
tr = truth(ds, 0, 0)
m2 = district_lp(Q2, V2, mu2, mud2, 0, 0, "M2")
m1 = district_lp(Q2, V2, mu2, mud2, 0, 0, "M1", zeta_c=zeta_from_Q(Q2))
l5 = goodman_within(ds, 0, 0)
X2 = np.vstack([d["X"] for d in ds]); T2 = np.vstack([d["T"] for d in ds])
m0 = district_lp([quad_fit(X2, T2)], [hull_vertices(X2)], mu2, [sum(mud2)], 0, 0, "M2")
print(f"     2x2 example (P1 section 7.3): truth {tr:.3f}  DD {fmt(dd_bound(ds, 0, 0))}  M0 {fmt(m0)}  "
      f"M2 {fmt(m2)}  M1 {fmt(m1)}  L5 point {l5:.3f}")
ok("2x2 example: M1 bound reproduces Phase 1 [0.618, 0.866]", abs(m1[0] - 0.618) < 5e-4 and abs(m1[1] - 0.866) < 5e-4)
ok("2x2 example: pooled M0 outer set is empty (as in Phase 1)", m0 is None)
