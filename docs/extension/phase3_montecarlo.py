# Monte Carlo for Phase 3 (PHASE3_inference.md): coverage of the finite-population
# rate B by (i) the plug-in bound and (ii) the inflated-constraint confidence
# interval of Theorem 5.
#   Design A: 2x2, 93 districts x 10 units, model M1 (district baselines, shared slopes).
#   Design B: 3x3, pooled, 1000 units.
# Unit rates are Beta/Dirichlet draws whose conditional means are exactly linear
# in composition, so the model holds. Run: python3 -I phase3_montecarlo.py [reps]
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
from scipy.stats import norm

REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 300
ALPHA1, ALPHA2 = 0.04, 0.01        # conservative method: estimation part, finite-population part
ALPHA_SIDE_EST, ALPHA_SIDE_E = 0.02, 0.005   # CLR method, per side (two sides: total 0.05)
ALPHA_SIDE_E_BD, ALPHA_V = 0.0025, 0.005     # Bhatia-Davis variant: E-term per side, variance-bound CI


def bd_sd(w, mlo, mhi):
    """Bhatia-Davis sd bound for E = sum w_i e_i given a CI [mlo, mhi] for the w^2-weighted mean rate."""
    lo, hi = min(max(mlo, 0.0), 1.0), min(max(mhi, 0.0), 1.0)
    vmax = 0.25 if lo <= 0.5 <= hi else max(lo * (1 - lo), hi * (1 - hi))
    return np.sqrt(np.sum(w ** 2) * vmax)
NDRAW = 4000


def hc3_cov(Z, E):
    """HC3 covariance of stacked OLS coefficients for several responses sharing design Z.
    E: n x C residuals. Returns (C*k) x (C*k) matrix, ordering column-major by response."""
    n, k = Z.shape
    ZtZi = np.linalg.inv(Z.T @ Z)
    h = np.einsum("ij,jk,ik->i", Z, ZtZi, Z)
    Eh = E / (1 - h)[:, None]
    C = E.shape[1]
    V = np.zeros((C * k, C * k))
    for a in range(C):
        for b in range(a, C):
            S = (Z * (Eh[:, a] * Eh[:, b])[:, None]).T @ Z
            blk = ZtZi @ S @ ZtZi
            V[a * k:(a + 1) * k, b * k:(b + 1) * k] = blk
            V[b * k:(b + 1) * k, a * k:(a + 1) * k] = blk.T
    return V


def crit(G, V, rng, two_sided_rows=()):
    """(1 - ALPHA1) quantile of max_j |or signed| (G_j Z) / se_j, Z ~ N(0, V)."""
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", G, V, G), 1e-300))
    L = np.linalg.cholesky(V + 1e-14 * np.eye(len(V)))
    Zs = rng.standard_normal((NDRAW, len(V))) @ L.T
    W = (Zs @ G.T) / se
    W[:, list(two_sided_rows)] = np.abs(W[:, list(two_sided_rows)])
    return np.quantile(W.max(1), 1 - ALPHA1), se


def clr_upper(theta_hat, Gth, V, rng, alpha, n):
    """CLR-type intersection bound for theta_U = min_j theta_j, theta_j linear in beta with gradients Gth.
    Returns an upper confidence bound theta_hat_U(1 - alpha) (Chernozhukov, Lee & Rosen 2013, with
    adaptive inequality selection; analytic Gaussian approximation)."""
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Gth, V, Gth), 1e-300))
    L = np.linalg.cholesky(V + 1e-14 * np.eye(len(V)))
    Zs = (rng.standard_normal((NDRAW, len(V))) @ L.T) @ Gth.T / se          # standardised errors
    gam_n = 1 - 0.1 / np.log(n)
    kp = np.quantile(Zs.max(1), gam_n)
    up = np.min(theta_hat + kp * se)
    sel = theta_hat <= up + 2 * kp * se                                      # near-binding set V_hat
    k = np.quantile(Zs[:, sel].max(1), 1 - alpha)
    return np.min(theta_hat[sel] + k * se[sel]), sel.sum()


# ======================= Design A: 2x2, M1 =======================
D, PD = 93, 10


def sim_A(rng):
    m = rng.uniform(0.1, 0.9, D)
    zw = 0.4 * (m - 0.5) + 0.05 * rng.normal(size=D)        # district shifts, correlated with mean x
    zb = 0.4 * (m - 0.5) + 0.05 * rng.normal(size=D)
    d = np.repeat(np.arange(D), PD)
    x = np.clip(m[d] + rng.uniform(-0.15, 0.15, D * PD), 0.01, 0.99)
    mw = 0.55 + zw[d] - 0.10 * x                               # E[beta_w | x, d]
    mb = 0.35 + zb[d] + 0.20 * x                               # E[beta_b | x, d]
    conc = 60.0
    assert 0 < mw.min() and mw.max() < 1 and 0 < mb.min() and mb.max() < 1
    bw = rng.beta(conc * mw, conc * (1 - mw)); bb = rng.beta(conc * mb, conc * (1 - mb))
    N = np.full(x.size, 400.0)
    return d, x, x * bb + (1 - x) * bw, bb, N, mb


def design_A(d, x):
    """M1 regression of T on homogeneous quadratic in X=(x,1-x) + district shifts zeta_d' X.
    Columns: Q11 (x^2), 2Q12 (x(1-x)), Q22 ((1-x)^2), then for d>=1: zeta1_d * x, zeta2_d * (1-x)."""
    n = x.size
    Z = np.zeros((n, 3 + 2 * (D - 1)))
    Z[:, 0], Z[:, 1], Z[:, 2] = x * x, x * (1 - x), (1 - x) ** 2
    for k in range(1, D):
        sel = d == k
        Z[sel, 3 + 2 * (k - 1)] = x[sel]
        Z[sel, 4 + 2 * (k - 1)] = 1 - x[sel]
    return Z


def lin_A(beta_len, d, x, N, w_override=None):
    """Linear maps from the coefficient vector to: constraint RHS terms s_j, kappa-coefficients a_j,
    and the constant of B~ (finite-population target with sample weights)."""
    k = beta_len
    def Qd_coef(dd):
        # returns coefficient vectors (len k) for Q11^(d), Q12^(d), Q22^(d) and shift sh_d = (z1 - z2)/2
        e = np.eye(k)
        q11, q12, q22 = e[0].copy(), e[1] / 2, e[2].copy()
        sh = np.zeros(k)
        if dd >= 1:
            z1, z2 = e[3 + 2 * (dd - 1)], e[4 + 2 * (dd - 1)]
            q11 = q11 + z1; q22 = q22 + z2; q12 = q12 + (z1 + z2) / 2
            sh = (z1 - z2) / 2
        return q11, q12, q22, sh
    S, A = [], []
    for dd in range(D):
        q11, q12, q22, sh = Qd_coef(dd)
        xs = x[d == dd]
        for t in (xs.min(), xs.max()):
            # column 1 (vote), K^(d)_12 = kappa + sh, H^(d) = Q^(d) + K^(d)
            # row 1: Q11 t + (Q12 + kappa + sh)(1-t) >= 0 ; row 2: (Q12 - kappa - sh) t + Q22 (1-t) >= 0
            S.append(q11 * t + (q12 + sh) * (1 - t)); A.append(1 - t)
            S.append((q12 - sh) * t + q22 * (1 - t)); A.append(-t)
            # column 2 = J - column 1 (rates 1 - ...): row r: 1 - [row r of column 1] >= 0
            S.append(-(q11 * t + (q12 + sh) * (1 - t))); A.append(-(1 - t))   # plus constant 1 below
            S.append(-((q12 - sh) * t + q22 * (1 - t))); A.append(t)
    S = np.array(S); A = np.array(A)
    one = np.zeros(len(S)); one[2::4] = 1; one[3::4] = 1
    # B~ = sum_i w_i [Q11^(d) x_i + (Q12^(d) + sh_d + kappa)(1 - x_i)],  w_i = N_i x_i / sum N x
    w = N * x / np.sum(N * x) if w_override is None else w_override
    c = np.zeros(k); gam = np.sum(w * (1 - x))
    for dd in range(D):
        q11, q12, q22, sh = Qd_coef(dd)
        sel = d == dd
        c += np.sum(w[sel] * x[sel]) * q11 + np.sum(w[sel] * (1 - x[sel])) * (q12 + sh)
    return S, A, one, c, gam, w


def kappa_range(sval, A, one):
    """{kappa: sval_j + one_j + A_j kappa >= 0 for all j}"""
    rhs = sval + one
    lo = max((-rhs[j] / A[j] for j in range(len(A)) if A[j] > 0), default=-np.inf)
    hi = min((-rhs[j] / A[j] for j in range(len(A)) if A[j] < 0), default=np.inf)
    return lo, hi


def sim_C(rng, n=1000):
    """Pooled 2x2 with the truth on the boundary of the identified set: group-b rate ~0 at the
    most group-b-heavy units (cf. JKST Example 3)."""
    x = rng.uniform(0.02, 0.98, n)
    mb = 0.03 * (1 - x); mw = 0.95 - 0.9 * x
    conc = 60.0
    bb = rng.beta(conc * mb, conc * (1 - mb)); bw = rng.beta(conc * mw, conc * (1 - mw))
    return np.zeros(n, int), x, x * bb + (1 - x) * bw, bb, np.full(n, 400.0), mb


def run_A(seed, sim=None):
    rng = np.random.default_rng(seed)
    d, x, t, bb, N, mb = (sim or sim_A)(rng)
    truth = np.sum(N * x * bb) / np.sum(N * x)
    truth_cm = np.sum(N * x * mb) / np.sum(N * x)          # conditional-mean rate B~
    Z = design_A(d, x)
    beta = np.linalg.lstsq(Z, t, rcond=None)[0]
    V = hc3_cov(Z, (t - Z @ beta)[:, None])
    S, A, one, c, gam, w = lin_A(len(beta), d, x, N)
    s_hat = S @ beta
    # plug-in
    lo, hi = kappa_range(s_hat, A, one)
    plug = None if lo > hi else (c @ beta + gam * lo, c @ beta + gam * hi)
    # inflated (Theorem 5): simultaneous over all constraint terms (one-sided) and the constant (two-sided)
    G = np.vstack([S, c])
    k, se = crit(G, V, rng, two_sided_rows=[len(S)])
    lo, hi = kappa_range(s_hat + k * se[:-1], A, one)
    S1 = 0.5 * np.sqrt(np.sum(w ** 2))
    z2 = norm.ppf(1 - ALPHA2 / 2)
    ci = (c @ beta - k * se[-1] + gam * lo - z2 * S1, c @ beta + k * se[-1] + gam * hi + z2 * S1)
    # CLR-type: candidate endpoints theta_j(beta) = c.beta + gam * (-(S_j.beta + one_j)/A_j)
    up_j = A < 0; lo_j = A > 0
    GU = c[None, :] - gam * S[up_j] / A[up_j, None]; iU = -gam * one[up_j] / A[up_j]
    GL = c[None, :] - gam * S[lo_j] / A[lo_j, None]; iL = -gam * one[lo_j] / A[lo_j]
    a_side, a_e = ALPHA_SIDE_EST, ALPHA_SIDE_E
    zE = norm.ppf(1 - a_e)
    U, _ = clr_upper(GU @ beta + iU, GU, V, rng, a_side, len(x))
    Lneg, _ = clr_upper(-(GL @ beta + iL), -GL, V, rng, a_side, len(x))
    clr = (-Lneg - zE * S1, U + zE * S1)
    # target B~ (conditional-mean rate): no finite-population term; per side alpha = 0.025
    U2, _ = clr_upper(GU @ beta + iU, GU, V, rng, 0.025, len(x))
    L2, _ = clr_upper(-(GL @ beta + iL), -GL, V, rng, 0.025, len(x))
    # Bhatia-Davis allowance: Var(E) <= sum w^2 * mbar(1 - mbar), mbar = pi-weighted mean rate
    # (pi = w^2 / sum w^2), bounded by its own CLR CI at total level ALPHA_V
    wpi = w ** 2 / np.sum(w ** 2)
    _, _, _, cp, gp, _ = lin_A(len(beta), d, x, N, w_override=wpi)
    GUp = cp[None, :] - gp * S[up_j] / A[up_j, None]; iUp = -gp * one[up_j] / A[up_j]
    GLp = cp[None, :] - gp * S[lo_j] / A[lo_j, None]; iLp = -gp * one[lo_j] / A[lo_j]
    mU, _ = clr_upper(GUp @ beta + iUp, GUp, V, rng, ALPHA_V / 2, len(x))
    mLn, _ = clr_upper(-(GLp @ beta + iLp), -GLp, V, rng, ALPHA_V / 2, len(x))
    S1bd = bd_sd(w, -mLn, mU)
    zE2 = norm.ppf(1 - ALPHA_SIDE_E_BD)
    clr_bd = (-Lneg - zE2 * S1bd, U + zE2 * S1bd)
    return truth, plug, ci, clr, truth_cm, (-L2, U2), clr_bd


# ======================= Design B: 3x3 pooled =======================
def M3(x):
    xP, xO, xN = x.T
    return np.stack([np.stack([0.08 - 0.03 * xP, 0.80 + 0.10 * xP, 0.12 - 0.07 * xP], -1),
                     np.stack([0.10 + 0.02 * xO, 0.10 - 0.05 * xO, 0.80 + 0.03 * xO], -1),
                     np.stack([0.75 - 0.10 * xN, 0.15 + 0.05 * xP, 0.10 + 0.10 * xN - 0.05 * xP], -1)], 1)


PAIRS = [(0, 1), (0, 2), (1, 2)]
MONO = [(r, s) for r in range(3) for s in range(r, 3)]
TARGETS = [(0, 1), (0, 2), (0, 0), (2, 1)]


def run_B(seed, p=1000, conc=40.0):
    rng = np.random.default_rng(seed)
    X = rng.dirichlet([4, 3, 2.5], p)
    M = M3(X)
    beta_u = np.stack([np.stack([rng.dirichlet(conc * M[i, r]) for r in range(3)]) for i in range(p)])
    T = np.einsum("ir,irc->ic", X, beta_u)
    N = rng.integers(300, 1500, p).astype(float)
    Z = np.column_stack([X[:, r] * X[:, s] for r, s in MONO])
    coef = np.linalg.lstsq(Z, T, rcond=None)[0]                 # 6 x 3
    V = hc3_cov(Z, T - Z @ coef)                               # stacked by column c: block c = coef[:, c]
    kz = len(MONO)
    def Qcoef(c, r, s):
        """coefficient vector (len 3*kz) giving Q_c[r, s] as a linear function of stacked coefs"""
        g = np.zeros(3 * kz)
        j = MONO.index((min(r, s), max(r, s)))
        g[c * kz + j] = 1.0 if r == s else 0.5
        return g
    # variables kk: K_c[pair] for c in {0,1} (K_2 = -K_0 - K_1): 6 variables
    def Kcoef(c, r, s):
        a = np.zeros(6)
        if r == s: return a
        sign = 1.0 if r < s else -1.0
        j = PAIRS.index((min(r, s), max(r, s)))
        if c < 2: a[c * 3 + j] = sign
        else: a[j] = -sign; a[3 + j] = -sign
        return a
    V3 = list(X[ConvexHull(X[:, :2]).vertices])
    S, A = [], []
    for v in V3:
        for c in range(3):
            for r in range(3):
                S.append(sum(Qcoef(c, r, s) * v[s] for s in range(3)))
                A.append(sum(Kcoef(c, r, s) * v[s] for s in range(3)))
    S = np.array(S); A = np.array(A)
    beta = coef.T.reshape(-1)                                  # stacked by column c
    res = []
    for (r, c) in TARGETS:
        w = N * X[:, r] / np.sum(N * X[:, r])
        cvec = sum(np.sum(w * X[:, r] * X[:, s]) / np.sum(w * X[:, r]) * 0 + np.sum(w * X[:, s]) * Qcoef(c, r, s)
                   for s in range(3))
        gam = sum(np.sum(w * X[:, s]) * Kcoef(c, r, s) for s in range(3))
        truth = np.sum(w * beta_u[:, r, c])
        truth_cm = np.sum(w * M[:, r, c])
        def lp(svals):
            # constraints: svals + A kk >= 0  ->  -A kk <= svals
            kw = dict(A_ub=-A, b_ub=svals, bounds=[(None, None)] * 6, method="highs")
            lo, hi = linprog(gam, **kw), linprog(-gam, **kw)
            return None if lo.status or hi.status else (lo.fun, -hi.fun)
        b0 = lp(S @ beta)
        plug = None if b0 is None else (cvec @ beta + b0[0], cvec @ beta + b0[1])
        G = np.vstack([S, cvec])
        k, se = crit(G, V, rng, two_sided_rows=[len(S)])
        b1 = lp(S @ beta + k * se[:-1])
        S1 = 0.5 * np.sqrt(np.sum(w ** 2))
        z2 = norm.ppf(1 - ALPHA2 / 2)
        ci = None if b1 is None else (cvec @ beta - k * se[-1] + b1[0] - z2 * S1,
                                      cvec @ beta + k * se[-1] + b1[1] + z2 * S1)
        # CLR-type via dual vertices collected from draws of the coefficients
        Lc = np.linalg.cholesky(V + 1e-14 * np.eye(len(V)))
        def duals(objsign):
            found = {}
            for b_ in [beta] + [beta + Lc @ rng.standard_normal(len(V)) for _ in range(150)]:
                res_ = linprog(-objsign * gam, A_ub=-A, b_ub=S @ b_, bounds=[(None, None)] * 6, method="highs")
                if res_.status: continue
                lam = -res_.ineqlin.marginals          # >= 0; objective value = -(lam . s)
                key = tuple(np.round(lam, 9))
                found[key] = lam
            return np.array(list(found.values()))
        LU = duals(+1.0)   # max gam.k = min_lam lam.s  (lam >= 0, A'lam = gam)
        LL = duals(-1.0)   # max -gam.k = min_lam lam.s  ->  min gam.k = -min_lam lam.s
        GU = cvec[None, :] + LU @ S
        GL = cvec[None, :] - LL @ S
        zE = norm.ppf(1 - ALPHA_SIDE_E)
        U, nU = clr_upper(GU @ beta, GU, V, rng, ALPHA_SIDE_EST, p)
        Ln, nL = clr_upper(-(GL @ beta), -GL, V, rng, ALPHA_SIDE_EST, p)
        clr = (-Ln - zE * S1, U + zE * S1)
        U2, _ = clr_upper(GU @ beta, GU, V, rng, 0.025, p)
        L2, _ = clr_upper(-(GL @ beta), -GL, V, rng, 0.025, p)
        # Bhatia-Davis variant
        wpi = w ** 2 / np.sum(w ** 2)
        cp = sum(np.sum(wpi * X[:, s]) * Qcoef(c, r, s) for s in range(3))
        gp = sum(np.sum(wpi * X[:, s]) * Kcoef(c, r, s) for s in range(3))
        gam_save = gam; gam = gp
        LUp, LLp = duals(+1.0), duals(-1.0)
        gam = gam_save
        GUp = cp[None, :] + LUp @ S; GLp = cp[None, :] - LLp @ S
        mU, _ = clr_upper(GUp @ beta, GUp, V, rng, ALPHA_V / 2, p)
        mLn, _ = clr_upper(-(GLp @ beta), -GLp, V, rng, ALPHA_V / 2, p)
        S1bd = bd_sd(w, -mLn, mU)
        zE2 = norm.ppf(1 - ALPHA_SIDE_E_BD)
        clr_bd = (-Ln - zE2 * S1bd, U + zE2 * S1bd)
        res.append((truth, plug, ci, clr, truth_cm, (-L2, U2), clr_bd))
    return res


def summarise(name, rows):
    cov = lambda iv, tr: iv is not None and iv[0] <= tr <= iv[1]
    wid = lambda ivs: np.median([iv[1] - iv[0] for iv in ivs if iv is not None]) if any(iv is not None for iv in ivs) else float("nan")
    out = [name]
    for k_, lab in ((1, "plug-in"), (2, "conservative"), (3, "CLR-type")):
        out.append(f"{lab}: cover {np.mean([cov(r[k_], r[0]) for r in rows]):.3f}, "
                   f"empty {np.mean([r[k_] is None for r in rows]):.2f}, width {wid([r[k_] for r in rows]):.3f}")
    out.append(f"CLR for B~: cover {np.mean([cov(r[5], r[4]) for r in rows]):.3f}, width {wid([r[5] for r in rows]):.3f}")
    out.append(f"CLR for B (Bhatia-Davis): cover {np.mean([cov(r[6], r[0]) for r in rows]):.3f}, width {wid([r[6] for r in rows]):.3f}")
    print("  " + f"{out[0]:26s} " + " | ".join(out[1:]))


if __name__ == "__main__":
    for D, PD in [(10, 93), (93, 30), (93, 10)]:
        print(f"Design A: 2x2, M1, {D} districts x {PD} units, {REPS} replications")
        summarise(f"B (D={D}, PD={PD})", [run_A(1000 + s) for s in range(REPS)])
    print(f"Design C: 2x2 pooled, 1000 units, truth on the boundary, {REPS} replications")
    D = 1
    summarise("B (boundary case)", [run_A(3000 + s, sim_C) for s in range(REPS)])
    print(f"Design B: 3x3 pooled, 1000 units, {REPS} replications")
    out = [run_B(5000 + s) for s in range(REPS)]
    for j, lab in enumerate(["loyalty P->P", "defection P->O", "demobilisation P->abstain", "mobilisation N->P"]):
        summarise(lab, [o[j] for o in out])
