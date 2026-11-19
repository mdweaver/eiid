# Monte Carlo for Phase 3 (PHASE3_inference.md, section 8), rerun for review item R-37 (WS6).
#
# Coverage of the realised rate B and of the conditional-mean rate B~ by
#   * the plug-in estimate of the identified set (not a CI),
#   * the conservative CI (Theorem P3.5(a)),
#   * the CLR CI as implemented (Conjecture P3.5(b3); Popoviciu and Bhatia-Davis allowances for B),
#   * the projection CI (Proposition P3.5(b2)(ii): k = sqrt(chi2_{q,1-alpha}), q = rank of Sigma-hat,
#     applied to the same candidates as CLR). q = dim(beta) in every 2x2 design (full rank); in B_3x3
#     the rank is 12 of 18 because of the exact adding-up reduction (case (ii-b)). For B the
#     Bhatia-Davis allowance takes mbar from the PROJECTION interval at level ALPHA_V (proj_B, covered
#     by P3.5(b2)(ii) + P3.5(c1)); the earlier CLR-based version is kept as proj_B_clrbd ("not covered").
#     (Changed after the round-2 review, Nov 2026. Neither uses random draws, so every other column
#     reproduces the earlier run exactly.)
#   * the finite-sample CI (Theorem P3.8; 2x2 M1 designs only), clipped to [0, 1].
# Every CI is nominal 95%. Widths are medians after clipping to [0, 1] (the logical range).
#
# Designs (unit rates are Beta/Dirichlet draws whose conditional means are exactly linear in
# composition given the district, so M1 holds; M0 holds only when the district shifts are zero):
#   A    2x2, M1, district shifts correlated with district mean composition; D districts x PD polls:
#        10x93, 93x30, and the R2 series 30x10, 93x10, 200x10 (with Theorem P3.8).
#   C    2x2 pooled, truth near the boundary (group-b rate ~ 0 at the most group-b-heavy units).
#   E    2x2 pooled, Slater fails: true rates exactly 0 at both range vertices, so the identified
#        set for the model-implied rate is a single point.
#   S    riding shocks, 2x2, 20 districts x 50 polls, fitted by M0 and by M1, with the pooling
#        test (Wald test of all district shifts = 0, P3.6(c)):
#          S0 no shifts (size of the pooling test), S1 shifts affine in district composition plus
#          noise, S2 a common shock to both groups, non-affine in district composition
#          (0.25 cos(2 pi m_d)).
#   B    3x3 pooled, 1000 units, four targets.
#
# Run: python3 -I phase3_montecarlo.py [reps] [design ...]     (default: 1000 reps, all designs)
# Writes a summary to phase3_montecarlo_summary.json and per-replication intervals to $MC_RAW_DIR
# (default /tmp/phase3_mc_raw). The 1000-replication run took about an hour on 4 cores.
# Uses all cores (one process per core, one BLAS thread each). Seeds are fixed per design and
# replication, so results do not depend on the number of processes.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import json
import gzip
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
from scipy.stats import norm, chi2

ALPHA = 0.05
ALPHA1, ALPHA2 = 0.04, 0.01        # conservative CI for B: estimation part, finite-population part
ALPHA_SIDE_EST, ALPHA_SIDE_E = 0.02, 0.005   # CLR for B, per side (two sides: total 0.05)
ALPHA_SIDE_E_BD, ALPHA_V = 0.0025, 0.005     # Bhatia-Davis variant: E-term per side, variance-bound CI
NDRAW = 4000


def bd_sd(w, mlo, mhi):
    """Bhatia-Davis sd bound for E = sum w_i e_i given a CI [mlo, mhi] for the w^2-weighted mean rate."""
    lo, hi = min(max(mlo, 0.0), 1.0), min(max(mhi, 0.0), 1.0)
    vmax = 0.25 if lo <= 0.5 <= hi else max(lo * (1 - lo), hi * (1 - hi))
    return np.sqrt(np.sum(w ** 2) * vmax)


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


def sqrt_psd(V):
    ev, U = np.linalg.eigh((V + V.T) / 2)
    return U * np.sqrt(np.clip(ev, 0, None))


def rank_psd(V):
    ev = np.linalg.eigvalsh((V + V.T) / 2)
    return max(1, int(np.sum(ev > 1e-10 * ev.max())))


def crit(G, V, rng, alpha, two_sided_rows=()):
    """(1 - alpha) quantile of max_j (G_j Z) / se_j (|.| for two-sided rows), Z ~ N(0, V)."""
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", G, V, G), 1e-300))
    L = sqrt_psd(V)
    Zs = rng.standard_normal((NDRAW, len(V))) @ L.T
    W = (Zs @ G.T) / se
    W[:, list(two_sided_rows)] = np.abs(W[:, list(two_sided_rows)])
    return np.quantile(W.max(1), 1 - alpha), se


def clr_upper(theta_hat, Gth, V, rng, alpha, n):
    """CLR-type intersection bound for theta_U = min_j theta_j, theta_j linear in beta with gradients Gth.
    Returns an upper confidence bound (Chernozhukov, Lee & Rosen 2013, adaptive inequality selection,
    analytic Gaussian approximation), as in the package's .rxc_clr_min."""
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Gth, V, Gth), 1e-300))
    if len(theta_hat) == 1:
        return theta_hat[0] + norm.ppf(1 - alpha) * se[0]
    L = sqrt_psd(V)
    Zs = (rng.standard_normal((NDRAW, len(V))) @ L.T) @ Gth.T / se
    gam_n = 1 - 0.1 / np.log(n)
    kp = np.quantile(Zs.max(1), gam_n)
    up = np.min(theta_hat + kp * se)
    sel = theta_hat <= up + 2 * kp * se
    k = np.quantile(Zs[:, sel].max(1), 1 - alpha)
    return np.min(theta_hat[sel] + k * se[sel])


def proj_upper(theta_hat, Gth, V, k):
    """Projection bound (P3.5(b2)(ii)): min_j (theta_j + k se_j)."""
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Gth, V, Gth), 0))
    return np.min(theta_hat + k * se)


# ======================= 2x2 designs (A, C, E, S) =======================
def sim_A(rng, D, PD, shift="affine"):
    """District shifts correlated with district mean composition m_d.
    shift: 'affine' (affine in m_d plus noise; design A, S1), 'none' (S0), 'nonaffine' (S2)."""
    m = rng.uniform(0.1, 0.9, D)
    if shift == "affine":
        zw = 0.4 * (m - 0.5) + 0.05 * rng.normal(size=D)
        zb = 0.4 * (m - 0.5) + 0.05 * rng.normal(size=D)
    elif shift == "none":
        zw = zb = np.zeros(D)
    elif shift == "nonaffine":
        zw = zb = 0.25 * np.cos(2 * np.pi * m)        # common shock to both groups, non-affine in m_d
    else:
        raise ValueError(shift)
    d = np.repeat(np.arange(D), PD)
    x = np.clip(m[d] + rng.uniform(-0.15, 0.15, D * PD), 0.01, 0.99)
    mw = 0.55 + zw[d] - 0.10 * x               # E[beta_w | x, d]
    mb = 0.35 + zb[d] + 0.20 * x               # E[beta_b | x, d]
    conc = 60.0
    assert 0 < mw.min() and mw.max() < 1 and 0 < mb.min() and mb.max() < 1
    bw = rng.beta(conc * mw, conc * (1 - mw)); bb = rng.beta(conc * mb, conc * (1 - mb))
    N = np.full(x.size, 400.0)
    return dict(d=d, x=x, t=x * bb + (1 - x) * bw, bb=bb, N=N, mb=mb, D=D)


def sim_C(rng, n=1000):
    """Pooled 2x2, truth near the boundary: group-b rate ~0 at the most group-b-heavy units (cf. JKST Ex. 3)."""
    x = rng.uniform(0.02, 0.98, n)
    mb = 0.03 * (1 - x); mw = 0.95 - 0.9 * x
    conc = 60.0
    bb = rng.beta(conc * mb, conc * (1 - mb)); bw = rng.beta(conc * mw, conc * (1 - mw))
    return dict(d=np.zeros(n, int), x=x, t=x * bb + (1 - x) * bw, bb=bb, N=np.full(n, 400.0), mb=mb, D=1)


def sim_E(rng, n=1000):
    """Pooled 2x2, Slater fails. Compositions range over [0.02, 0.98] with both ends observed.
    E[beta_b | x] = 0.30 (0.98 - x) / 0.96 is 0 at x = 0.98 and E[beta_w | x] = 0.50 (x - 0.02) / 0.96 is
    0 at x = 0.02, so the lower bound on kappa (from beta_b >= 0 at the right vertex) meets the upper
    bound (from beta_w >= 0 at the left vertex): the identified set is the single true value and no
    kappa satisfies every constraint strictly."""
    x = rng.uniform(0.02, 0.98, n); x[0], x[1] = 0.02, 0.98
    mb = 0.30 * (0.98 - x) / 0.96; mw = 0.50 * (x - 0.02) / 0.96
    conc = 60.0
    def draw(m):
        out = np.zeros_like(m); pos = m > 1e-12
        out[pos] = rng.beta(conc * m[pos], conc * (1 - m[pos]))
        return out
    bb = draw(mb); bw = draw(mw)
    return dict(d=np.zeros(n, int), x=x, t=x * bb + (1 - x) * bw, bb=bb, N=np.full(n, 400.0), mb=mb, D=1)


def design22(d, x, D):
    """M1 regression of T on the homogeneous quadratic in X=(x,1-x) plus district shifts zeta_d' X.
    Columns: Q11 (x^2), 2Q12 (x(1-x)), Q22 ((1-x)^2), then for d>=1: zeta1_d * x, zeta2_d * (1-x).
    D = 1 gives M0."""
    n = x.size
    Z = np.zeros((n, 3 + 2 * (D - 1)))
    Z[:, 0], Z[:, 1], Z[:, 2] = x * x, x * (1 - x), (1 - x) ** 2
    for k in range(1, D):
        sel = d == k
        Z[sel, 3 + 2 * (k - 1)] = x[sel]
        Z[sel, 4 + 2 * (k - 1)] = 1 - x[sel]
    return Z


def lin22(k, d, x, N, D, w_override=None):
    """Linear maps from the coefficient vector to the constraint terms s_j = S_j beta + one_j, the
    kappa-coefficients a_j, and the target B~ = c'beta + gam kappa (sample weights)."""
    def Qd_coef(dd):
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
            S.append(q11 * t + (q12 + sh) * (1 - t)); A.append(1 - t)          # beta_b(t) >= 0
            S.append((q12 - sh) * t + q22 * (1 - t)); A.append(-t)             # beta_w(t) >= 0
            S.append(-(q11 * t + (q12 + sh) * (1 - t))); A.append(-(1 - t))    # beta_b(t) <= 1
            S.append(-((q12 - sh) * t + q22 * (1 - t))); A.append(t)           # beta_w(t) <= 1
    S = np.array(S); A = np.array(A)
    one = np.zeros(len(S)); one[2::4] = 1; one[3::4] = 1
    w = N * x / np.sum(N * x) if w_override is None else w_override
    c = np.zeros(k); gam = np.sum(w * (1 - x))
    for dd in range(D):
        q11, q12, q22, sh = Qd_coef(dd)
        sel = d == dd
        c += np.sum(w[sel] * x[sel]) * q11 + np.sum(w[sel] * (1 - x[sel])) * (q12 + sh)
    return S, A, one, c, gam, w


def kappa_range(sval, A, one):
    rhs = sval + one
    lo = max((-rhs[j] / A[j] for j in range(len(A)) if A[j] > 0), default=-np.inf)
    hi = min((-rhs[j] / A[j] for j in range(len(A)) if A[j] < 0), default=np.inf)
    return lo, hi


def dd_width22(x, t, w):
    lo = np.maximum(0, (t - (1 - x)) / x); hi = np.minimum(1, t / x)
    return np.sum(w * (hi - lo))


def run22(seed, simf, simargs, models=("M1",), finite=False):
    rng = np.random.default_rng(seed)
    s = simf(rng, *simargs)
    d, x, t, bb, N, mb, D = s["d"], s["x"], s["t"], s["bb"], s["N"], s["mb"], s["D"]
    truth = np.sum(N * x * bb) / np.sum(N * x)
    truth_cm = np.sum(N * x * mb) / np.sum(N * x)
    w0 = N * x / np.sum(N * x)
    out = dict(B=truth, Bt=truth_cm, dd=dd_width22(x, t, w0))
    for model in models:
        dfit, Dfit = (d, D) if model == "M1" else (np.zeros_like(d), 1)
        Z = design22(dfit, x, Dfit)
        beta = np.linalg.lstsq(Z, t, rcond=None)[0]
        V = hc3_cov(Z, (t - Z @ beta)[:, None])
        S, A, one, c, gam, w = lin22(len(beta), dfit, x, N, Dfit)
        s_hat = S @ beta
        r = {}
        lo, hi = kappa_range(s_hat, A, one)
        r["plug"] = None if lo > hi else (c @ beta + gam * lo, c @ beta + gam * hi)
        S1 = 0.5 * np.sqrt(np.sum(w ** 2))
        # conservative (Theorem P3.5(a)), for B (alpha 0.04 + 0.01) and B~ (alpha 0.05)
        G = np.vstack([S, c])
        for key, a_est, aE in (("cons_B", ALPHA1, ALPHA2), ("cons_Bt", ALPHA, 0.0)):
            k, se = crit(G, V, rng, a_est, two_sided_rows=[len(S)])
            lo, hi = kappa_range(s_hat + k * se[:-1], A, one)
            zz = norm.ppf(1 - aE / 2) * S1 if aE > 0 else 0.0
            r[key] = None if lo > hi else (c @ beta - k * se[-1] + gam * lo - zz, c @ beta + k * se[-1] + gam * hi + zz)
        # candidates: in 2x2 every dual vertex is a single constraint, so all are used
        up_j = A < 0; lo_j = A > 0
        GU = c[None, :] - gam * S[up_j] / A[up_j, None]; iU = -gam * one[up_j] / A[up_j]
        GL = c[None, :] - gam * S[lo_j] / A[lo_j, None]; iL = -gam * one[lo_j] / A[lo_j]
        thU, thL = GU @ beta + iU, GL @ beta + iL
        n = len(x)
        U = clr_upper(thU, GU, V, rng, ALPHA_SIDE_EST, n)
        Ln = clr_upper(-thL, -GL, V, rng, ALPHA_SIDE_EST, n)
        zE = norm.ppf(1 - ALPHA_SIDE_E)
        r["clr_B_pop"] = (-Ln - zE * S1, U + zE * S1)
        U2 = clr_upper(thU, GU, V, rng, ALPHA / 2, n)
        L2 = clr_upper(-thL, -GL, V, rng, ALPHA / 2, n)
        r["clr_Bt"] = (-L2, U2)
        # Bhatia-Davis allowance (as the package): Var(E) <= sum w^2 mbar(1 - mbar), mbar bounded by CLR
        wpi = w ** 2 / np.sum(w ** 2)
        _, _, _, cp, gp, _ = lin22(len(beta), dfit, x, N, Dfit, w_override=wpi)
        GUp = cp[None, :] - gp * S[up_j] / A[up_j, None]; iUp = -gp * one[up_j] / A[up_j]
        GLp = cp[None, :] - gp * S[lo_j] / A[lo_j, None]; iLp = -gp * one[lo_j] / A[lo_j]
        mU = clr_upper(GUp @ beta + iUp, GUp, V, rng, ALPHA_V / 2, n)
        mLn = clr_upper(-(GLp @ beta + iLp), -GLp, V, rng, ALPHA_V / 2, n)
        S1bd = bd_sd(w, -mLn, mU)
        zE2 = norm.ppf(1 - ALPHA_SIDE_E_BD)
        r["clr_B_bd"] = (-Ln - zE2 * S1bd, U + zE2 * S1bd)
        # projection (P3.5(b2)(ii)), same candidates; for B as the package: est. part at 0.04, BD allowance
        q = rank_psd(V)
        kq = np.sqrt(chi2.ppf(1 - ALPHA, q))
        r["proj_Bt"] = (-proj_upper(-thL, -GL, V, kq), proj_upper(thU, GU, V, kq))
        kq4 = np.sqrt(chi2.ppf(1 - 2 * ALPHA_SIDE_EST, q))
        r["proj_B_clrbd"] = (-proj_upper(-thL, -GL, V, kq4) - zE2 * S1bd, proj_upper(thU, GU, V, kq4) + zE2 * S1bd)
        kv = np.sqrt(chi2.ppf(1 - ALPHA_V, q))       # projection interval for mbar, level ALPHA_V (P3.5(c1))
        S1bdp = bd_sd(w, -proj_upper(-(GLp @ beta + iLp), -GLp, V, kv), proj_upper(GUp @ beta + iUp, GUp, V, kv))
        r["proj_B"] = (-proj_upper(-thL, -GL, V, kq4) - zE2 * S1bdp, proj_upper(thU, GU, V, kq4) + zE2 * S1bdp)
        r["q"] = q
        if finite and model == "M1":
            r["fin_Bt"] = p38(beta, Z, t, S, A, one, c, gam, w, ALPHA, 0.0)
            r["fin_B"] = p38(beta, Z, t, S, A, one, c, gam, w, 2 * ALPHA_SIDE_EST, ALPHA2)
        if model == "M1" and D > 1:
            idx = np.arange(3, len(beta))
            dz = beta[idx]; Vz = V[np.ix_(idx, idx)]
            ev, Uv = np.linalg.eigh((Vz + Vz.T) / 2); keep = ev > 1e-10 * ev.max()
            W = np.sum((Uv[:, keep].T @ dz) ** 2 / ev[keep])
            r["pool_p"] = chi2.sf(W, keep.sum()); r["pool_df"] = int(keep.sum()); r["pool_W"] = W
        out[model] = r
    return out


def p38(beta, Z, t, S, A, one, c, gam, w, a_s, a_E):
    """Theorem P3.8 for 2x2 (one free column): Hoeffding radius h(f) = sqrt(0.5 log(1/delta) sum_i w_i(f)^2)."""
    L = np.linalg.solve(Z.T @ Z, Z.T)            # beta_hat = L t
    J = len(S); delta = a_s / (J + 2)
    lf = np.sqrt(0.5 * np.log(1 / delta))
    hS = lf * np.sqrt(np.sum((S @ L) ** 2, axis=1))
    hc = lf * np.sqrt(np.sum((c @ L) ** 2))
    lo, hi = kappa_range(S @ beta + hS, A, one)
    if lo > hi:
        return None
    tE = np.sqrt(0.5 * np.log(2 / a_E) * np.sum(w ** 2)) if a_E > 0 else 0.0
    return (c @ beta - hc + gam * lo - tE, c @ beta + hc + gam * hi + tE)


# ======================= Design B: 3x3 pooled =======================
def M3(x):
    xP, xO, xN = x.T
    return np.stack([np.stack([0.08 - 0.03 * xP, 0.80 + 0.10 * xP, 0.12 - 0.07 * xP], -1),
                     np.stack([0.10 + 0.02 * xO, 0.10 - 0.05 * xO, 0.80 + 0.03 * xO], -1),
                     np.stack([0.75 - 0.10 * xN, 0.15 + 0.05 * xP, 0.10 + 0.10 * xN - 0.05 * xP], -1)], 1)


PAIRS = [(0, 1), (0, 2), (1, 2)]
MONO = [(r, s) for r in range(3) for s in range(r, 3)]
TARGETS = [(0, 1), (0, 2), (0, 0), (2, 1)]
TARGET_LABELS = ["loyalty P->P", "defection P->O", "demobilisation P->abstain", "mobilisation N->P"]


def run_B(seed, p=1000, conc=40.0):
    rng = np.random.default_rng(seed)
    X = rng.dirichlet([4, 3, 2.5], p)
    M = M3(X)
    beta_u = np.stack([np.stack([rng.dirichlet(conc * M[i, r]) for r in range(3)]) for i in range(p)])
    T = np.einsum("ir,irc->ic", X, beta_u)
    N = rng.integers(300, 1500, p).astype(float)
    Z = np.column_stack([X[:, r] * X[:, s] for r, s in MONO])
    coef = np.linalg.lstsq(Z, T, rcond=None)[0]                 # 6 x 3
    V = hc3_cov(Z, T - Z @ coef)                               # stacked by column c (rank 12: adding up)
    kz = len(MONO)
    def Qcoef(c, r, s):
        g = np.zeros(3 * kz)
        j = MONO.index((min(r, s), max(r, s)))
        g[c * kz + j] = 1.0 if r == s else 0.5
        return g
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
    beta = coef.T.reshape(-1)
    q = rank_psd(V)
    kq = np.sqrt(chi2.ppf(1 - ALPHA, q)); kq4 = np.sqrt(chi2.ppf(1 - 2 * ALPHA_SIDE_EST, q))
    Lc = sqrt_psd(V)
    res = []
    for (r, c) in TARGETS:
        w = N * X[:, r] / np.sum(N * X[:, r])
        cvec = sum(np.sum(w * X[:, s]) * Qcoef(c, r, s) for s in range(3))
        gam = sum(np.sum(w * X[:, s]) * Kcoef(c, r, s) for s in range(3))
        truth = np.sum(w * beta_u[:, r, c])
        truth_cm = np.sum(w * M[:, r, c])
        lo_dd = np.maximum(0, (T[:, c] - (1 - X[:, r])) / X[:, r]); hi_dd = np.minimum(1, T[:, c] / X[:, r])
        out = dict(B=truth, Bt=truth_cm, dd=np.sum(w * (hi_dd - lo_dd)))
        def lp(svals, g=gam):
            kw = dict(A_ub=-A, b_ub=svals, bounds=[(None, None)] * 6, method="highs")
            lo, hi = linprog(g, **kw), linprog(-g, **kw)
            return None if lo.status or hi.status else (lo.fun, -hi.fun)
        b0 = lp(S @ beta)
        rr = {}
        rr["plug"] = None if b0 is None else (cvec @ beta + b0[0], cvec @ beta + b0[1])
        G = np.vstack([S, cvec])
        S1 = 0.5 * np.sqrt(np.sum(w ** 2))
        for key, a_est, aE in (("cons_B", ALPHA1, ALPHA2), ("cons_Bt", ALPHA, 0.0)):
            k, se = crit(G, V, rng, a_est, two_sided_rows=[len(S)])
            b1 = lp(S @ beta + k * se[:-1])
            zz = norm.ppf(1 - aE / 2) * S1 if aE > 0 else 0.0
            rr[key] = None if b1 is None else (cvec @ beta - k * se[-1] + b1[0] - zz,
                                               cvec @ beta + k * se[-1] + b1[1] + zz)
        draws = [beta] + [beta + Lc @ rng.standard_normal(len(V)) for _ in range(150)]
        def duals(objsign, g):
            found = {}
            for b_ in draws:
                res_ = linprog(-objsign * g, A_ub=-A, b_ub=S @ b_, bounds=[(None, None)] * 6, method="highs")
                if res_.status: continue
                lam = -res_.ineqlin.marginals
                found[tuple(np.round(lam, 9))] = lam
            return np.array(list(found.values()))
        LU, LL = duals(+1.0, gam), duals(-1.0, gam)
        GU = cvec[None, :] + LU @ S
        GL = cvec[None, :] - LL @ S
        thU, thL = GU @ beta, GL @ beta
        zE = norm.ppf(1 - ALPHA_SIDE_E)
        U = clr_upper(thU, GU, V, rng, ALPHA_SIDE_EST, p)
        Ln = clr_upper(-thL, -GL, V, rng, ALPHA_SIDE_EST, p)
        rr["clr_B_pop"] = (-Ln - zE * S1, U + zE * S1)
        U2 = clr_upper(thU, GU, V, rng, ALPHA / 2, p)
        L2 = clr_upper(-thL, -GL, V, rng, ALPHA / 2, p)
        rr["clr_Bt"] = (-L2, U2)
        wpi = w ** 2 / np.sum(w ** 2)
        cp = sum(np.sum(wpi * X[:, s]) * Qcoef(c, r, s) for s in range(3))
        gp = sum(np.sum(wpi * X[:, s]) * Kcoef(c, r, s) for s in range(3))
        LUp, LLp = duals(+1.0, gp), duals(-1.0, gp)
        GUp = cp[None, :] + LUp @ S; GLp = cp[None, :] - LLp @ S
        mU = clr_upper(GUp @ beta, GUp, V, rng, ALPHA_V / 2, p)
        mLn = clr_upper(-(GLp @ beta), -GLp, V, rng, ALPHA_V / 2, p)
        S1bd = bd_sd(w, -mLn, mU)
        zE2 = norm.ppf(1 - ALPHA_SIDE_E_BD)
        rr["clr_B_bd"] = (-Ln - zE2 * S1bd, U + zE2 * S1bd)
        rr["proj_Bt"] = (-proj_upper(-thL, -GL, V, kq), proj_upper(thU, GU, V, kq))
        rr["proj_B_clrbd"] = (-proj_upper(-thL, -GL, V, kq4) - zE2 * S1bd, proj_upper(thU, GU, V, kq4) + zE2 * S1bd)
        kv = np.sqrt(chi2.ppf(1 - ALPHA_V, q))
        S1bdp = bd_sd(w, -proj_upper(-(GLp @ beta), -GLp, V, kv), proj_upper(GUp @ beta, GUp, V, kv))
        rr["proj_B"] = (-proj_upper(-thL, -GL, V, kq4) - zE2 * S1bdp, proj_upper(thU, GU, V, kq4) + zE2 * S1bdp)
        rr["q"] = q
        out["M0"] = rr
        res.append(out)
    return res


# ======================= driver and summaries =======================
METHODS = [("plug", "B"), ("cons_B", "B"), ("clr_B_pop", "B"), ("clr_B_bd", "B"), ("proj_B", "B"), ("proj_B_clrbd", "B"),
           ("fin_B", "B"), ("cons_Bt", "Bt"), ("clr_Bt", "Bt"), ("proj_Bt", "Bt"), ("fin_Bt", "Bt")]


def _job(args):
    kind, seed, spec = args
    if kind == "22":
        simname, simargs, models, finite = spec
        return run22(seed, SIMS[simname], simargs, models, finite)
    return run_B(seed)


SIMS = {"A": sim_A, "C": sim_C, "E": sim_E}

# name -> (kind, seed base, spec)
DESIGNS = {
    "A_10x93":  ("22", 1000, ("A", (10, 93), ("M1",), True)),
    "A_93x30":  ("22", 11000, ("A", (93, 30), ("M1",), True)),
    "R2_30x10": ("22", 21000, ("A", (30, 10), ("M1",), True)),
    "R2_93x10": ("22", 31000, ("A", (93, 10), ("M1",), True)),
    "R2_200x10": ("22", 41000, ("A", (200, 10), ("M1",), True)),
    "C_boundary": ("22", 51000, ("C", (), ("M0",), False)),
    "E_slater":  ("22", 61000, ("E", (), ("M0",), False)),
    "S0_noshift": ("22", 71000, ("A", (20, 50, "none"), ("M0", "M1"), False)),
    "S1_affine":  ("22", 81000, ("A", (20, 50, "affine"), ("M0", "M1"), False)),
    "S2_nonaffine": ("22", 91000, ("A", (20, 50, "nonaffine"), ("M0", "M1"), False)),
    "B_3x3":    ("B", 101000, None),
}


def mcse(p, R):
    return np.sqrt(p * (1 - p) / R)


def clipw(iv):
    return min(iv[1], 1.0) - max(iv[0], 0.0)


def summarise_rows(rows, model):
    """rows: list of per-replication dicts (B, Bt, dd, model -> {method: interval})."""
    R = len(rows)
    out = {"reps": R, "dd_width": float(np.median([r["dd"] for r in rows]))}
    for meth, tgt in METHODS:
        if meth not in rows[0][model]:
            continue
        ivs = [r[model][meth] for r in rows]
        ivs = [None if iv is None or not iv[0] <= iv[1] else iv for iv in ivs]   # crossed ends: empty
        cov = np.array([iv is not None and iv[0] - 1e-12 <= r[tgt] <= iv[1] + 1e-12 for iv, r in zip(ivs, rows)])
        empty = np.mean([iv is None for iv in ivs])
        ws = [clipw(iv) for iv in ivs if iv is not None]
        p = float(cov.mean())
        out[meth] = dict(cover=p, se=float(mcse(p, R)), lb95=float(0.05 ** (1 / R)) if p == 1 else None,
                         empty=float(empty), width=float(np.median(ws)) if ws else float("nan"))
    if "pool_p" in rows[0][model]:
        pv = np.array([r[model]["pool_p"] for r in rows])
        rej = float(np.mean(pv < 0.05))
        out["pool"] = dict(reject=rej, se=float(mcse(rej, R)), df=rows[0][model]["pool_df"])
    out["q"] = int(np.median([r[model]["q"] for r in rows]))
    return out


def fmt_cov(m):
    s = f"{m['cover']:.3f} ({m['se']:.3f})"
    if m["lb95"] is not None:
        s += f" [>= {m['lb95']:.3f}]"
    return s


def print_summary(name, model, sm):
    print(f"  {name} / {model}: reps {sm['reps']}, DD width {sm['dd_width']:.3f}, q = rank Sigma = {sm['q']}")
    for meth, tgt in METHODS:
        if meth in sm:
            m = sm[meth]
            print(f"    {meth:12s} target {tgt:2s}: cover {fmt_cov(m)}, empty {m['empty']:.3f}, median width {m['width']:.3f}")
    if "pool" in sm:
        print(f"    pooling test (Wald, zeta = 0, df {sm['pool']['df']}): rejection rate at 5% {sm['pool']['reject']:.3f} ({sm['pool']['se']:.3f})")


if __name__ == "__main__":
    REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    names = sys.argv[2:] if len(sys.argv) > 2 else list(DESIGNS)
    procs = os.cpu_count() or 1
    print(f"phase3_montecarlo.py: {REPS} replications per design, {procs} processes. "
          f"Coverage: share of replications (MC standard error in parentheses; when every replication covers, "
          f"[>= b] is the exact one-sided 95% lower confidence bound 0.05^(1/R)). Widths: median, clipped to [0, 1].")
    summaries = {}
    RAW_DIR = os.environ.get("MC_RAW_DIR", "/tmp/phase3_mc_raw")
    os.makedirs(RAW_DIR, exist_ok=True)
    with Pool(procs) as pool:
        for name in names:
            kind, base, spec = DESIGNS[name]
            t0 = time.time()
            rows = pool.map(_job, [(kind, base + s, spec) for s in range(REPS)], chunksize=1)
            print(f"\nDesign {name} ({time.time() - t0:.0f} s)")
            with gzip.open(os.path.join(RAW_DIR, f"{name}.json.gz"), "wt") as f:
                json.dump(rows, f, default=lambda o: o.item() if hasattr(o, "item") else list(o))
            if kind == "B":
                for j, lab in enumerate(TARGET_LABELS):
                    sm = summarise_rows([r[j] for r in rows], "M0")
                    summaries[f"{name}:{lab}:M0"] = sm
                    print_summary(f"{name} {lab}", "M0", sm)
            else:
                for model in spec[2]:
                    sm = summarise_rows(rows, model)
                    summaries[f"{name}:{model}"] = sm
                    print_summary(name, model, sm)
            sys.stdout.flush()
    tag = ("_" + "_".join(names)) if len(sys.argv) > 2 else ""
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"phase3_montecarlo_summary{tag}.json"), "w") as f:
        json.dump(summaries, f, indent=1)
