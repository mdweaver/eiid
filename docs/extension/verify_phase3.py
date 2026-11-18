# Deterministic and small-simulation checks of Phase 3 results (PHASE3_inference.md).
# Run: python3 -I verify_phase3.py
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

rng = np.random.default_rng(7)
ok = lambda name, cond: print(("PASS " if cond else "FAIL ") + name)

# ---------- shared 3x3 pieces (pooled outer set; same parametrisation as Phase 2a Prop. 7) ----------
PAIRS = [(0, 1), (0, 2), (1, 2)]


def Kcoef(c, r, s):
    a = np.zeros(6)
    if r == s: return a
    sign = 1.0 if r < s else -1.0
    j = PAIRS.index((min(r, s), max(r, s)))
    if c < 2: a[c * 3 + j] = sign
    else: a[j] = -sign; a[3 + j] = -sign
    return a


def H_true():
    def M(x):
        xP, xO, xN = x
        return np.array([[0.08 - 0.03 * xP, 0.80 + 0.10 * xP, 0.12 - 0.07 * xP],
                         [0.10 + 0.02 * xO, 0.10 - 0.05 * xO, 0.80 + 0.03 * xO],
                         [0.75 - 0.10 * xN, 0.15 + 0.05 * xP, 0.10 + 0.10 * xN - 0.05 * xP]])
    return [np.array([[M(np.eye(3)[s])[r, c] for s in range(3)] for r in range(3)]) for c in range(3)]


H = H_true(); Q = [(h + h.T) / 2 for h in H]
X = rng.dirichlet([4, 3, 2.5], 1500); Nn = rng.integers(300, 1500, 1500).astype(float)


def bounds(Qs, V, r, c, eta=0.0, X=X, N=Nn):
    """outer bounds for B_rc (sample weights), constraints (Q_c + K_c) v >= -eta at v in V,
    then widened by eta (Proposition S1)."""
    w = N * X[:, r] / np.sum(N * X[:, r])
    const = sum(np.sum(w * X[:, s]) * Qs[c][r, s] for s in range(3))
    gam = sum(np.sum(w * X[:, s]) * Kcoef(c, r, s) for s in range(3))
    Aub, bub = [], []
    for v in V:
        for cc in range(3):
            for rr in range(3):
                a = sum(Kcoef(cc, rr, s) * v[s] for s in range(3))
                Aub.append(-a); bub.append((Qs[cc] @ v)[rr] + eta)
    kw = dict(A_ub=np.array(Aub), b_ub=np.array(bub), bounds=[(None, None)] * 6, method="highs")
    lo, hi = linprog(gam, **kw), linprog(-gam, **kw)
    if lo.status or hi.status: return None
    return const + lo.fun - eta, const - hi.fun + eta


V_all = list(X[ConvexHull(X[:, :2]).vertices])

# Lemma P3.1 / Cor: the upper endpoint is concave and piecewise linear along any line in Q-space
Hdir = [rng.normal(0, 0.01, (3, 3)) for _ in range(3)]; Hdir = [(h + h.T) / 2 for h in Hdir]
Hdir[2] = -Hdir[0] - Hdir[1]                     # keep sum_c Q_c = J
ts = np.linspace(-1, 1, 41)
U = np.array([bounds([Q[c] + t * Hdir[c] for c in range(3)], V_all, 0, 1)[1] for t in ts])
second = U[2:] - 2 * U[1:-1] + U[:-2]
ok("Lemma P3.1: upper endpoint is concave along a line in Q (second differences <= 0)", np.all(second <= 1e-9))
ok("Lemma P3.1: ... and piecewise linear (most second differences exactly 0)", np.mean(np.abs(second) < 1e-9) > 0.7)

# Proposition P3.2(c): using the hull of a subsample (inner hull) widens the bounds (conservative)
sub = rng.choice(len(X), 200, replace=False)
V_sub = list(X[sub][ConvexHull(X[sub][:, :2]).vertices])
wider = True
for r, c in [(0, 1), (0, 2), (0, 0), (2, 1)]:
    b_all, b_sub = bounds(Q, V_all, r, c), bounds(Q, V_sub, r, c)
    wider &= b_sub[0] <= b_all[0] + 1e-12 and b_sub[1] >= b_all[1] - 1e-12
ok("Prop P3.2(c): sample hull inside the true hull => bounds weakly wider (4 targets)", wider)

# Theorem P3.4: plug-in endpoint biased toward the interior at ties, unbiased with a unique active constraint
def endpoint_bias(gap, reps=20000):
    """theta_U = min(theta1, theta2) with theta2 = theta1 + gap; estimates have N(0,1) noise (scaled)."""
    th = np.array([0.0, gap]); e = rng.standard_normal((reps, 2))
    return np.mean(np.min(th + e, 1)) - th.min()
b_tie, b_far = endpoint_bias(0.0), endpoint_bias(6.0)
ok(f"Thm P3.4: at a tie the plug-in endpoint is biased inward (mean error {b_tie:.3f} ~ -1/sqrt(pi) = {-1/np.sqrt(np.pi):.3f})",
   abs(b_tie + 1 / np.sqrt(np.pi)) < 0.02)
ok(f"Thm P3.4: with a unique active constraint the bias vanishes (mean error {b_far:.4f})", abs(b_far) < 0.02)

# Proposition S1: eta-sensitivity bounds are nested in eta and contain a worst-case undetectable violation
etas = [0.0, 0.01, 0.02, 0.05]
nested = True; prev = None
for eta in etas:
    b = bounds(Q, V_all, 0, 1, eta)
    if prev is not None: nested &= b[0] <= prev[0] + 1e-12 and b[1] >= prev[1] - 1e-12
    prev = b
ok("Prop S1: eta-bounds nested (wider as eta grows)", nested)
# undetectable violation (Phase 2a Prop. 9): Delta_rc(x) with x'Delta = 0, Delta 1 = 0, |Delta| <= eta
eta = 0.02
def Delta(x):
    """rows P, O: shift loyalty-type cells in opposite directions, weighted so x'Delta = 0;
    nonlinear in x (proportional to x_O)."""
    D = np.zeros((3, 3))
    b = eta * x[1]                       # |b| <= eta
    D[0, 1], D[0, 2] = +b * x[1], -b * x[1]           # row P: +to P, -to O  (row sums 0)
    D[1, 1], D[1, 2] = -b * x[0], +b * x[0]           # row O: offsets so that x'Delta = 0
    return D
assert all(np.allclose(X[i] @ Delta(X[i]), 0) and np.allclose(Delta(X[i]).sum(1), 0) for i in range(50))
w = Nn * X[:, 0] / np.sum(Nn * X[:, 0])
B_lin = sum(np.sum(w * X[:, s]) * H[1][0, s] for s in range(3))
B_viol = B_lin + np.sum(w * np.array([Delta(xi)[0, 1] for xi in X]))
b0, b_eta = bounds(Q, V_all, 0, 1, 0.0), bounds(Q, V_all, 0, 1, eta)
ok(f"Prop S1: violated-model truth {B_viol:.4f} lies in the eta={eta} bound [{b_eta[0]:.4f}, {b_eta[1]:.4f}]",
   b_eta[0] <= B_viol <= b_eta[1])

# Breakdown value for a claim (monotone in eta, found by bisection)
claim = lambda b: b is not None and b[0] > 0.30      # "loyalty is above 30%"
lo_e, hi_e = 0.0, 0.5
if claim(bounds(Q, V_all, 0, 1, 0.0)):
    for _ in range(40):
        mid = (lo_e + hi_e) / 2
        lo_e, hi_e = (mid, hi_e) if claim(bounds(Q, V_all, 0, 1, mid)) else (lo_e, mid)
    print(f"     breakdown eta* for 'loyalty > 0.30': {lo_e:.4f}")
    ok("Breakdown: claim holds just below eta* and fails just above",
       claim(bounds(Q, V_all, 0, 1, lo_e * 0.99)) and not claim(bounds(Q, V_all, 0, 1, hi_e * 1.01 + 1e-6)))

# ---------- Correction after review: Theorem P3.4(ii) needs distinct projected gradients ----------
reps = 200000
Sig = np.diag([1.0, 0.0])                                    # singular limiting covariance
Zs = rng.multivariate_normal(np.zeros(2), Sig, reps)
g1, g2 = np.array([1.0, 0.0]), np.array([1.0, 1.0])          # g1 != g2, but (g1-g2)' Sig (g1-g2) = 0
m_same = np.minimum(Zs @ g1, Zs @ g2).mean()
Zi = rng.standard_normal((reps, 2))
m_dist = np.minimum(Zi @ g1, Zi @ g2).mean()                 # Sig = I: Var(Y - Y') = 1, E min = -1/sqrt(2 pi)
ok(f"Thm P3.4(ii): two active gradients equal on supp(Z) give mean {m_same:.4f} (normal limit; the earlier "
   f"'non-singleton => negative mean' fails), distinct ones give {m_dist:.4f} ~ {-1/np.sqrt(2*np.pi):.4f}",
   abs(m_same) < 0.01 and abs(m_dist + 1 / np.sqrt(2 * np.pi)) < 0.01)

# ---------- Workstream 3: Proposition P3.5(b2) ----------
from scipy.stats import chi2, norm
# A small LP family: phi_U(b) = max{gam'k : b + A k >= 0}, k in R^2, m = 8 constraints (s(b) = b).
m_, A_ = 8, None
ang = np.linspace(0, 2 * np.pi, m_, endpoint=False)
A_ = -np.column_stack([np.cos(ang), np.sin(ang)])           # b_j - a_j'k >= 0: a bounded polygon
gam_ = np.array([1.0, 0.3])
b0 = np.full(m_, 1.0) + 0.02 * np.arange(m_)

def primal(b):
    r = linprog(-gam_, A_ub=-A_, b_ub=b, bounds=[(None, None)] * 2, method="highs")
    return -r.fun, -r.ineqlin.marginals                      # value, optimal dual (>= 0)

phiU0, _ = primal(b0)
# (1) weak duality: every dual-feasible lambda (A'lambda = -gam, lambda >= 0) gives lambda'b0 >= phi_U(b0)
viol = 0
for _ in range(500):
    r = linprog(rng.exponential(size=m_), A_eq=A_.T, b_eq=-gam_, bounds=[(0, None)] * m_, method="highs")
    lam = r.x
    viol += lam @ b0 < phiU0 - 1e-9
ok("Prop P3.5(b2) step 1: every dual-feasible lambda gives an upper bound at beta_0 (500 random duals)", viol == 0)

# (2) coverage of the projection version with a data-dependent candidate set (the plug-in optimal dual only)
n_, reps_ = 400, 3000
kq = np.sqrt(chi2.ppf(0.95, m_))
cover_proj = cover_naive = 0
for _ in range(reps_):
    bh = b0 + rng.standard_normal(m_) / np.sqrt(n_)
    val, lam = primal(bh)
    se = np.linalg.norm(lam) / np.sqrt(n_)
    cover_proj += phiU0 <= lam @ bh + kq * se
    cover_naive += phiU0 <= lam @ bh + norm.ppf(0.95) * se
ok(f"Prop P3.5(b2)(ii): projection upper bound with a data-dependent candidate covers {cover_proj / reps_:.3f} "
   f">= 0.95 (one-candidate z-bound: {cover_naive / reps_:.3f})", cover_proj / reps_ >= 0.95)

# ---------- Workstream 4: Theorem P3.8 (finite-sample CI), regime R2: M1, 2x2, 30 districts x 8 polls ----------
def p38_once(rng_, D=30, nd=8, a_s=0.05, a_E=0.01, target="Btilde"):
    Hr = np.array([[0.7, 0.5], [0.3, 0.2]])                    # column-1 rates: (H x)_r = Hr[r,0] x1 + Hr[r,1] x2
    zeta = np.vstack([np.zeros(2), rng_.uniform(-0.08, 0.08, (D - 1, 2))])
    X1, Dl, Nn_ = [], [], []
    for d in range(D):
        a = rng_.uniform(0.05, 0.5); b = min(0.95, a + rng_.uniform(0.2, 0.45))
        X1.append(rng_.uniform(a, b, nd)); Dl += [d] * nd; Nn_.append(rng_.integers(200, 1000, nd))
    X1 = np.concatenate(X1); X2 = 1 - X1; Dl = np.array(Dl); Nn_ = np.concatenate(Nn_).astype(float)
    m = np.column_stack([Hr[0, 0] * X1 + Hr[0, 1] * X2 + zeta[Dl, 0], Hr[1, 0] * X1 + Hr[1, 1] * X2 + zeta[Dl, 1]])
    conc = 15.0
    beta = rng_.beta(m * conc, (1 - m) * conc)
    T1 = X1 * beta[:, 0] + X2 * beta[:, 1]
    # M1 regression: pooled quadratic + district x X_r shifts (d >= 2)
    Zs = [X1 ** 2, 2 * X1 * X2, X2 ** 2]
    for d in range(1, D):
        Zs += [(Dl == d) * X1, (Dl == d) * X2]
    Z = np.column_stack(Zs); L = np.linalg.solve(Z.T @ Z, Z.T)  # beta_hat = L @ T1
    bh = L @ T1
    def zidx(d, r): return None if d == 0 else 3 + 2 * (d - 1) + r
    def rad(b, delta):                                          # Hoeffding radius of b'beta_hat (column 1 only)
        w = b @ L
        return np.sqrt(0.5 * np.log(1 / delta) * np.sum(w ** 2))
    cons = []                                                   # (b vector, f0, a): s = f0 + b'beta, constraint s + a kappa >= 0
    for d in range(D):
        xs_d = X1[Dl == d]
        for t in (xs_d.min(), xs_d.max()):
            b1 = np.zeros(Z.shape[1]); b1[0] = t; b1[1] = 1 - t
            if d: b1[zidx(d, 0)] = 1
            b2 = np.zeros(Z.shape[1]); b2[1] = t; b2[2] = 1 - t
            if d: b2[zidx(d, 1)] = 1
            cons += [(b1, 0.0, 1 - t), (b2, 0.0, -t), (-b1, 1.0, -(1 - t)), (-b2, 1.0, t)]
    J = len(cons); delta = a_s / (J + 2)
    lo, hi = -np.inf, np.inf
    for b, f0, a in cons:
        sh = f0 + b @ bh + rad(b, delta)
        if a > 0: lo = max(lo, -sh / a)
        else: hi = min(hi, sh / (-a))
    om = Nn_ * X1 / np.sum(Nn_ * X1)
    bc = np.zeros(Z.shape[1]); bc[0] = om @ X1; bc[1] = om @ X2
    for d in range(1, D): bc[zidx(d, 0)] = om[Dl == d].sum()
    gam_k = om @ X2
    hc = rad(bc, delta)
    Bt = om @ m[:, 0]; Br = om @ beta[:, 0]
    if lo > hi: return None, Bt, Br
    tE = np.sqrt(0.5 * np.log(2 / a_E) * np.sum(om ** 2)) if target == "B" else 0.0
    ci = (bc @ bh - hc + gam_k * lo - tE, bc @ bh + hc + gam_k * hi + tE)
    return ci, Bt, Br

r8 = np.random.default_rng(88)
cov_t = cov_b = empty = 0; wid = []
reps8 = 300
for _ in range(reps8):
    ci, Bt, Br = p38_once(r8, target="Btilde")
    if ci is None: empty += 1; continue
    cov_t += ci[0] <= Bt <= ci[1]; wid.append(ci[1] - ci[0])
for _ in range(reps8):
    ci2, Bt2, Br2 = p38_once(r8, a_s=0.04, a_E=0.01, target="B")
    if ci2 is not None: cov_b += ci2[0] <= Br2 <= ci2[1]
print(f"     Thm P3.8, M1 2x2, 30 districts x 8 polls: CI for B~ covers {cov_t / reps8:.3f}, for B {cov_b / reps8:.3f}; "
      f"empty {empty}; median unclipped width {np.median(wid):.2f} (valid but uninformative; see P3 section 5.1)")
ok("Thm P3.8: finite-sample CI covers B~ and B at >= 95% in an R2 design (300 reps each); emptiness rate <= 5%",
   cov_t / reps8 >= 0.95 and cov_b / reps8 >= 0.95 and empty / reps8 <= 0.05)

# Prop P3.6(c): with affinely spanning districts, district shifts that are constant within district cannot be
# written as G X_i unless they are all zero (conditional-on-design pooling test).
Xd = [rng.dirichlet(np.ones(3), 4) for _ in range(3)]
zs = [np.zeros(3), np.array([0.05, -0.02, 0.01]), np.array([-0.03, 0.04, 0.0])]
Xs_ = np.vstack(Xd); Zt = np.vstack([np.tile(z, (4, 1)) for z in zs])
G_, res_, *_ = np.linalg.lstsq(Xs_, Zt, rcond=None)
Zt0 = np.zeros_like(Zt); _, res0, *_ = np.linalg.lstsq(Xs_, Zt0, rcond=None)
ok("Prop P3.6(c): nonzero district shifts are not of the form G X_i (residual > 0); zero shifts are",
   res_.sum() > 1e-6 and (len(res0) == 0 or res0.sum() < 1e-20))

# ---------- Batch B review (Nov 2026): two corrections in P3.5 ----------
# (1) P3.5(b2)(ii) null direction. Error (Z/sqrt(n), -1/n): limit covariance diag(1, 0), rank 1. The old rule
#     (q = rank Sigma, pseudoinverse) gives the second coordinate standard error 0, so its upper bound
#     beta_hat_2 + 0 misses beta_0,2 in every sample. Under the corrected (ii-b) the null-direction error is
#     exactly zero (beta_hat = L eta_hat + b), and the same projection bound covers.
from scipy.stats import chi2 as _chi2, norm as _norm
n_ = 10_000; reps_ = 2000; kq = np.sqrt(_chi2.ppf(0.95, 1)); beta0 = np.array([0.3, 0.6])
Zs = rng.standard_normal(reps_)
Sig_hat = np.diag([1.0, 0.0])                                # estimated limit covariance, singular
g2 = np.array([0.0, 1.0]); se2 = np.sqrt(g2 @ Sig_hat @ g2 / n_)   # = 0
bh_old = beta0[None, :] + np.column_stack([Zs / np.sqrt(n_), np.full(reps_, -1.0 / n_)])
cov_old = np.mean(bh_old[:, 1] + kq * se2 >= beta0[1])
bh_new = beta0[None, :] + np.column_stack([Zs / np.sqrt(n_), np.zeros(reps_)])   # (ii-b): null error exactly 0
g1 = np.array([1.0, 0.0]); se1 = np.sqrt(g1 @ Sig_hat @ g1 / n_)
cov_new = np.mean((bh_new[:, 1] + kq * se2 >= beta0[1]) & (bh_new[:, 0] + kq * se1 >= beta0[0]))
print(f"     P3.5(b2)(ii) null direction: old rule covers coordinate 2 in {cov_old:.3f} of {reps_}; "
      f"under (ii-b) both coordinates covered in {cov_new:.3f}")
ok("P3.5(b2)(ii) null direction (Batch B): rank-q pseudoinverse rule misses surely (old statement fails); "
   "reduced full-rank parameterisation covers (corrected statement holds)", cov_old == 0.0 and cov_new >= 0.95)

# (2) P3.5(a) step 5, residual CLT. Weights w_i = 1/n (max w^2 / sum w^2 = 1/n -> 0) but only unit 1 has a
#     positive residual variance (e_1 = +-1/2). E / sd(E) is Rademacher, not normal: the old "CLT under the
#     weight condition" fails. The corrected step (Chebyshev branch, sd(E)/S_1 -> 0) still gives the bound.
n_r = 400; w_r = np.full(n_r, 1.0 / n_r); S1_r = 0.5 * np.sqrt(np.sum(w_r ** 2)); sd_r = 0.5 / n_r
E_r = w_r[0] * rng.choice([-0.5, 0.5], 20000)
z_r = _norm.ppf(0.975)
clt_tail = np.mean(np.abs(E_r / sd_r) > z_r)                 # would be ~0.05 if the old CLT claim held
bound_tail = np.mean(np.abs(E_r) > z_r * S1_r)
print(f"     P3.5(a) residual: P(|E/sd| > 1.96) = {clt_tail:.3f} (normal: 0.05); P(|E| > 1.96 S_1) = {bound_tail:.3f}; "
      f"sd/S_1 = {sd_r / S1_r:.3f}")
ok("P3.5(a) step 5 (Batch B): weight condition alone does not give a standardised CLT (old step fails); "
   "the allowance z*S_1 still holds (Chebyshev branch)", clt_tail == 0.0 and bound_tail <= 0.05)
