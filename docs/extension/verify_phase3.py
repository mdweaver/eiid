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
