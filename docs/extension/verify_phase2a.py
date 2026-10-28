# Numerical / symbolic checks of the Phase 2a results (PHASE2a_rxc.md).
# Each check prints PASS/FAIL. Run: python3 -I verify_phase2a.py
# Needs numpy, scipy, sympy (verification only; not package dependencies).
import itertools
import numpy as np
import sympy as sp
from scipy.optimize import linprog

rng = np.random.default_rng(2026)
ok = lambda name, cond: print(("PASS " if cond else "FAIL ") + name)


# ---------- helpers: parametrisation of K and the outer-set LP ----------
def pairs(R):
    return [(r, s) for r in range(R) for s in range(r + 1, R)]


def K_from_vec(k, R, C):
    """k has (C-1)*R(R-1)/2 entries; returns list of C skew matrices summing to 0."""
    Ks, P, i = [], pairs(R), 0
    for c in range(C - 1):
        K = np.zeros((R, R))
        for (r, s) in P:
            K[r, s], K[s, r] = k[i], -k[i]
            i += 1
        Ks.append(K)
    Ks.append(-sum(Ks))
    return Ks


def outer_constraints(Q, V):
    """A_ub k <= b_ub encoding (Q_c + K_c) v >= 0 for all c, rows r, v in V."""
    R, C = Q[0].shape[0], len(Q)
    nk = (C - 1) * len(pairs(R))
    A, b = [], []
    for v in V:
        for c in range(C):
            for r in range(R):
                # e_r^T (Q_c + K_c) v >= 0  <=>  -(e_r^T K_c v) <= e_r^T Q_c v
                row = np.zeros(nk)
                for j in range(nk):
                    e = np.zeros(nk); e[j] = 1
                    row[j] = (K_from_vec(e, R, C)[c] @ v)[r]
                A.append(-row); b.append((Q[c] @ v)[r])
    return np.array(A), np.array(b)


def lp_range(obj, A, b):
    lo = linprog(obj, A_ub=A, b_ub=b, bounds=[(None, None)] * len(obj), method="highs")
    hi = linprog(-obj, A_ub=A, b_ub=b, bounds=[(None, None)] * len(obj), method="highs")
    return lo.fun, -hi.fun


def random_H(R, C):
    """Random valid homogeneous coefficients: each H_c[r, s] >= 0, sum_c H_c = J."""
    G = rng.gamma(1.0, size=(C, R, R))
    return [G[c] / G.sum(0) for c in range(C)]


# ---------- Lemma 1: identification decomposition ----------
for R, C in [(2, 2), (2, 3), (3, 3), (4, 3)]:
    # linear map from free H entries (with sum_c H_c = J) to quadratic-form coefficients
    nfree = (C - 1) * R * R
    def coeffs(h):
        H = [h[c * R * R:(c + 1) * R * R].reshape(R, R) for c in range(C - 1)]
        H.append(np.ones((R, R)) - sum(H))
        out = []
        for Hc in H:
            Qc = (Hc + Hc.T) / 2
            out += [Qc[r, s] for r in range(R) for s in range(r, R)]
        return np.array(out)
    base = coeffs(np.zeros(nfree))
    M = np.array([coeffs(np.eye(nfree)[j]) - base for j in range(nfree)]).T
    null_dim = nfree - np.linalg.matrix_rank(M)
    ok(f"Lemma 1: dim of unidentified part for {R}x{C} = {null_dim} == (C-1)R(R-1)/2 = {(C-1)*R*(R-1)//2}",
       null_dim == (C - 1) * R * (R - 1) // 2)

# symbolic: x'Kx = 0 for skew K and quadratic form uniqueness on the simplex (R = 3)
x1, x2 = sp.symbols("x1 x2", positive=True)
x = sp.Matrix([x1, x2, 1 - x1 - x2])
Ksym = sp.Matrix(3, 3, lambda r, s: 0 if r == s else (sp.Symbol(f"k{min(r,s)}{max(r,s)}") * (1 if r < s else -1)))
ok("Lemma 1(b): x'Kx == 0 for skew K (symbolic, R=3)", sp.simplify((x.T * Ksym * x)[0]) == 0)
Dsym = sp.Matrix(3, 3, lambda r, s: sp.Symbol(f"d{min(r,s)}{max(r,s)}"))
poly = sp.Poly(sp.expand((x.T * Dsym * x)[0]), x1, x2)
sol = sp.solve(poly.coeffs(), list(Dsym.free_symbols), dict=True)
ok("Lemma 1(a): x'Dx == 0 on the simplex forces symmetric D = 0 (symbolic, R=3)",
   len(sol) == 1 and all(v == 0 for v in sol[0].values()))

# ---------- Proposition 4: JKST Prop. 3 (A = whole simplex) for 2x2 ----------
w0, c1, w1 = sp.symbols("w0 c1 w1")
Q12 = (c1 + 2 * w0) / 2            # Q_{1,12} (vote column)
K = c1 / 2 - w1                    # K_{1,12}
# Prop 5: |K_c| <= Q_c for both columns (vote: Q12, non-vote: 1 - Q12, K -> -K)
conds = [sp.simplify(Q12 - K), sp.simplify(Q12 + K), sp.simplify(1 - Q12 + K), sp.simplify(1 - Q12 - K)]
jkst = [c1 + w0 - w1, w1 + w0, w1 - (c1 + w0 - 1), 1 - w0 - w1]
ok("Prop 4/5: |K_c,12| <= Q_c,12 reproduces JKST Prop. 3 constraints",
   {sp.simplify(a) for a in conds} == {sp.simplify(b) for b in jkst})

# ---------- Proposition 5: closed form when A = simplex ----------
fails = 0
for _ in range(200):
    R, C = int(rng.integers(2, 5)), int(rng.integers(2, 5))
    H = random_H(R, C)
    Q = [(h + h.T) / 2 for h in H]
    A, b = outer_constraints(Q, list(np.eye(R)))
    for (r, s) in pairs(R):
        for c in range(C - 1):
            j = c * len(pairs(R)) + pairs(R).index((r, s))
            obj = np.zeros(A.shape[1]); obj[j] = 1
            lo, hi = lp_range(obj, A, b)
            m = min(Q[c][r, s], 1 - Q[c][r, s])
            fails += not (abs(lo + m) < 1e-8 and abs(hi - m) < 1e-8)
ok("Prop 5: K_c,rs in [-min(Q,1-Q), min(Q,1-Q)] matches LP (200 random tables)", fails == 0)


# ---------- Proposition 6: R = 2 closed form, projection and greedy LP ----------
def r2_intervals(Q, V):
    C = len(Q)
    lo = np.array([max(-(Q[c][0, 0] * v[0] + Q[c][0, 1] * v[1]) / v[1] for v in V) for c in range(C)])
    hi = np.array([min((Q[c][0, 1] * v[0] + Q[c][1, 1] * v[1]) / v[0] for v in V) for c in range(C)])
    return lo, hi


def greedy_max(g, lo, hi):
    """max sum g_c k_c s.t. lo <= k <= hi, sum k = 0 (assumes feasible)."""
    k = lo.copy(); slack = -lo.sum()
    for c in np.argsort(-g):
        add = min(hi[c] - lo[c], slack); k[c] += add; slack -= add
    return g @ k


fails = fails_g = 0
for _ in range(300):
    C = int(rng.integers(2, 6))
    H = random_H(2, C)
    Q = [(h + h.T) / 2 for h in H]
    l, u = np.sort(rng.uniform(0.02, 0.98, 2))
    V = [np.array([l, 1 - l]), np.array([u, 1 - u])]
    A, b = outer_constraints(Q, V)
    lo, hi = r2_intervals(Q, V)
    for c in range(C - 1):
        obj = np.zeros(C - 1); obj[c] = 1
        L, U = lp_range(obj, A, b)
        pl = max(lo[c], -(hi.sum() - hi[c])); pu = min(hi[c], -(lo.sum() - lo[c]))
        fails += not (abs(L - pl) < 1e-8 and abs(U - pu) < 1e-8)
    g = rng.normal(size=C)
    obj = g[:C - 1] - g[C - 1]                     # k_C = -sum others
    _, U = lp_range(obj, A, b)
    fails_g += abs(U - greedy_max(g, lo, hi)) > 1e-8
ok("Prop 6: projected interval formula matches LP (R=2, random C, A=[l,u])", fails == 0)
ok("Prop 6: greedy solution of the functional LP matches LP", fails_g == 0)

# ---------- Corollary 6.1: joint R x C vs collapsed (column c vs rest) ----------
def H_valid_on(V, C):
    """H whose implied rates are valid probabilities at the R affinely independent
    points V (hence on conv V), but may extrapolate outside [0,1] elsewhere."""
    Vm = np.array(V); R = len(V)
    Ms = [rng.dirichlet(np.ones(C), R) for _ in V]          # Ms[k][r, c] = rate at V[k]
    return [np.linalg.solve(Vm, np.array([Ms[k][:, c] for k in range(R)])).T for c in range(C)]


def k_range(Q, V, c, rs):
    R = Q[0].shape[0]; A, b = outer_constraints(Q, V)
    obj = np.zeros(A.shape[1]); obj[c * len(pairs(R)) + pairs(R).index(rs)] = 1
    return lp_range(obj, A, b)


for R, label in [(2, "2x3, H valid only on A"), (2, "2x3, H >= 0 everywhere"), (3, "3x3, H >= 0 everywhere")]:
    viol = strict = 0
    for _ in range(150):
        V = list(np.sort(rng.uniform(.05, .95, 2))[:, None] * np.array([1, -1]) + np.array([0, 1])) if R == 2 \
            else list(rng.dirichlet([1, 1, 1], 5))
        H = H_valid_on(V, 3) if "only" in label else random_H(R, 3)
        Q = [(h + h.T) / 2 for h in H]
        for rs in pairs(R):
            jl, ju = k_range(Q, V, 0, rs)
            cl, cu = k_range([Q[0], np.ones((R, R)) - Q[0]], V, 0, rs)
            viol += (jl < cl - 1e-8) or (ju > cu + 1e-8)
            strict += (jl > cl + 1e-6) or (ju < cu - 1e-6)
    ok(f"Cor 6.1 ({label}): joint never looser than collapsed; strictly tighter in {strict} cases", viol == 0)
    if label == "2x3, H >= 0 everywhere":
        ok("Cor 6.1: for R=2 with nonnegative pure-poll rates, joint == collapsed", strict == 0)

# ---------- Proposition 3 / Remark: Theta_X vs Theta_O in 2x2 with dispersed T|X ----------
def theta_X_2x2(tau_pts, xs):
    """JKST Prop. 1 bounds on w1 when T|X=x is two-point {tau-e, tau+e} (prob 1/2)."""
    lo, hi = -np.inf, np.inf
    for x, (w0, ts) in zip(xs, tau_pts):
        EL = np.mean([max(0, (t - x) / (1 - x)) for t in ts])
        EU = np.mean([min(1, t / (1 - x)) for t in ts])
        lo, hi = max(lo, (EL - w0) / x), min(hi, (EU - w0) / x)
    return lo, hi


# Random valid designs: tau(x) = w0 + c1 x + d1 x^2 in (0,1) on [l,u], T | X=x two-point tau +/- e.
xs = np.linspace(0.05, 0.95, 91)
n_ok = n_strict = n_bad = n_lin = n_feas = 0
while n_ok < 400:
    w0_, c1_, d1_ = rng.uniform(0, 1), rng.uniform(-1, 1), rng.uniform(-1, 1) * (rng.uniform() < 0.5)
    tau = w0_ + c1_ * xs + d1_ * xs**2
    e = rng.uniform(0, 0.2)
    if tau.min() - e < 0 or tau.max() + e > 1:
        continue                                          # not a valid dataset
    n_ok += 1
    lx, ux = theta_X_2x2([(w0_, (t - e, t + e)) for t in tau], xs)
    lo_o, hi_o = theta_X_2x2([(w0_, (t, t)) for t in tau], xs)    # Jensen at the mean = Theta_O
    n_bad += (lx < lo_o - 1e-10) or (ux > hi_o + 1e-10)
    n_strict += (lx > lo_o + 1e-6) or (ux < hi_o - 1e-6)
    if d1_ == 0:   # then beta_w = beta_b = T is a feasible structure with w1 = c1
        n_lin += 1
        n_feas += lx - 1e-12 <= c1_ <= ux + 1e-12
ok("Prop 1: Theta_X inside Theta_O in all 400 valid 2x2 designs", n_bad == 0)
ok(f"Prop 1: inclusion is strict in {n_strict}/400 designs (dispersion straddling a binding kink)", n_strict > 0)
ok(f"Thm 1 (=>): a known feasible structure lies in Theta_X ({n_feas}/{n_lin} linear designs)", n_feas == n_lin)


# ---------- 3x3 worked example (Section 9 of PHASE2a_rxc.md) ----------
# rows/cols: 0 = in-party (P), 1 = out-party (O), 2 = non-voter (N)
def M_true(x):
    """True conditional-mean transition matrix: linear contextual effects."""
    xP, xO, xN = x
    m = np.array([
        [0.08 - 0.03 * xP, 0.80 + 0.10 * xP, 0.12 - 0.07 * xP],   # P-voters
        [0.10 + 0.02 * xO, 0.10 - 0.05 * xO, 0.80 + 0.03 * xO],   # O-voters
        [0.75 - 0.10 * xN, 0.15 + 0.05 * xP, 0.10 + 0.10 * xN - 0.05 * xP],   # non-voters
    ])
    return m  # columns: 0 = abstain, 1 = vote P, 2 = vote O


def H_from_M():
    """Recover homogeneous H_c from the affine specification by evaluating at vertices e_s."""
    return [np.array([[M_true(np.eye(3)[s])[r, c] for s in range(3)] for r in range(3)]) for c in range(3)]


H = H_from_M()
assert all(np.allclose(np.array([[(H[c] @ v)[r] for c in range(3)] for r in range(3)]), M_true(v))
           for v in rng.dirichlet([2, 2, 2], 5))
p = 2000
X = rng.dirichlet([4, 3, 2.5], p)
Nn = rng.integers(300, 1500, p).astype(float)
Q = [(h + h.T) / 2 for h in H]
mu_r = (Nn[:, None] * X).sum(0)
mu_rs = np.einsum("i,ir,is->rs", Nn, X, X)
Bstar_true = np.array([[(H[c][r] @ mu_rs[r]) / mu_r[r] for c in range(3)] for r in range(3)])

# outer set with A = convex hull of the observed X (constraints at all points = at hull vertices)
from scipy.spatial import ConvexHull
hull = ConvexHull(X[:, :2])
V = list(X[hull.vertices])
A, b = outer_constraints(Q, V)
names = ["P", "O", "N"]; cols = ["abstain", "P", "O"]
print("     3x3 example (population Q, A = hull of 2000 simulated polls):")
for r, c, lab in [(0, 1, "loyalty P->P"), (0, 2, "defection P->O"), (0, 0, "demobilisation P->abstain"),
                  (2, 1, "mobilisation N->P")]:
    # B*_rc = (1/mu_r)[sum_s Q_c,rs mu_rs + sum_{s != r} K_c,rs mu_rs]; K_c,rs is +/- a coordinate of k
    const = (Q[c][r] @ mu_rs[r]) / mu_r[r]
    obj = np.zeros(A.shape[1])
    for j in range(A.shape[1]):
        e = np.zeros(A.shape[1]); e[j] = 1
        obj[j] = (K_from_vec(e, 3, 3)[c][r] @ mu_rs[r]) / mu_r[r]
    L, U = lp_range(obj, A, b)
    Tc = np.array([M_true(xi)[:, c] @ xi for xi in X])             # noise-free column shares
    dd = (np.sum(Nn * np.maximum(0, X[:, r] + Tc - 1)) / mu_r[r], np.sum(Nn * np.minimum(X[:, r], Tc)) / mu_r[r])
    print(f"       {lab:27s} truth {Bstar_true[r, c]:.3f}   outer [{const + L:.3f}, {const + U:.3f}]   DD [{dd[0]:.3f}, {dd[1]:.3f}]")
    ok(f"3x3 example: truth inside outer bound for {lab}", const + L - 1e-9 <= Bstar_true[r, c] <= const + U + 1e-9)
