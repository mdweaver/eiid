"""
Numerical checks of two simple structural claims in NEW_AVENUES.md.

(C1) Joint (cast-vote-record) margins add nothing to Duncan-Davis.  For one unit with
     rows r (shares x), contest-1 choice y1 and contest-2 choice y2, the sharp bounds on
     a cell n(r, y1) are the same whether the column information is the one-way
     margin t1 or the joint two-way margin t12 of (y1, y2).
(C2) A known aggregate table pins only (R-1)(C-1) directions of K.  The map
     K -> B*(K) (2a Lemma 2) has rank (R-1)(C-1), while dim K = (C-1)R(R-1)/2, so for
     R >= 3 the true table leaves (C-1)(R-1)(R-2)/2 directions of K free.
"""
import itertools
import numpy as np
from scipy.optimize import linprog
rng = np.random.default_rng(11)
out = ["# Structural checks", ""]

# ---- C1
def cell_bounds_oneway(x, t1, r, c):
    R, C = len(x), len(t1)
    # variables n[r,c] >= 0, row sums x, col sums t1
    A = []; b = []
    for i in range(R):
        row = np.zeros(R * C); row[i * C:(i + 1) * C] = 1; A.append(row); b.append(x[i])
    for j in range(C):
        col = np.zeros(R * C); col[j::C] = 1; A.append(col); b.append(t1[j])
    obj = np.zeros(R * C); obj[r * C + c] = 1
    lo = linprog(obj, A_eq=A, b_eq=b, bounds=(0, None)).fun
    hi = -linprog(-obj, A_eq=A, b_eq=b, bounds=(0, None)).fun
    return lo, hi

def cell_bounds_joint(x, t12, r, c):
    R, C1, C2 = len(x), *t12.shape
    idx = lambda i, j, k: (i * C1 + j) * C2 + k
    nv = R * C1 * C2; A = []; b = []
    for i in range(R):
        row = np.zeros(nv)
        for j, k in itertools.product(range(C1), range(C2)): row[idx(i, j, k)] = 1
        A.append(row); b.append(x[i])
    for j, k in itertools.product(range(C1), range(C2)):
        row = np.zeros(nv)
        for i in range(R): row[idx(i, j, k)] = 1
        A.append(row); b.append(t12[j, k])
    obj = np.zeros(nv)
    for k in range(C2): obj[idx(r, c, k)] = 1
    lo = linprog(obj, A_eq=A, b_eq=b, bounds=(0, None)).fun
    hi = -linprog(-obj, A_eq=A, b_eq=b, bounds=(0, None)).fun
    return lo, hi

maxdiff = 0.0
for _ in range(300):
    R, C1, C2 = rng.integers(2, 5), rng.integers(2, 5), rng.integers(2, 4)
    n = rng.dirichlet(np.ones(R * C1 * C2) * 0.5).reshape(R, C1, C2)
    x, t12 = n.sum((1, 2)), n.sum(0); t1 = t12.sum(1)
    r, c = rng.integers(R), rng.integers(C1)
    a = cell_bounds_oneway(x, t1, r, c); bb = cell_bounds_joint(x, t12, r, c)
    maxdiff = max(maxdiff, abs(a[0] - bb[0]), abs(a[1] - bb[1]))
out.append(f"C1: 300 random unit tables (R in 2..4, C1 in 2..4, C2 in 2..3): max |DD(one-way) - DD(joint)| = {maxdiff:.2e}  -> {'PASS' if maxdiff < 1e-7 else 'FAIL'}")

# ---- C2
def bstar_map(R, C, mu):
    """Linear part of K -> B*_{rc} (2a Lemma 2) in the parametrisation of 2a Prop. 7."""
    pairs = [(r, s) for r in range(R) for s in range(r + 1, R)]
    d = (C - 1) * len(pairs)
    M = np.zeros((R * C, d))
    for ci in range(C - 1):
        for pi, (r, s) in enumerate(pairs):
            col = ci * len(pairs) + pi
            # K_c has +k at (r,s), -k at (s,r); K_C = -sum K_c
            for cc, sign in ((ci, 1.0), (C - 1, -1.0)):
                M[r * C + cc, col] += sign * mu[r, s] / mu.sum(1)[r]
                M[s * C + cc, col] -= sign * mu[r, s] / mu.sum(1)[s]
    return M, d

ok = True
for R, C in [(2, 2), (2, 3), (3, 2), (3, 3), (4, 3), (4, 4), (5, 3)]:
    X = rng.dirichlet(np.ones(R), 500); N = rng.integers(100, 1000, 500)
    mu = np.einsum("i,ir,is->rs", N, X, X) / 500
    M, d = bstar_map(R, C, mu)
    rk = np.linalg.matrix_rank(M, tol=1e-10)
    pred = (R - 1) * (C - 1)
    ok &= (rk == min(pred, d))
    out.append(f"C2: R={R}, C={C}: dim K = {d}, rank of K -> B* = {rk}, predicted min(d,(R-1)(C-1)) = {min(pred, d)}; free directions left by the true table = {d - rk}")
out.append(f"C2: {'PASS' if ok else 'FAIL'}")
txt = "\n".join(out); print(txt)
open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
