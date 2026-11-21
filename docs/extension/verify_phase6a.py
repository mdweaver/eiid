"""Checks for PHASE6a_count_layer.md (M-COUNT: the count layer; proof items R19, R27).

Each check prints PASS/FAIL, its tag and the numbers it relies on. Tags follow 2c/5a:
  [F] fails on an earlier (or tempting) statement and passes on the statement of 6a;
  [V] verifies a statement of 6a;
  [I] illustrates.
Run:  python3 -I verify_phase6a.py     (numpy, scipy >= 1.9 for milp)
"""
import itertools, math
from fractions import Fraction as Fr
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.stats import chi2, norm

rng = np.random.default_rng(20261121)
TOL = 1e-7
results = []
WITNESS_LOG = {"n": 0, "bad": 0}


def report(label, tag, ok, msg):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'} {label} [{tag}]: {msg}")


# ---------------------------------------------------------------------------
# The count layer: unknown cells n[i,r,c] >= 0 with unit row sums a[i,r] and
# column sums b[i,c]; extra rows; bounds per cell. Witness = the LP solution.
# ---------------------------------------------------------------------------
class CountLP:
    def __init__(self, tables=None, a=None, b=None):
        if tables is not None:
            a = [np.asarray(t).sum(1) for t in tables]; b = [np.asarray(t).sum(0) for t in tables]
        self.a = [np.asarray(x, float) for x in a]; self.b = [np.asarray(x, float) for x in b]
        self.I, self.R, self.C = len(self.a), len(self.a[0]), len(self.b[0])
        self.nv = self.I * self.R * self.C
        self.Aeq, self.beq, self.Aub, self.bub = [], [], [], []
        self.lb = np.zeros(self.nv); self.ub = np.full(self.nv, np.inf)
        for i in range(self.I):
            for r in range(self.R):
                row = np.zeros(self.nv); row[[self.ix(i, r, c) for c in range(self.C)]] = 1
                self.Aeq.append(row); self.beq.append(self.a[i][r])
            for c in range(self.C):
                row = np.zeros(self.nv); row[[self.ix(i, r, c) for r in range(self.R)]] = 1
                self.Aeq.append(row); self.beq.append(self.b[i][c])

    def ix(self, i, r, c):
        return (i * self.R + r) * self.C + c

    def copy(self):
        o = CountLP(a=self.a, b=self.b)
        o.Aeq, o.beq = list(self.Aeq), list(self.beq)
        o.Aub, o.bub = list(self.Aub), list(self.bub)
        o.lb, o.ub = self.lb.copy(), self.ub.copy()
        return o

    def ge(self, row, rhs): self.Aub.append(-np.asarray(row, float)); self.bub.append(-float(rhs))
    def le(self, row, rhs): self.Aub.append(np.asarray(row, float)); self.bub.append(float(rhs))
    def eq(self, row, rhs): self.Aeq.append(np.asarray(row, float)); self.beq.append(float(rhs))

    def cellrow(self, units, r, c):
        row = np.zeros(self.nv); row[[self.ix(i, r, c) for i in units]] = 1; return row

    def total(self, units, r, c, lo=None, hi=None):
        row = self.cellrow(units, r, c)
        if lo is not None and hi is not None and lo == hi: self.eq(row, lo); return
        if lo is not None: self.ge(row, lo)
        if hi is not None: self.le(row, hi)

    def box(self, i, r, c, lo=None, hi=None):
        k = self.ix(i, r, c)
        if lo is not None: self.lb[k] = max(self.lb[k], lo)
        if hi is not None: self.ub[k] = min(self.ub[k], hi)

    def feasible_point(self, x, tol=1e-6):
        """every margin, row and bound holds at x (the witness check)"""
        ok = np.all(np.abs(np.array(self.Aeq) @ x - np.array(self.beq)) <= tol)
        if self.Aub: ok &= np.all(np.array(self.Aub) @ x <= np.array(self.bub) + tol)
        ok &= np.all(x >= self.lb - tol) and np.all(x <= self.ub + tol)
        return bool(ok)

    def _bounds(self):
        return list(zip(self.lb, [None if np.isinf(u) else u for u in self.ub]))

    def lp(self, c, sense=1):
        """sense=+1: min c.n ; sense=-1: max c.n. Dual simplex, so the solution is a vertex."""
        kw = dict(A_eq=np.array(self.Aeq), b_eq=np.array(self.beq), bounds=self._bounds(), method="highs-ds")
        if self.Aub: kw.update(A_ub=np.array(self.Aub), b_ub=np.array(self.bub))
        r = linprog(sense * np.asarray(c, float), **kw)
        if r.status != 0: return None, None
        WITNESS_LOG["n"] += 1; WITNESS_LOG["bad"] += (not self.feasible_point(r.x))
        return sense * r.fun, r.x

    def ip(self, c, sense=1):
        cons = [LinearConstraint(np.array(self.Aeq), np.array(self.beq), np.array(self.beq))]
        if self.Aub: cons.append(LinearConstraint(np.array(self.Aub), -np.inf, np.array(self.bub)))
        m = milp(sense * np.asarray(c, float), constraints=cons, integrality=np.ones(self.nv),
                 bounds=Bounds(self.lb, self.ub))
        if m.status != 0: return None, None
        WITNESS_LOG["n"] += 1; WITNESS_LOG["bad"] += (not self.feasible_point(m.x))
        return sense * m.fun, m.x

    def rounded(self):
        """Round every row and bound to its integer hull where that is valid for integer
        tables: bounds to ceil/floor; rows with rational coefficients scaled to integer
        coefficients, then rhs rounded (Chvatal); equality rows with a non-integer
        right-hand side after scaling make the integer problem infeasible."""
        o = self.copy()
        o.lb = np.ceil(o.lb - 1e-9); o.ub = np.where(np.isinf(o.ub), o.ub, np.floor(o.ub + 1e-9))
        Aub, bub = [], []
        for row, rhs in zip(o.Aub, o.bub):      # row.n <= rhs
            ci, ri = integerise(row, rhs)
            Aub.append(ci); bub.append(math.floor(ri + 1e-9))
        o.Aub, o.bub = Aub, bub
        return o


def integerise(row, rhs):
    fr = [Fr(float(x)).limit_denominator(10 ** 6) for x in row]
    L = 1
    for x in fr: L = L * x.denominator // math.gcd(L, x.denominator)
    ci = [int(x * L) for x in fr]; g = 0
    for x in ci: g = math.gcd(g, abs(x))
    g = g or 1
    return np.array(ci, float) / g, float(Fr(float(rhs)).limit_denominator(10 ** 9) * L / g)


def frac(x):
    return 0.0 if x is None else float(np.max(np.abs(x - np.round(x))))


def tables_with_margins(a, b):
    """all nonnegative integer R x C tables with row sums a and column sums b"""
    R, C = len(a), len(b)
    out = []
    def rec(k, cur, rr, cc):
        if k == R * C:
            if all(v == 0 for v in rr) and all(v == 0 for v in cc): out.append(np.array(cur).reshape(R, C))
            return
        r, c = divmod(k, C)
        hi = min(rr[r], cc[c])
        lo = 0
        if c == C - 1: lo = hi = rr[r] if rr[r] <= cc[c] else -1
        if r == R - 1: lo = max(lo, cc[c]); hi = min(hi, cc[c])
        for v in range(max(lo, 0), hi + 1):
            rr2 = list(rr); cc2 = list(cc); rr2[r] -= v; cc2[c] -= v
            rec(k + 1, cur + [v], rr2, cc2)
    rec(0, [], list(a), list(b))
    return out


def dd_count(N, V, Y):
    return max(0, V + Y - N), min(V, Y)


print("=== Section 2: R19 (bound, DD, integrality, witnesses) ===")

# ---------------------------------------------------------------------------
# D1: margins only, single-cell rate targets: LP endpoint = Duncan-Davis (Frechet)
# ---------------------------------------------------------------------------
maxd = 0.0; ncmp = 0
for _ in range(150):
    I, R, C = rng.integers(1, 5), rng.integers(2, 5), rng.integers(2, 5)
    tabs = [rng.integers(0, 12, size=(R, C)) for _ in range(I)]
    L = CountLP(tabs)
    for r in range(R):
        A_r = sum(L.a[i][r] for i in range(I))
        if A_r == 0: continue
        for c in range(C):
            tgt = L.cellrow(range(I), r, c) / A_r
            lo, _ = L.lp(tgt, 1); hi, _ = L.lp(tgt, -1)
            ddlo = sum(max(0, L.a[i][r] + L.b[i][c] - L.a[i].sum()) for i in range(I)) / A_r
            ddhi = sum(min(L.a[i][r], L.b[i][c]) for i in range(I)) / A_r
            maxd = max(maxd, abs(lo - ddlo), abs(hi - ddhi)); ncmp += 1
report("D1 DD = margins-only LP, single cells (integer data)", "V", maxd < 1e-9,
       f"{ncmp} cell targets in random R x C (R,C = 2..4, 1-4 units): max |LP - Frechet| = {maxd:.1e}")


# the package's reference datasets (shares and sizes; counts N*X are not integers)
def load_csv(path):
    import csv
    with open(path) as f:
        rd = list(csv.DictReader(f))
    return rd


def package_dd(X, T, N, r, c):
    den = np.sum(N * X[:, r])
    return np.sum(N * np.maximum(0, X[:, r] + T[:, c] - 1)) / den, np.sum(N * np.minimum(X[:, r], T[:, c])) / den


import os
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "package", "eiPartialID", "tests", "testthat")
maxd = 0.0; ncmp = 0; names = []
for fn, R, C in (("rxc_example_3x3.csv", 3, 3), ("msign_example_3x2.csv", 3, 2)):
    p = os.path.join(base, fn)
    if not os.path.exists(p): continue
    rd = load_csv(p); names.append(f"{fn} ({len(rd)} units)")
    X = np.array([[float(row[f"X{r+1}"]) for r in range(R)] for row in rd])
    T = np.array([[float(row[f"T{c+1}"]) for c in range(C)] for row in rd])
    N = np.array([float(row["N"]) for row in rd])
    X = X / X.sum(1, keepdims=True); T = T / T.sum(1, keepdims=True)   # exact adding-up
    a = [N[i] * X[i] for i in range(len(N))]; b = [N[i] * T[i] for i in range(len(N))]
    for r in range(R):
        for c in range(C):
            # the LP separates over units; solve each unit's transportation LP, then sum
            lo = hi = 0.0
            for i in range(len(N)):
                Li = CountLP(a=[a[i]], b=[b[i]]); t = Li.cellrow([0], r, c)
                lo += Li.lp(t, 1)[0]; hi += Li.lp(t, -1)[0]
            den = np.sum(N * X[:, r]); plo, phi = package_dd(X, T, N, r, c)
            maxd = max(maxd, abs(lo / den - plo), abs(hi / den - phi)); ncmp += 1
report("D1b DD = margins-only LP on the package reference datasets", "V", ncmp > 0 and maxd < 1e-9,
       f"{', '.join(names)}: {ncmp} cell targets, max |LP - package formula| = {maxd:.1e} "
       "(non-integer counts N*X: continuous tables; exact for the stated margins, R19(e))")

# ---------------------------------------------------------------------------
# D2: within one column, a contrast of two rows is exact by cellwise interval arithmetic
# D3: [F] targets with several cells of one row, or same-sign cells of one column, are not
# ---------------------------------------------------------------------------
maxd = 0.0
for _ in range(150):
    I, R, C = rng.integers(1, 5), rng.integers(3, 5), rng.integers(2, 4)
    tabs = [rng.integers(1, 10, size=(R, C)) for _ in range(I)]
    L = CountLP(tabs); r, s = rng.choice(R, 2, replace=False); c = rng.integers(C)
    Ar = sum(x[r] for x in L.a); As = sum(x[s] for x in L.a)
    tgt = L.cellrow(range(I), r, c) / Ar - L.cellrow(range(I), s, c) / As
    lo, _ = L.lp(tgt, 1); hi, _ = L.lp(tgt, -1)
    lor = sum(max(0, L.a[i][r] + L.b[i][c] - L.a[i].sum()) for i in range(I)) / Ar
    hir = sum(min(L.a[i][r], L.b[i][c]) for i in range(I)) / Ar
    los = sum(max(0, L.a[i][s] + L.b[i][c] - L.a[i].sum()) for i in range(I)) / As
    his = sum(min(L.a[i][s], L.b[i][c]) for i in range(I)) / As
    maxd = max(maxd, abs(lo - (lor - his)), abs(hi - (hir - los)))
report("D2 one-column gap B_rc - B_sc: LP = cellwise interval arithmetic", "V", maxd < 1e-9,
       f"150 random R x C designs (R = 3, 4): max |diff| = {maxd:.1e}")

# merged-group rate (rows 1+2, one column) and a row total over two columns
tab = [np.array([[3, 2, 0], [2, 3, 0], [0, 0, 5]])]      # one unit, N = 15
L = CountLP(tab)
t1 = (L.cellrow([0], 0, 0) + L.cellrow([0], 1, 0)) / (L.a[0][0] + L.a[0][1])
lo1, _ = L.lp(t1, 1); hi1, _ = L.lp(t1, -1)
cw = lambda r, c: (max(0, L.a[0][r] + L.b[0][c] - 15), min(L.a[0][r], L.b[0][c]))
cwlo1 = (cw(0, 0)[0] + cw(1, 0)[0]) / 10; cwhi1 = (cw(0, 0)[1] + cw(1, 0)[1]) / 10
t2 = (L.cellrow([0], 0, 0) + L.cellrow([0], 0, 1)) / L.a[0][0]
lo2, _ = L.lp(t2, 1); hi2, _ = L.lp(t2, -1)
cwlo2 = (cw(0, 0)[0] + cw(0, 1)[0]) / 5; cwhi2 = (cw(0, 0)[1] + cw(0, 1)[1]) / 5
report("D3 multi-cell targets: cellwise DD is outer, the LP is exact", "F",
       (cwhi1 - cwlo1) > (hi1 - lo1) + 1e-9 and (cwhi2 - cwlo2) > (hi2 - lo2) + 1e-9,
       f"rows 1+2 merged, column 1: LP [{lo1:.3f},{hi1:.3f}] vs cellwise [{cwlo1:.3f},{cwhi1:.3f}]; "
       f"row 1, columns 1+2: LP [{lo2:.3f},{hi2:.3f}] vs cellwise [{cwlo2:.3f},{cwhi2:.3f}] "
       "(fails on 'the DD column is the exact count bound for every linear target')")

# ---------------------------------------------------------------------------
# B1: brute force over all integer tables in tiny designs, TU-class rows
# ---------------------------------------------------------------------------
nchk = 0; nbad = 0; nfrac = 0
for trial in range(60):
    I = rng.integers(2, 4); R, C = (2, 2) if trial % 3 else (2, 3) if trial % 2 else (3, 2)
    tabs = [rng.integers(0, 4, size=(R, C)) for _ in range(I)]
    L = CountLP(tabs)
    r0, c0 = rng.integers(R), rng.integers(C)
    # rows of TU class (i): a total on (r0,c0) over all units (lower), one over units {0,1},
    # a fractional box on another cell, a known cell.
    tot_all = sum(t[r0, c0] for t in tabs) - rng.integers(0, 2) + 0.4
    L.total(range(I), r0, c0, lo=tot_all)
    t01 = tabs[0][r0, c0] + tabs[1][r0, c0]; L.total([0, 1], r0, c0, hi=t01 + 0.7)
    r1, c1 = rng.integers(R), rng.integers(C); L.box(I - 1, r1, c1, lo=max(0, tabs[-1][r1, c1] - 1.5))
    L.box(0, R - 1, C - 1, lo=tabs[0][R - 1, C - 1], hi=tabs[0][R - 1, C - 1])
    Lr = L.rounded()
    tgt = rng.integers(-2, 3, size=L.nv)
    # enumerate
    per_unit = [tables_with_margins(t.sum(1), t.sum(0)) for t in tabs]
    vals = []
    for combo in itertools.product(*per_unit):
        x = np.concatenate([m.ravel() for m in combo]).astype(float)
        if L.feasible_point(x, tol=1e-9): vals.append(float(tgt @ x))
    for sense in (1, -1):
        v, x = Lr.lp(tgt, sense); nchk += 1
        best = (min(vals) if sense == 1 else max(vals)) if vals else None
        if (v is None) != (best is None) or (v is not None and abs(v - best) > 1e-7): nbad += 1
        nfrac += frac(x) > 1e-7
report("B1 brute force: integer endpoints = LP over the rounded TU-class system", "V", nbad == 0 and nfrac == 0,
       f"60 tiny designs (2-3 units of 2x2, 2x3, 3x2; all integer tables enumerated), {nchk} endpoints: "
       f"{nbad} mismatches, {nfrac} fractional witnesses")

# B2: every integer value of a TU-compatible aggregate between its extremes is attained
nmiss = 0; ntot = 0
for trial in range(40):
    I = rng.integers(2, 4); tabs = [rng.integers(0, 4, size=(2, 3)) for _ in range(I)]
    L = CountLP(tabs); L.box(0, 1, 2, lo=tabs[0][1, 2])
    per_unit = [tables_with_margins(t.sum(1), t.sum(0)) for t in tabs]
    agg = set()
    for combo in itertools.product(*per_unit):
        x = np.concatenate([m.ravel() for m in combo]).astype(float)
        if L.feasible_point(x, 1e-9): agg.add(int(round(sum(m[0, 0] for m in combo))))
    lo, hi = min(agg), max(agg); ntot += 1
    nmiss += len(agg) != hi - lo + 1
report("B2 the integer values of an aggregate cell total form an integer interval", "V", nmiss == 0,
       f"{ntot} designs (2-3 units of 2x3 with a box row): {nmiss} with a gap")

# ---------------------------------------------------------------------------
# T1-T3: TU classes (i)-(iii) on random instances (LP vertex integrality and = MILP)
# ---------------------------------------------------------------------------
def laminar_partition(rng, units):
    units = list(units); perm = rng.permutation(units)
    k = rng.integers(0, max(1, len(units) - 1))
    cuts = sorted(rng.choice(range(1, len(units)), size=min(k, len(units) - 1), replace=False)) if len(units) > 1 else []
    blocks = [list(b) for b in np.split(perm, cuts)]
    # refine one block into a nested sub-block (keeps the family laminar)
    out = blocks + [b[: max(1, len(b) // 2)] for b in blocks if len(b) > 1]
    return out


def tu_trial(cls):
    I = rng.integers(2, 6)
    if cls == "iii": R, C = 2, rng.integers(2, 5)
    else: R, C = rng.integers(2, 4), rng.integers(2, 4)
    tabs = [rng.integers(0, 7, size=(R, C)) for _ in range(I)]
    L = CountLP(tabs)
    def add_family(fam, r, c):
        for S in fam:
            tv = sum(tabs[i][r, c] for i in S)
            L.total(S, r, c, lo=tv - rng.integers(0, 2) - 0.3 * rng.random(), hi=tv + rng.integers(0, 2) + 0.3 * rng.random())
    if cls == "i":
        r0, c0 = rng.integers(R), rng.integers(C)
        add_family(laminar_partition(rng, range(I)) + [list(range(I))], r0, c0)
        add_family(laminar_partition(rng, range(I)), r0, c0)
    elif cls == "ii":
        c0 = rng.integers(C)
        for r in range(R):
            add_family(laminar_partition(rng, range(I)) + [list(range(I))], r, c0)
    else:
        for c in range(C):
            fam = laminar_partition(rng, range(I)) + [list(range(I))]
            for S in fam:
                r = rng.integers(2); tv = sum(tabs[i][r, c] for i in S)
                L.total(S, r, c, lo=tv - rng.integers(0, 2) - 0.3 * rng.random(), hi=tv + rng.integers(0, 2))
    for _ in range(2):   # boxes (fractional) and a known cell
        i, r, c = rng.integers(I), rng.integers(R), rng.integers(C)
        L.box(i, r, c, lo=max(0, tabs[i][r, c] - 1.5 * rng.random()), hi=tabs[i][r, c] + 1.5 * rng.random())
    i, r, c = rng.integers(I), rng.integers(R), rng.integers(C); L.box(i, r, c, lo=tabs[i][r, c], hi=tabs[i][r, c])
    Lr = L.rounded(); tgt = rng.integers(-3, 4, size=L.nv); bad = 0; fr_ = 0
    for sense in (1, -1):
        v, x = Lr.lp(tgt, sense); vi, _ = Lr.ip(tgt, sense)
        if v is None or vi is None or abs(v - vi) > 1e-7: bad += 1
        fr_ += frac(x) > 1e-7
    return bad, fr_


for cls, lab, txt in (("i", "T1", "class (i): totals on one cell over two laminar families, boxes, known cells"),
                      ("ii", "T2", "class (ii): totals on the cells of one column, one laminar family per cell"),
                      ("iii", "T3", "class (iii): R = 2, totals on any cells, one laminar family per column")):
    bad = fr_ = 0; n = 400
    for _ in range(n):
        b_, f_ = tu_trial(cls); bad += b_; fr_ += f_
    report(f"{lab} TU {txt}", "V", bad == 0 and fr_ == 0,
           f"{n} random designs, {2*n} endpoints: {bad} LP != MILP, {fr_} fractional dual-simplex vertices")


# exhaustive determinant check of total unimodularity on small constraint matrices
def is_tu(M, maxk=None):
    m, n = M.shape; maxk = maxk or min(m, n)
    for k in range(1, maxk + 1):
        for rows in itertools.combinations(range(m), k):
            sub = M[list(rows)]
            for cols in itertools.combinations(range(n), k):
                d = round(np.linalg.det(sub[:, cols]))
                if d not in (-1, 0, 1): return False, (rows, cols, d)
    return True, None


L2u = CountLP(a=[[1, 1], [1, 1]], b=[[1, 1], [1, 1]])
M_i = np.array(L2u.Aeq + [L2u.cellrow([0, 1], 0, 0), L2u.cellrow([0], 0, 0), L2u.cellrow([1], 0, 0)])  # A: {01},{0}; B: {1}
okA, _ = is_tu(M_i)
L = CountLP(a=[[1, 1]] * 3, b=[[1, 1]] * 3)
M_x = np.array(L.Aeq + [L.cellrow([0, 1], 0, 0), L.cellrow([1, 2], 0, 0), L.cellrow([0, 2], 0, 0)])  # three crossing sets
okX, witX = is_tu(M_x, maxk=3)
report("T4 exhaustive determinants: a class-(i) matrix is TU; three pairwise-crossing totals are not", "V",
       okA and not okX,
       f"2 units 2x2 + totals on one cell over {{01}},{{0}} (family A) and {{1}} (family B): all square minors in "
       f"{{-1,0,1}} = {okA}; 3 units + totals over {{01}},{{12}},{{02}}: TU = {okX} (a 3x3 minor = {witX[2] if witX else '-'})")

# ---------------------------------------------------------------------------
# T5 [F]: rounding the LP endpoint is not enough; round the rows, then solve
# ---------------------------------------------------------------------------
# two 2x2 units N=5, V=1, yes=1: the realised ordering gives n >= 1/5 in each unit
L = CountLP(a=[[1, 4], [1, 4]], b=[[1, 4], [1, 4]])
for i in range(2): L.box(i, 0, 0, lo=1 * 1 / 5)
tgt = L.cellrow([0, 1], 0, 0)
v_lp, _ = L.lp(tgt, 1); v_r, xr = L.rounded().lp(tgt, 1); v_ip, _ = L.ip(tgt, 1)
report("T5 round the rows, not the endpoint", "F",
       abs(v_lp - 0.4) < 1e-9 and math.ceil(v_lp - 1e-9) == 1 and abs(v_r - 2) < 1e-9 and abs(v_ip - 2) < 1e-9,
       f"two units (N=5, V=1, yes=1) with the realised ordering: LP lower total {v_lp:.1f}, its ceiling "
       f"{math.ceil(v_lp-1e-9)}; LP over the rounded rows {v_r:.0f} = integer {v_ip:.0f} "
       "(fails on 'the integer endpoint is the LP endpoint rounded')")

# ---------------------------------------------------------------------------
# N10a-c (from verify_next_steps.py, reproduced; N10c with explicit, named rows)
# ---------------------------------------------------------------------------
nonint = 0
for _ in range(300):
    R, C = rng.integers(2, 4), rng.integers(2, 5)
    n = rng.integers(0, 30, size=(R, C)); L = CountLP([n])
    v, x = L.lp(rng.integers(-3, 4, size=R * C), -1); nonint += frac(x) > 1e-7
report("N10a margins only: vertices integer", "V", nonint == 0,
       f"300 random integer tables (2-3 x 2-4), random integer targets: {nonint} fractional vertices")

L = CountLP(a=[[1, 4]], b=[[1, 4]]); L.box(0, 0, 0, lo=1 * 1 / 5)
v_lp, _ = L.lp(L.cellrow([0], 0, 0), 1); v_r, _ = L.rounded().lp(L.cellrow([0], 0, 0), 1)
report("N10b fractional LP endpoint (one unit, realised ordering)", "V", abs(v_lp - .2) < 1e-9 and abs(v_r - 1) < 1e-9,
       f"N=5, V=1, yes=1: LP lower veteran-yes {v_lp:.3f}; after rounding the ordering row (a box): {v_r:.0f}")

# N10c: three 2x2 units, published veteran-yes total 11, and two average-gap rows with
# different unit weights: sum_i w_i (y_i/V_i - (Y_i-y_i)/(N_i-V_i)) >= 0 with w = 1 (unweighted
# mean of unit gaps) and w_i = N_i - V_i (non-veteran weighted). These are NOT the realised
# associations of 2c M1 (those are rows on the aggregate count: S1, S2 below).
Ns, Vs, Ys = [13, 10, 14], [2, 8, 9], [13, 8, 4]
L = CountLP(a=[[Vs[i], Ns[i] - Vs[i]] for i in range(3)], b=[[Ys[i], Ns[i] - Ys[i]] for i in range(3)])
L.total(range(3), 0, 0, lo=11, hi=11)
for w in ([1, 1, 1], [Ns[i] - Vs[i] for i in range(3)]):
    row = np.zeros(L.nv); rhs = 0.0
    for i in range(3):
        row[L.ix(i, 0, 0)] = w[i] * (1 / Vs[i] + 1 / (Ns[i] - Vs[i])); rhs += w[i] * Ys[i] / (Ns[i] - Vs[i])
    L.ge(row, rhs)
t1 = L.cellrow([0], 0, 0)
v_lp, x_lp = L.lp(t1, 1); v_r, _ = L.rounded().lp(t1, 1); v_ip, _ = L.ip(t1, 1)
cnt_int = 0
for y in itertools.product(*[range(dd_count(Ns[i], Vs[i], Ys[i])[0], dd_count(Ns[i], Vs[i], Ys[i])[1] + 1) for i in range(3)]):
    if sum(y) != 11: continue
    if all(sum(w[i] * (Fr(y[i], Vs[i]) - Fr(Ys[i] - y[i], Ns[i] - Vs[i])) for i in range(3)) >= 0
           for w in ([1, 1, 1], [Ns[i] - Vs[i] for i in range(3)])): cnt_int += 1
yv = [x_lp[L.ix(i, 0, 0)] for i in range(3)]
report("N10c LP-feasible, no integer table (coupled rows, even after rounding each row)", "V",
       v_lp is not None and v_r is not None and v_ip is None and cnt_int == 0,
       f"N={Ns}, V={Vs}, yes={Ys}, total 11, average-gap rows with w=1 and w=N-V: LP witness veteran-yes "
       f"({yv[0]:.3f}, {yv[1]:.3f}, {yv[2]:.3f}); LP over Chvatal-rounded rows feasible: {v_r is not None}; "
       f"MILP: infeasible; enumeration: {cnt_int} integer tables")

# ---------------------------------------------------------------------------
# C1: rows outside the TU classes can break integrality (R19(c)): four row types
# ---------------------------------------------------------------------------
def search_gap(build, n=1500):
    for _ in range(n):
        L, tabs = build()
        Lr = L.rounded(); tgt = rng.integers(-2, 3, size=L.nv)
        for sense in (1, -1):
            v, x = Lr.lp(tgt, sense)
            if v is None: continue
            vi, _ = Lr.ip(tgt, sense)
            if vi is None: return ("LP feasible, integer infeasible", v, None, tabs)
            rnd = math.ceil(v - 1e-9) if sense == 1 else math.floor(v + 1e-9)
            if abs(vi - rnd) > 1e-7: return ("integer endpoint beyond the rounded LP endpoint", v, vi, tabs)
            if abs(vi - v) > 1e-7: return ("fractional LP endpoint", v, vi, tabs)
    return None


def build_pair_order():
    I = rng.integers(1, 3); tabs = [rng.integers(0, 7, size=(3, rng.integers(2, 4))) for _ in range(I)]
    L = CountLP(tabs)
    for i in range(I):
        row = np.zeros(L.nv); row[L.ix(i, 0, 0)] = L.a[i][1]; row[L.ix(i, 1, 0)] = -L.a[i][0]; L.ge(row, 0)
    return L, tabs


def build_pair_dB():
    I = rng.integers(2, 4); tabs = [rng.integers(1, 6, size=(3, 2)) for _ in range(I)]
    L = CountLP(tabs); row = np.zeros(L.nv)
    for i in range(I):
        Ni = L.a[i].sum(); row[L.ix(i, 0, 0)] += L.a[i][1] / Ni; row[L.ix(i, 1, 0)] -= L.a[i][0] / Ni
    L.ge(row, 0.0); return L, tabs


def build_two_cells():
    I = rng.integers(2, 4); tabs = [rng.integers(0, 5, size=(3, 3)) for _ in range(I)]
    L = CountLP(tabs); tv = sum(t[0, 0] for t in tabs); L.total(range(I), 0, 0, lo=tv, hi=tv)
    tv2 = tabs[0][1, 1] + tabs[1][1, 1]; L.total([0, 1], 1, 1, lo=tv2, hi=tv2); return L, tabs


g_ord = search_gap(build_pair_order); g_dB = search_gap(build_pair_dB); g_two = search_gap(build_two_cells, n=4000)
# crossing totals over three pairwise-crossing unit sets (e.g. totals published by three geographies)
L = CountLP(a=[[1, 1]] * 3, b=[[1, 1]] * 3)
for S in ([0, 1], [1, 2], [0, 2]): L.total(S, 0, 0, lo=1, hi=1)
v_x, x_x = L.lp(L.cellrow([0], 0, 0), 1); vi_x, _ = L.ip(L.cellrow([0], 0, 0), 1)
report("C1 rows outside the TU classes: integrality fails for each of four row types", "V",
       all(g is not None for g in (g_ord, g_dB, g_two)) and v_x is not None and vi_x is None,
       f"pairwise unit ordering (R=3): {g_ord[0]} (LP {g_ord[1]:.3f}, integer {g_ord[2]}); "
       f"pairwise realised delta_B(r,s) (R=3): {g_dB[0]} (LP {g_dB[1]:.3f}, integer {g_dB[2]}); "
       f"totals on cells (1,1) and (2,2) over overlapping units (R=3): {g_two[0]} (LP {g_two[1]:.3f}, integer {g_two[2]}); "
       f"three crossing totals (units N=2, V=1, yes=1; each pair total 1): LP witness y = {x_x[L.ix(0,0,0)]:.1f} each, integer infeasible")

print(f"(witness check so far: {WITNESS_LOG['n']} LP/MILP solutions, {WITNESS_LOG['bad']} violating a margin, row or bound)")

print("=== Section 3: realised signs and CR in count coordinates ===")


# ---------------------------------------------------------------------------
# S1: 2c Def. 1 associations as functions of the counts
# ---------------------------------------------------------------------------
def assoc_from_rates(tabs, r, c):
    """2c Definition 1 with p_i = N_i/N, X = a/N_i, beta = n/a (shares and rates)"""
    N_i = np.array([t.sum() for t in tabs], float); N = N_i.sum(); p = N_i / N
    X = np.array([t.sum(1) / t.sum() for t in tabs]); beta = np.array([t[:, c] / t.sum(1) for t in tabs])
    T = np.array([t[:, c].sum() / t.sum() for t in tabs])
    R = X.shape[1]; mu = p @ X
    dBpair = {s: float(np.sum(p * X[:, r] * X[:, s] * (beta[:, r] - beta[:, s]))) for s in range(R) if s != r}
    dBr = sum(dBpair.values())
    wr = p * X[:, r] / mu[r]; wn = p * (1 - X[:, r]) / (1 - mu[r])
    bn = (T - X[:, r] * beta[:, r]) / (1 - X[:, r])
    cov = lambda w, f, g: float(np.sum(w * f * g) - np.sum(w * f) * np.sum(w * g))
    dWr = mu[r] * cov(wr, beta[:, r], X[:, r]) + (1 - mu[r]) * cov(wn, bn, X[:, r])
    return dBpair, dBr, dWr


maxd = 0.0
for _ in range(200):
    I, R, C = rng.integers(2, 7), rng.integers(2, 5), rng.integers(2, 4)
    tabs = [rng.integers(1, 15, size=(R, C)) for _ in range(I)]
    r, c = rng.integers(R), rng.integers(C)
    dBpair, dBr, dWr = assoc_from_rates(tabs, r, c)
    a = np.array([t.sum(1) for t in tabs], float); bc = np.array([t[:, c].sum() for t in tabs], float)
    N_i = a.sum(1); N = N_i.sum(); Ar = a[:, r].sum(); Yc = bc.sum()
    Y = sum(t[r, c] for t in tabs)
    # count forms
    dBr_cnt = (Y - np.sum(a[:, r] * bc / N_i)) / N
    pair_cnt = {s: float(np.sum((a[:, s] * np.array([t[r, c] for t in tabs]) - a[:, r] * np.array([t[s, c] for t in tabs])) / N_i) / N)
                for s in range(R) if s != r}
    p = N_i / N; Xr = a[:, r] / N_i; Tc = bc / N_i
    varX = np.sum(p * Xr ** 2) - np.sum(p * Xr) ** 2
    ER = (np.sum(p * Xr * Tc) - np.sum(p * Xr) * np.sum(p * Tc)) / varX
    D = Y / Ar - (Yc - Y) / (N - Ar)
    dWr_cnt = varX * (ER - D)
    maxd = max(maxd, abs(dBr - dBr_cnt), abs(dWr - dWr_cnt), max(abs(dBpair[s] - pair_cnt[s]) for s in pair_cnt))
report("S1 realised associations in counts (one-vs-rest rows depend on n only through Y_rc)", "V", maxd < 1e-12,
       f"200 random R x C populations: delta_B^r = (Y_rc - sum_i a_ir b_ic/N_i)/N, "
       f"delta_W^r = Var_p(X_r)(ER^r - D^r(Y_rc)), pairwise delta_B(r,s) = sum_i (a_is n_irc - a_ir n_isc)/(N N_i)/N: "
       f"max |diff| = {maxd:.1e}")


# helpers: thresholds Y^NM and Y^ER(eps) for the one-vs-rest rows on cell (r,c)
def thresholds(a, b, r, c, eps_B=0.0, eps_W=0.0):
    a = np.asarray(a, float); b = np.asarray(b, float)
    N_i = a.sum(1); N = N_i.sum(); Ar = a[:, r].sum(); Yc = b[:, c].sum()
    p = N_i / N; Xr = a[:, r] / N_i; Tc = b[:, c] / N_i
    varX = np.sum(p * Xr ** 2) - np.sum(p * Xr) ** 2
    ER = (np.sum(p * Xr * Tc) - np.sum(p * Xr) * np.sum(p * Tc)) / varX
    YNM = np.sum(a[:, r] * b[:, c] / N_i) - N * eps_B          # delta_B^r >= -eps_B  <=>  Y >= YNM
    kap = 1 / Ar + 1 / (N - Ar)
    YER = (ER + eps_W / varX + Yc / (N - Ar)) / kap             # delta_W^r >= -eps_W  <=>  Y <= YER
    return YNM, YER, kap, Ar, N, Yc


# ---------------------------------------------------------------------------
# S2 [F]: one-vs-rest signs and CR are integer-exact by the rounded LP (no MILP needed)
# ---------------------------------------------------------------------------
nbad = 0; nfr = 0; nch = 0; nbr = 0
for _ in range(250):
    I, R, C = rng.integers(2, 6), rng.integers(2, 4), rng.integers(2, 4)
    tabs = [rng.integers(0, 9, size=(R, C)) for _ in range(I)]
    a = [t.sum(1) for t in tabs]; b = [t.sum(0) for t in tabs]
    if any(x.sum() == 0 for x in a): continue
    r, c = rng.integers(R), rng.integers(C)
    try:
        YNM, YER, *_ = thresholds(a, b, r, c)
    except ZeroDivisionError:
        continue
    if not np.isfinite(YER): continue
    for branch in (+1, -1):                       # CR = (dB>=0, dW>=0) U (dB<=0, dW<=0)
        L = CountLP(tabs)
        if branch == 1: L.total(range(I), r, c, lo=YNM, hi=YER)
        else: L.total(range(I), r, c, lo=YER, hi=YNM)
        i, rr, cc = rng.integers(I), rng.integers(R), rng.integers(C)       # a box and a known cell
        L.box(i, rr, cc, lo=max(0, tabs[i][rr, cc] - 1.3))
        j = rng.integers(I); L.box(j, R - 1, 0, lo=tabs[j][R - 1, 0], hi=tabs[j][R - 1, 0])
        Lr = L.rounded(); tgt = rng.integers(-3, 4, size=L.nv)
        for sense in (1, -1):
            v, x = Lr.lp(tgt, sense); vi, _ = L.ip(tgt, sense); nch += 1
            if (v is None) != (vi is None) or (v is not None and abs(v - vi) > 1e-7): nbad += 1
            nfr += frac(x) > 1e-7
        nbr += 1
report("S2 one-vs-rest realised signs and CR: rounded LP = MILP (TU class (i))", "F", nbad == 0 and nfr == 0,
       f"{nbr} CR branches (R x C up to 3 x 3, with a box and a known cell), {nch} endpoints: {nbad} mismatches, "
       f"{nfr} fractional vertices (fails on 'pooled realised signs and CR break TU, so realised CR is two MILPs "
       "whenever integrality matters'; that holds for pairwise/W2 rows, C1)")


# ---------------------------------------------------------------------------
# S3: 2c M3 in count coordinates; integer sets; the referee examples of 2c Remark M2.1
# ---------------------------------------------------------------------------
def gap_sets(Ns, Vs, Ys, r=0, c=0):
    a = [[Vs[i], Ns[i] - Vs[i]] for i in range(len(Ns))]; b = [[Ys[i], Ns[i] - Ys[i]] for i in range(len(Ns))]
    YNM, YER, kap, Ar, N, Yc = thresholds(a, b, 0, 0)
    Ylo = sum(dd_count(Ns[i], Vs[i], Ys[i])[0] for i in range(len(Ns))); Yhi = sum(dd_count(Ns[i], Vs[i], Ys[i])[1] for i in range(len(Ns)))
    D = lambda Y: Y * kap - Yc / (N - Ar)
    return dict(YNM=YNM, YER=YER, Ylo=Ylo, Yhi=Yhi, D=D, kap=kap)


# referee example 1: three units of size 4, group-1 counts 1,2,3, positives 0,1,0
g = gap_sets([4, 4, 4], [1, 2, 3], [0, 1, 0])
cr_cont = (max(g["Ylo"], min(g["YNM"], g["YER"])), min(g["Yhi"], max(g["YNM"], g["YER"])))
cr_int = [Y for Y in range(g["Ylo"], g["Yhi"] + 1) if min(g["YNM"], g["YER"]) - 1e-9 <= Y <= max(g["YNM"], g["YER"]) + 1e-9]
# referee example 2: two four-voter precincts with 1 and 2 A voters, one P vote each
g2 = gap_sets([4, 4], [1, 2], [1, 1])
Bp_int_lo = math.ceil(g2["YNM"] - 1e-9); gaps2 = sorted(set(g2["D"](Y) for Y in range(g2["Ylo"], g2["Yhi"] + 1)))
report("S3 2c M3 and Remark M2.1 in counts: thresholds Y^NM, Y^ER; integer sets = lattice", "V",
       abs(g["D"](cr_cont[0])) < 1e-12 and abs(g["D"](cr_cont[1])) < 1e-12 and cr_int == []
       and abs(g2["D"](g2["YNM"])) < 1e-12 and abs(g2["D"](Bp_int_lo) - 2 / 15) < 1e-12
       and np.allclose(gaps2, [-2 / 5, 2 / 15, 2 / 3]),
       f"example 1: Y^NM = {g['YNM']:.3f}, Y^ER = {g['YER']:.3f}, continuous CR gap image "
       f"{{{g['D'](cr_cont[0]):.3f}}}, integer CR set empty; example 2: integer gaps "
       f"{[round(float(x),4) for x in gaps2]}, between sign: continuous lower gap {g2['D'](g2['YNM']):.3f}, "
       f"integer (Y >= ceil(Y^NM) = {Bp_int_lo}) {g2['D'](Bp_int_lo):.4f} = 2/15")

# S3b: closed forms against the LP on random 2x2 populations (continuous), and the lattice
maxd = 0.0; nlat = 0
for _ in range(200):
    I = rng.integers(2, 6); Ns = rng.integers(4, 30, size=I).tolist()
    Vs = [int(rng.integers(1, n)) for n in Ns]; Ys = [int(rng.integers(0, n + 1)) for n in Ns]
    g = gap_sets(Ns, Vs, Ys)
    if not np.isfinite(g["YER"]): continue
    L = CountLP(a=[[Vs[i], Ns[i] - Vs[i]] for i in range(I)], b=[[Ys[i], Ns[i] - Ys[i]] for i in range(I)])
    L.total(range(I), 0, 0, lo=g["YNM"])          # delta_B >= 0
    tgt = L.cellrow(range(I), 0, 0)
    lo, _ = L.lp(tgt, 1); hi, _ = L.lp(tgt, -1)
    maxd = max(maxd, abs(lo - max(g["Ylo"], g["YNM"])), abs(hi - g["Yhi"]))
    lo_i, _ = L.rounded().lp(tgt, 1); nlat += abs(lo_i - max(g["Ylo"], math.ceil(g["YNM"] - 1e-9))) > 1e-9
report("S3b between sign in counts: LP = [max(Y^-, Y^NM), Y^+]; integer = ceil", "V", maxd < 1e-8 and nlat == 0,
       f"200 random 2x2 populations: max |LP - closed form| = {maxd:.1e}; integer lower end = "
       f"max(Y^-, ceil Y^NM) in all ({nlat} mismatches); Y^NM in [Y^-, Y^+] always (never refuted without unit rows)")

# S4: a known cell can remove the equal-rates witness, so the between sign becomes refutable
Ns, Vs, Ys = [10, 10], [5, 5], [5, 5]       # equal-rates witness y = (2.5, 2.5); Y^NM = 5
g = gap_sets(Ns, Vs, Ys)
L = CountLP(a=[[5, 5], [5, 5]], b=[[5, 5], [5, 5]]); L.total([0, 1], 0, 0, lo=g["YNM"])
L.box(0, 0, 0, lo=0, hi=0); L.box(1, 0, 0, hi=4)     # unit 1: no veteran said yes; unit 2: at most 4
feas, _ = L.lp(L.cellrow([0, 1], 0, 0), 1)
report("S4 unit rows can make a between sign refutable (2c M3(c) needs no unit rows)", "V", feas is None,
       f"two units N=10, V=5, yes=5; Y^NM = {g['YNM']:.1f}; known cell y_1 = 0 and box y_2 <= 4: "
       f"between-sign system infeasible = {feas is None}")

print("=== Worked examples (count layer) ===")
# E1: realised ordering, two precincts (P1: 100, 20 V, 30 yes; P2: 100, 60 V, 50 yes)
L = CountLP(a=[[20, 80], [60, 40]], b=[[30, 70], [50, 50]])
tgt = L.cellrow([0, 1], 0, 0) / 80
dd_lo, _ = L.lp(tgt, 1); dd_hi, _ = L.lp(tgt, -1)
for i, (V, Y, N) in enumerate(((20, 30, 100), (60, 50, 100))): L.box(i, 0, 0, lo=V * Y / N)   # one-vs-rest unit ordering = box
lo, wlo = L.rounded().lp(tgt, 1); hi, whi = L.rounded().lp(tgt, -1)
ok_w = all(L.feasible_point(w) for w in (wlo, whi)) and frac(wlo) < 1e-9 and frac(whi) < 1e-9
report("E1 realised ordering example", "V", (dd_lo, dd_hi) == (0.125, 0.875) and abs(lo - .45) < 1e-12 and abs(hi - .875) < 1e-12 and ok_w,
       f"DD [{dd_lo:.3f}, {dd_hi:.3f}]; ordering rows y_1 >= 6, y_2 >= 30: [{lo:.3f}, {hi:.3f}]; witnesses: "
       f"lower veteran-yes ({wlo[L.ix(0,0,0)]:.0f}, {wlo[L.ix(1,0,0)]:.0f}), upper ({whi[L.ix(0,0,0)]:.0f}, {whi[L.ix(1,0,0)]:.0f}), integer, every row holds")

# E2: published total for P1+P2 veteran-yes = 45 (E5)
L = CountLP(a=[[20, 80], [60, 40]], b=[[30, 70], [50, 50]]); L.total([0, 1], 0, 0, lo=45, hi=45)
lo2, _ = L.lp(L.cellrow([1], 0, 0), 1); hi2, _ = L.lp(L.cellrow([1], 0, 0), -1)
loB, _ = L.lp(L.cellrow([0, 1], 0, 0) / 80, 1); hiB, _ = L.lp(L.cellrow([0, 1], 0, 0) / 80, -1)
report("E2 published total example", "V", (round(lo2), round(hi2)) == (25, 45) and abs(loB - hiB) < 1e-12 and abs(loB - .5625) < 1e-12,
       f"total 45: P2 veteran-yes in [{lo2:.0f}, {hi2:.0f}] (rate [{lo2/60:.3f}, {hi2/60:.3f}]) vs DD [10, 50]; "
       f"the two-precinct veteran rate is pinned at {loB:.4f}")


# E3: odds-ratio dial on a 2x2 unit is a box on the cell (P2: N=100, V=60, yes=50)
def or_cell(N, V, Y, psi):
    """n11 with odds ratio n11*n00/(n10*n01) = psi, given the margins (monotone in n11)"""
    lo, hi = max(0, V + Y - N), min(V, Y)
    f = lambda y: y * (N - V - Y + y) - psi * (V - y) * (Y - y)
    a_, b_ = lo + 1e-12, hi - 1e-12
    for _ in range(200):
        m = (a_ + b_) / 2
        if f(m) < 0: a_ = m
        else: b_ = m
    return (a_ + b_) / 2


L = CountLP(a=[[60, 40]], b=[[50, 50]])
L.box(0, 0, 0, lo=or_cell(100, 60, 50, 1.0), hi=or_cell(100, 60, 50, 4.0))
lo3, _ = L.lp(L.cellrow([0], 0, 0) / 60, 1); hi3, _ = L.lp(L.cellrow([0], 0, 0) / 60, -1)
report("E3 odds-ratio dial [1, 4] = a box on the cell", "V", abs(lo3 - .5) < 1e-6 and abs(hi3 - .6333) < 1e-3,
       f"P2 veteran rate: DD [0.167, 0.833]; OR in [1, 4]: [{lo3:.3f}, {hi3:.3f}] (N17)")

# E4: slack dial on the between sign: nested, piecewise-linear endpoints (2c M2(d))
Ns, Vs, Ys = [100, 100], [20, 60], [30, 50]
a = [[20, 80], [60, 40]]; b = [[30, 70], [50, 50]]
prev = None; nested = True; los = []
for eps in (0, .005, .01, .02, .04):
    YNM, *_ = thresholds(a, b, 0, 0, eps_B=eps)
    L = CountLP(a=a, b=b); L.total([0, 1], 0, 0, lo=YNM)
    lo_e, _ = L.lp(L.cellrow([0, 1], 0, 0) / 80, 1); los.append(lo_e)
    nested &= prev is None or lo_e <= prev + 1e-12; prev = lo_e
report("E4 slack dial on delta_B >= -eps: lower end nonincreasing in eps", "V", nested,
       "veteran rate lower end for eps = 0, .005, .01, .02, .04: " + ", ".join(f"{x:.4f}" for x in los))

# non-integer margins: exact for the stated margins (R19(e)); the N20 first number
L = CountLP(a=[[60.0, 40.0]], b=[[50.0, 50.0]])
lo_, _ = L.lp(L.cellrow([0], 0, 0) / 60, 1); hi_, _ = L.lp(L.cellrow([0], 0, 0) / 60, -1)
L2 = CountLP(a=[[59.4, 40.6]], b=[[50.0, 50.0]])
lo2_, w2 = L2.lp(L2.cellrow([0], 0, 0) / 59.4, 1)
report("E5 non-integer (estimated) margins: the LP is exact for the stated margins", "V",
       abs(lo_ - 1 / 6) < 1e-9 and abs(hi_ - 5 / 6) < 1e-9 and L2.feasible_point(w2),
       f"expected count 60 treated as exact: [{lo_:.3f}, {hi_:.3f}] (N20); margin 59.4: lower end {lo2_:.4f}, "
       "attained by a real table with exactly those margins; no integer table has them, so the integer check is off")

print(f"(witness check: {WITNESS_LOG['n']} LP/MILP solutions in Sections 2-3, {WITNESS_LOG['bad']} violating a margin, row or bound)")
report("W1 every witness satisfies every margin, row and bound", "V", WITNESS_LOG["bad"] == 0,
       f"{WITNESS_LOG['n']} LP/MILP solutions checked: {WITNESS_LOG['bad']} violations")

print("=== Section 4: R27, the cross-layer combination ===")
C_cnt = (0.45, 0.875)
inter = lambda A, B: (max(A[0], B[0]), min(A[1], B[1])) if max(A[0], B[0]) <= min(A[1], B[1]) else None
report("X1 cross-layer worked example", "I", inter(C_cnt, (.30, .60)) == (.45, .60) and inter(C_cnt, (.20, .40)) is None,
       f"count set [0.45, 0.875] with model interval [0.30, 0.60]: {inter(C_cnt, (.30, .60))}; "
       f"with [0.20, 0.40]: disjoint, a level-alpha rejection of (model and its rows) given the count rows")


# ---------------------------------------------------------------------------
# A small 2x2 pooled-M0 (LCE) model layer with the projection interval P3.5(b2)(ii) and
# the Popoviciu allowance of P3.5(a). m_0(x) = q0 + k x, m_1(x) = q0 + q1 + q2 x + k (x - 1),
# tau(x) = q0 + q1 x + q2 x^2; k (the slope of m_0) is the free coordinate.
# ---------------------------------------------------------------------------
def model_interval(X, T, w, alpha_est=0.04, alpha_E=0.01, allowance=True, plug_in=False):
    Z = np.column_stack([np.ones_like(X), X, X ** 2])
    ZtZi = np.linalg.inv(Z.T @ Z); q = ZtZi @ Z.T @ T; e = T - Z @ q
    n = len(X); V = ZtZi @ (Z.T * e ** 2) @ Z @ ZtZi * n / (n - 3)          # HC1
    xbar = w @ X; gam = xbar - 1.0
    vs = (X.min(), X.max())
    # candidate endpoints theta_j(q) = g_j . q + h_j
    up, lo = [], []
    cgrad = np.array([1.0, 1.0, xbar])
    for v in vs:
        # lower bounds on k (give upper ends of B~ since gam < 0)
        up.append((cgrad + gam * np.array([-1 / v, 0, 0]), 0.0))                           # m_0(v) >= 0
        up.append((cgrad + gam * np.array([1, 1, v]) / (1 - v), -gam / (1 - v)))             # m_1(v) <= 1
        # upper bounds on k (give lower ends)
        lo.append((cgrad + gam * np.array([-1 / v, 0, 0]), gam / v))                         # m_0(v) <= 1
        lo.append((cgrad + gam * np.array([1, 1, v]) / (1 - v), 0.0))                        # m_1(v) >= 0
    kap = 0.0 if plug_in else math.sqrt(chi2.ppf(1 - alpha_est, 3))
    U = min(g @ q + h + kap * math.sqrt(g @ V @ g) for g, h in up)
    Lo = max(g @ q + h - kap * math.sqrt(g @ V @ g) for g, h in lo)
    S1 = 0.5 * math.sqrt(np.sum(w ** 2)); zE = norm.ppf(1 - alpha_E / 2)
    a_ = zE * S1 if (allowance and not plug_in) else 0.0
    return Lo - a_, U + a_


def simulate(rep, known_frac=0.9, false_row=False, mode="ci", p=80, het=0.0, design="interior", seed=7):
    r_ = np.random.default_rng(seed)
    out = dict(cov_ci=0, cov_int=0, rej=0, n=0, w_int=[], w_joint=[], cov_joint=0, narrower=0)
    for _ in range(rep):
        X = r_.uniform(0.05, 0.95, size=p); N = r_.integers(50, 151, size=p)
        V = np.clip(np.rint(N * X).astype(int), 1, N - 1); X = V / N
        m1 = 0.30 + 0.40 * X
        m0 = 0.20 + 0.20 * X if design == "interior" else 0.30 * (X - 0.05)   # boundary: m_0 near 0 at the hull vertex
        if het > 0:            # unit heterogeneity in group 1 around its mean (E[beta | V] = m still holds)
            m1 = m1 + het * r_.uniform(-1, 1, size=p)
        y1 = r_.binomial(V, m1); y0 = r_.binomial(N - V, np.clip(m0, 0, 1))
        T = (y1 + y0) / N; w = V / V.sum(); B = y1.sum() / V.sum()
        ci = model_interval(X, T, w, allowance=(mode != "no_allowance"), plug_in=(mode == "plug_in"))
        # count layer: margins, plus known veteran-yes cells (true values) in a fraction of units
        known = r_.random(p) < known_frac
        dlo = np.maximum(0, V + y1 + y0 - N); dhi = np.minimum(V, y1 + y0)
        lo_i = np.where(known, y1, dlo); hi_i = np.where(known, y1, dhi)
        if false_row:          # a false realised-table belief: 'every known unit's veteran-yes count is at its DD minimum'
            lo_i = np.where(known, dlo, lo_i); hi_i = np.where(known, dlo, hi_i)
        Cset = (lo_i.sum() / V.sum(), hi_i.sum() / V.sum())   # the count LP separates (one cell per unit)
        I_ = inter(Cset, ci)
        out["n"] += 1; out["cov_ci"] += ci[0] <= B <= ci[1]
        out["cov_int"] += I_ is not None and I_[0] - 1e-12 <= B <= I_[1] + 1e-12
        out["rej"] += I_ is None
        if mode == "ci" and not false_row and I_ is not None and (~known).sum() >= 1:
            # joint construction (Remark 4.4): the known part exactly, the model interval only for the rest
            U_ = ~known; wU = V * U_ / (V * U_).sum()
            ciU = model_interval_subset(X, T, wU)
            AU = (V * U_).sum(); Akn = y1[known].sum()
            J = inter(((Akn + AU * ciU[0]) / V.sum(), (Akn + AU * ciU[1]) / V.sum()), Cset)
            out["cov_joint"] += J is not None and J[0] - 1e-12 <= B <= J[1] + 1e-12
            wj = 0 if J is None else J[1] - J[0]; out["w_joint"].append(wj); out["w_int"].append(I_[1] - I_[0])
            out["narrower"] += wj < (I_[1] - I_[0]) - 1e-12
    return out


def model_interval_subset(X, T, wU):
    return model_interval(X, T, wU)


REP = 1000
s_dd = simulate(REP, known_frac=0.0, seed=11)
s_kn = simulate(REP, known_frac=0.9, seed=12)
s_big = simulate(REP, known_frac=0.9, p=1000, seed=18)
mc = lambda k, n: math.sqrt(max(k, 1) / n * (1 - k / n) / n)
report("X2 size of the disjointness test and coverage of the intersection, model holds", "V",
       max(s_kn["rej"], s_big["rej"]) / REP <= 0.05 + 2 * mc(max(s_kn["rej"], s_big["rej"]), REP)
       and all(s_["cov_int"] == s_["cov_ci"] for s_ in (s_dd, s_kn, s_big)),
       f"pooled 2x2 LCE, 80 units, projection interval + Popoviciu allowance at alpha = .05, {REP} reps: "
       f"margins only: coverage of CI {s_dd['cov_ci']/REP:.3f}, of the intersection {s_dd['cov_int']/REP:.3f}, rejection rate "
       f"{s_dd['rej']/REP:.3f}; known veteran-yes cells in 90% of units: coverage {s_kn['cov_ci']/REP:.3f} / "
       f"{s_kn['cov_int']/REP:.3f}, rejection rate {s_kn['rej']/REP:.3f}; 1000 units, 90% known: coverage "
       f"{s_big['cov_ci']/REP:.3f} / {s_big['cov_int']/REP:.3f}, rejection rate {s_big['rej']/REP:.3f} (MC s.e. {mc(s_big['rej'],REP):.3f})")

s_pi = simulate(REP, known_frac=1.0, mode="plug_in", het=0.15, design="boundary", seed=13)
s_na = simulate(REP, known_frac=1.0, mode="no_allowance", het=0.15, design="boundary", seed=14)
s_ok = simulate(REP, known_frac=1.0, mode="ci", het=0.15, design="boundary", seed=15)
report("X3 the model interval must cover the actual B (P3.5 with its allowance)", "F",
       s_pi["rej"] / REP > 0.05 + 3 * mc(s_pi["rej"], REP) and s_ok["rej"] / REP <= 0.05 + 2 * mc(s_ok["rej"], REP),
       f"count layer pins B (all cells known), truth on the boundary of the LCE set (m_0 = 0 at x = .05), group-1 unit heterogeneity +/- .15: rejection rate with the "
       f"plug-in model set {s_pi['rej']/REP:.3f}; with the interval for B~ (no residual allowance) {s_na['rej']/REP:.3f}; "
       f"with the P3.5 interval for B {s_ok['rej']/REP:.3f} (fails on 'any model-layer set may be intersected')")

s_fr = simulate(REP, known_frac=0.9, false_row=True, p=1000, seed=16)
report("X4 disjointness tests the conjunction: a false count row is rejected too", "F", s_fr["rej"] / REP > 0.2,
       f"model true; false count rows (the known cells of 90% of units set to their DD minima), 1000 units: rejection rate "
       f"{s_fr['rej']/REP:.3f} (fails on 'a disjoint answer rejects the model, not the realised-table rows')")

s_j = simulate(REP, known_frac=0.5, p=1000, seed=17)
wj = np.array(s_j["w_joint"]); wi = np.array(s_j["w_int"])
report("X5 the intersection is not the smallest valid set (R27(iii))", "I",
       s_j["narrower"] > 0 and s_j["cov_joint"] / max(1, len(wj)) >= 0.95 - 2 * mc(len(wj) - s_j["cov_joint"], max(1, len(wj))),
       f"1000 units, known cells in 50% of units: joint set (known part exact + model interval for the rest, own allowance) "
       f"narrower than C ∩ CI in {s_j['narrower']}/{len(wj)} reps; median widths {np.median(wj):.4f} vs {np.median(wi):.4f}; "
       f"coverage: joint {s_j['cov_joint']/max(1,len(wj)):.3f}, C ∩ CI {s_j['cov_int']/REP:.3f}")

# G1: the residual allowance for a general linear count target (Lemma 4.1): Var(E_phi) <= S_phi^2.
# Gap target phi = Y_1/A_1 - Y_0/A_0; within each unit the two group rates have an arbitrary
# (dependent) joint law on [0,1]^2; exact variances from random discrete laws.
maxratio = 0.0; attained = 0.0
for _ in range(200):
    p = 15; N = rng.integers(10, 60, size=p); V = np.minimum(rng.integers(1, 9, size=p), N - 1)
    c1 = V / V.sum(); c0 = -(N - V) / (N - V).sum()          # coefficients on the unit rates beta_i1, beta_i0
    rho = np.abs(c1) + np.abs(c0); S2 = 0.25 * np.sum(rho ** 2)
    var = 0.0
    for i in range(p):
        k = rng.integers(2, 5); pts = rng.random((k, 2)); pr = rng.dirichlet(np.ones(k))
        s_ = c1[i] * pts[:, 0] + c0[i] * pts[:, 1]; var += pr @ s_ ** 2 - (pr @ s_) ** 2
    maxratio = max(maxratio, var / S2)
    # extremal law: (beta_1, beta_0) = (1, 0) or (0, 1) with probability 1/2 each
    var_ex = np.sum(0.25 * (c1 - c0) ** 2); attained = max(attained, var_ex / S2)
report("G1 general-target allowance: Var(E_phi) <= S_phi^2 under within-unit dependence; the bound is attained", "V",
       maxratio <= 1 + 1e-12 and abs(attained - 1) < 1e-12,
       f"200 random designs (gap target, random dependent joint laws of the two group rates per unit): max Var/S^2 = "
       f"{maxratio:.4f}; perfectly opposite rates attain {attained:.4f}")

print()
print(f"{sum(results)}/{len(results)} PASS")
