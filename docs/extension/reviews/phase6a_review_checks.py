"""Independent checks for the proof review of PHASE6a_count_layer.md (M-COUNT, step 1).

Written by the reviewer, not the author. Each check is labelled RC<k> and is cited in
PHASE6a_proof_review.md. Tags:
  [CE]  counterexample to a statement as written in 6a
  [V]   independent verification of a 6a claim (different generator / method from verify_phase6a.py)
  [B]   boundary: shows a hypothesis of 6a is needed (not an error in 6a)
Run: python3 phase6a_review_checks.py   (numpy, scipy >= 1.9)
"""
import itertools, math
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

rng = np.random.default_rng(20261121)
RESULTS = []


def report(label, tag, ok, msg):
    RESULTS.append(ok)
    print(f"{'PASS' if ok else 'FAIL'} {label} [{tag}]: {msg}")


# ---------------------------------------------------------------------------
# a small count-layer builder (independent of verify_phase6a.py)
# ---------------------------------------------------------------------------
class Sys:
    def __init__(self, a, b):
        self.a = [np.asarray(x, float) for x in a]; self.b = [np.asarray(x, float) for x in b]
        self.I, self.R, self.C = len(a), len(a[0]), len(b[0]); self.nv = self.I * self.R * self.C
        self.E, self.e, self.U, self.u = [], [], [], []
        self.lb = np.zeros(self.nv); self.ub = np.full(self.nv, np.inf)
        for i in range(self.I):
            for r in range(self.R): self.E.append(self.vec([(i, r, c) for c in range(self.C)])); self.e.append(self.a[i][r])
            for c in range(self.C): self.E.append(self.vec([(i, r, c) for r in range(self.R)])); self.e.append(self.b[i][c])

    def ix(self, i, r, c): return (i * self.R + r) * self.C + c

    def vec(self, cells, w=None):
        v = np.zeros(self.nv)
        for k, (i, r, c) in enumerate(cells): v[self.ix(i, r, c)] += 1 if w is None else w[k]
        return v

    def total(self, S, r, c, lo=-np.inf, hi=np.inf):
        v = self.vec([(i, r, c) for i in S])
        if lo > -np.inf: self.U.append(-v); self.u.append(-lo)
        if hi < np.inf: self.U.append(v); self.u.append(hi)

    def box(self, i, r, c, lo=0, hi=np.inf):
        k = self.ix(i, r, c); self.lb[k] = max(self.lb[k], lo); self.ub[k] = min(self.ub[k], hi)

    def rounded(self):
        o = Sys(self.a, self.b); o.E, o.e = list(self.E), list(self.e)
        o.lb = np.ceil(self.lb - 1e-9); o.ub = np.where(np.isinf(self.ub), np.inf, np.floor(self.ub + 1e-9))
        o.U = list(self.U); o.u = [math.floor(x + 1e-9) for x in self.u]   # all U rows here have 0/+-1 coefficients
        return o

    def _bnds(self): return [(l, None if np.isinf(h) else h) for l, h in zip(self.lb, self.ub)]

    def lp(self, cvec, sense=1, extra_eq=None):
        E, e = list(self.E), list(self.e)
        if extra_eq is not None: E.append(extra_eq[0]); e.append(extra_eq[1])
        kw = dict(A_eq=np.array(E), b_eq=np.array(e), bounds=self._bnds(), method="highs-ds")
        if self.U: kw.update(A_ub=np.array(self.U), b_ub=np.array(self.u))
        r = linprog(sense * np.asarray(cvec, float), **kw)
        return (None, None) if r.status != 0 else (sense * r.fun, r.x)

    def ip(self, cvec, sense=1, extra_eq=None):
        E, e = list(self.E), list(self.e)
        if extra_eq is not None: E.append(extra_eq[0]); e.append(extra_eq[1])
        cons = [LinearConstraint(np.array(E), np.array(e), np.array(e))]
        if self.U: cons.append(LinearConstraint(np.array(self.U), -np.inf, np.array(self.u)))
        m = milp(sense * np.asarray(cvec, float), constraints=cons, integrality=np.ones(self.nv),
                 bounds=Bounds(self.lb, self.ub))
        return (None, None) if m.status != 0 else (sense * m.fun, m.x)


def fracpart(x): return float(np.max(np.abs(x - np.round(x))))


def random_laminar(units, depth=3):
    """a random laminar family on `units` by recursive splitting (sets of size >= 2 only)"""
    units = list(units); fam = []
    def rec(S, d):
        if len(S) >= 2: fam.append(sorted(S))
        if d == 0 or len(S) < 2: return
        perm = list(rng.permutation(S)); k = rng.integers(1, len(S) + 1)
        cuts = sorted(rng.choice(range(1, len(S)), size=min(k - 1, len(S) - 1), replace=False)) if len(S) > 1 else []
        for blk in np.split(np.array(perm), cuts):
            if rng.random() < 0.8: rec(list(blk), d - 1)
    rec(units, depth)
    out = []
    for S in fam:
        if S not in out: out.append(S)
    return out


def is_laminar(fam):
    for A, B in itertools.combinations(fam, 2):
        A, B = set(A), set(B)
        if A & B and not (A <= B or B <= A): return False
    return True


# ===========================================================================
print("=== R19(b): total unimodularity ===")
# RC1 [CE]: class (iii) as stated ("the constraint matrix of the rounded system is TU") is false for the
# original (unsubstituted) matrix; the polytope is nevertheless integral (the proof's substitution).
I, R, C = 3, 2, 2
S0 = Sys(a=[[1, 1]] * I, b=[[1, 1]] * I)
fam_rows = [((0, 1, 2), 0, 0), ((0, 1), 1, 0), ((1, 2), 0, 1)]   # column 1: {012} on (1,1), {01} on (2,1); column 2: {12} on (1,2)
M = np.array(S0.E + [S0.vec([(i, r, c) for i in S]) for S, r, c in fam_rows])
assert is_laminar([S for S, r, c in fam_rows if c == 0]) and is_laminar([S for S, r, c in fam_rows if c == 1])
# odd cycle: rows T1, R(2,row1), T3, C(1,col2), R(1,row2), T2, C(0,col1); vars n0_11, n2_11, n2_12, n1_12, n1_22, n1_21, n0_21
ri = lambda i, r: i * (R + C) + r; ci = lambda i, c: i * (R + C) + R + c; T = len(S0.E)
rows = [T + 0, ri(2, 0), T + 2, ci(1, 1), ri(1, 1), T + 1, ci(0, 0)]
cols = [S0.ix(0, 0, 0), S0.ix(2, 0, 0), S0.ix(2, 0, 1), S0.ix(1, 0, 1), S0.ix(1, 1, 1), S0.ix(1, 1, 0), S0.ix(0, 1, 0)]
d = round(np.linalg.det(M[np.ix_(rows, cols)]))
# integrality nevertheless: random integer rhs / boxes, LP over rounded system = MILP, vertices integer
bad = frs = feas = 0
for _ in range(300):
    tabs = [rng.integers(0, 4, size=(R, C)) for _ in range(I)]
    L = Sys([t.sum(1) for t in tabs], [t.sum(0) for t in tabs])
    for S, r, c in fam_rows:
        tv = sum(tabs[i][r, c] for i in S); L.total(S, r, c, lo=tv - rng.integers(0, 3) - 0.4 * rng.random(), hi=tv + rng.integers(0, 3) + 0.4 * rng.random())
    Lr = L.rounded(); tg = rng.integers(-3, 4, size=L.nv)
    for s in (1, -1):
        v, x = Lr.lp(tg, s); vi, _ = Lr.ip(tg, s)
        if (v is None) != (vi is None) or (v is not None and abs(v - vi) > 1e-7): bad += 1
        if x is not None: feas += 1; frs += fracpart(x) > 1e-7
report("RC1 class (iii): the ORIGINAL constraint matrix is not TU, but F_rd is integral", "CE",
       abs(d) == 2 and bad == 0 and frs == 0,
       f"3 units 2x2; column 1 totals {{012}} on (1,1), {{01}} on (2,1); column 2 total {{12}} on (1,2) (one laminar family per "
       f"column): 7x7 minor of the margins+totals matrix = {d} (odd cycle), so 'the constraint matrix ... is TU' is false here; "
       f"300 random rhs: {bad} LP != MILP, {frs} fractional vertices out of {feas} (integrality holds via the substitution y = n_1c)")


# RC2 [V]: classes (i)-(iii) with deeper random laminar families, random rhs including INFEASIBLE systems
def trial(cls):
    I = int(rng.integers(3, 7))
    R, C = (2, int(rng.integers(2, 4))) if cls == "iii" else (int(rng.integers(2, 4)), int(rng.integers(2, 4)))
    tabs = [rng.integers(0, 5, size=(R, C)) for _ in range(I)]
    L = Sys([t.sum(1) for t in tabs], [t.sum(0) for t in tabs])
    def add(S, r, c):
        tv = sum(tabs[i][r, c] for i in S); shift = rng.integers(-2, 3)     # may be infeasible
        L.total(S, r, c, lo=tv + shift - rng.integers(0, 2) - 0.5 * rng.random(), hi=tv + shift + rng.integers(0, 2) + 0.5 * rng.random())
    if cls == "i":
        r0, c0 = int(rng.integers(R)), int(rng.integers(C))
        for S in random_laminar(range(I)) + random_laminar(range(I)): add(S, r0, c0)
    elif cls == "ii":
        c0 = int(rng.integers(C))
        for r in range(R):
            for S in random_laminar(range(I)): add(S, r, c0)
    else:
        for c in range(C):
            for S in random_laminar(range(I)): add(S, int(rng.integers(2)), c)
    for _ in range(3):
        i, r, c = int(rng.integers(I)), int(rng.integers(R)), int(rng.integers(C))
        lo = max(0, tabs[i][r, c] - 2 * rng.random()); L.box(i, r, c, lo=lo, hi=lo + 2.5 * rng.random())
    Lr = L.rounded(); tg = rng.integers(-3, 4, size=L.nv); out = [0, 0, 0, 0]   # feas-mismatch, value-mismatch, fractional, infeasible count
    for s in (1, -1):
        v, x = Lr.lp(tg, s); vi, _ = Lr.ip(tg, s)
        if (v is None) != (vi is None): out[0] += 1
        elif v is None: out[3] += 1
        else:
            out[1] += abs(v - vi) > 1e-7; out[2] += fracpart(x) > 1e-7
    return out


for cls in ("i", "ii", "iii"):
    tot = np.zeros(4, int)
    for _ in range(300): tot += trial(cls)
    report(f"RC2 class ({cls}) with deep laminar families and random (possibly infeasible) rhs", "V",
           tot[0] == 0 and tot[1] == 0 and tot[2] == 0,
           f"300 designs, 600 endpoints: feasibility LP(F_rd) vs MILP disagree {tot[0]}, values disagree {tot[1]}, "
           f"fractional vertices {tot[2]}; infeasible endpoints (both agree) {tot[3]}")

# RC3 [B]: R19(b) item 3 needs S addable: a crossing target set can have a gap in its integer values
L = Sys(a=[[1, 1]] * 3, b=[[1, 1]] * 3)
L.total([0, 1], 0, 0, lo=1, hi=1); L.total([1, 2], 0, 0, lo=1, hi=1)   # family A: {01}, family B: {12} (class (i))
tgt = L.vec([(0, 0, 0), (2, 0, 0)])
vals = [k for k in range(0, 3) if L.ip(tgt, 1, extra_eq=(tgt, k))[0] is not None]
lo_, _ = L.rounded().lp(tgt, 1); hi_, _ = L.rounded().lp(tgt, -1)
report("RC3 item 3 needs 'S addable': crossing target set in class (i) has integer values with a gap", "B",
       vals == [0, 2] and abs(lo_ - 0) < 1e-9 and abs(hi_ - 2) < 1e-9,
       f"3 units (N=2,V=1,yes=1), totals y0+y1=1 (family A), y1+y2=1 (family B); target y0+y2: rounded LP [{lo_:.0f},{hi_:.0f}] "
       f"(integer endpoints, as R19(b) item 2 says), integer values {vals}: 1 is not attained (consistent with 6a, which claims "
       f"every integer only for addable S)")


# RC4 [B]: the class boundaries are real: two families on one cell plus another cell of the column; R = 2
# with two crossing families on one column
def search_boundary(kind, n=4000):
    for _ in range(n):
        I = int(rng.integers(3, 6)); R, C = (3, 2) if kind == "ii2" else (2, 2)
        tabs = [rng.integers(0, 4, size=(R, C)) for _ in range(I)]
        L = Sys([t.sum(1) for t in tabs], [t.sum(0) for t in tabs])
        def add(S, r, c):
            tv = sum(tabs[i][r, c] for i in S); sh = rng.integers(-1, 2); L.total(S, r, c, lo=tv + sh, hi=tv + sh)
        if kind == "ii2":   # two crossing families on (0,0) and one more cell (1,0) of column 0
            for S in random_laminar(range(I), 1) + [sorted(rng.choice(I, size=2, replace=False).tolist())]: add(S, 0, 0)
            for S in [sorted(rng.choice(I, size=int(rng.integers(2, I + 1)), replace=False).tolist())]: add(S, 1, 0)
        else:               # R = 2, column 0 carries two crossing sets on its two cells
            S1 = sorted(rng.choice(I, size=2, replace=False).tolist()); S2 = sorted(rng.choice(I, size=2, replace=False).tolist())
            S3 = sorted(rng.choice(I, size=2, replace=False).tolist())
            add(S1, 0, 0); add(S2, 1, 0); add(S3, 0, 0)
        tg = rng.integers(-2, 3, size=L.nv)
        v, x = L.lp(tg, 1); vi, _ = L.ip(tg, 1)
        if v is not None and (vi is None or abs(v - vi) > 1e-7): return True, (I, v, vi)
    return False, None


f1, w1 = search_boundary("ii2"); f2, w2 = search_boundary("iii2")
report("RC4 class boundaries: (ii) with two families on one cell, and (iii) with crossing sets on a column, can lose integrality", "B",
       f1 and f2, f"(ii)+second family: {'LP %.3f vs integer %s (I=%d)' % (w1[1], w1[2], w1[0]) if f1 else 'none found'}; "
       f"(iii) with crossing sets: {'LP %.3f vs integer %s (I=%d)' % (w2[1], w2[2], w2[0]) if f2 else 'none found'}")

# RC10 [CE for the qualifier 'R = 3']: totals on cells in different rows and columns break integrality only
# when R >= 3 AND C >= 3; with C = 2 the transposed class (iii) applies (6a's own C1 example is 3 x 3).
bad = 0
for _ in range(1500):
    I = int(rng.integers(2, 4)); tabs = [rng.integers(0, 5, size=(3, 2)) for _ in range(I)]
    L = Sys([t.sum(1) for t in tabs], [t.sum(0) for t in tabs])
    tv = sum(t[0, 0] for t in tabs); L.total(range(I), 0, 0, lo=tv, hi=tv)
    tv2 = tabs[0][1, 1] + tabs[1][1, 1]; L.total([0, 1], 1, 1, lo=tv2, hi=tv2)
    tg = rng.integers(-2, 3, size=L.nv)
    for s_ in (1, -1):
        v, x = L.lp(tg, s_); vi, _ = L.ip(tg, s_)
        if v is not None and (vi is None or abs(v - vi) > 1e-7 or fracpart(x) > 1e-7): bad += 1
report("RC10 'totals on cells in different rows and columns' with R = 3, C = 2 stay integral (transposed class (iii))", "CE",
       bad == 0, f"1500 designs of 6a's C1 type with C = 2: {bad} fractional or LP != MILP endpoints; the breaking type needs R >= 3 and C >= 3")

# ===========================================================================
print("=== Prop. C2 ===")
# RC5 [CE]: C2(b) 'the integer gap set is the lattice ... inside the continuous interval' is false with
# fractional boxes on several units (T5-type); the formula with Y^-_rd is right, the sentence is not.
L = Sys(a=[[1, 4], [1, 4]], b=[[1, 4], [1, 4]])
for i in range(2): L.box(i, 0, 0, lo=0.2)                     # realised one-vs-rest ordering boxes a*b/N = 0.2
YNM = sum(1 * 1 / 5 for _ in range(2))                       # between sign over both units: Y >= Y^NM = 0.4
L.total([0, 1], 0, 0, lo=YNM)
tgt = L.vec([(0, 0, 0), (1, 0, 0)])
clo, _ = L.lp(tgt, 1); chi, _ = L.lp(tgt, -1)
lattice = [k for k in range(0, 3) if clo - 1e-9 <= k <= chi + 1e-9]
ints = [k for k in range(0, 3) if L.ip(tgt, 1, extra_eq=(tgt, k))[0] is not None]
report("RC5 C2(b): integer set is NOT 'the lattice inside the continuous interval' with fractional boxes", "CE",
       lattice == [1, 2] and ints == [2],
       f"two units (N=5,V=1,yes=1), ordering boxes n>=0.2, between sign Y>=Y^NM=0.4: continuous Y range [{clo:.1f},{chi:.1f}], "
       f"integers in it {lattice}, integer values {ints}; C2(b)'s formula max(Y^-_rd, ceil Y^NM) = 2 is right, the 'So ...' sentence is not")

# RC11 [CE, round 2]: the corrected C2(b) says the integer gap set is the full lattice "without unit rows, or with
# integer boxes". Integer boxes are not enough: fractional TOTALS on disjoint multi-unit sets (class (i) unit rows)
# have the same effect as T5's fractional boxes. The right condition is integer right-hand sides for all unit rows.
L = Sys(a=[[1, 4]] * 4, b=[[1, 4]] * 4)                      # four units N=5, V=1, yes=1; no boxes at all
L.total([0, 1], 0, 0, lo=0.2); L.total([2, 3], 0, 0, lo=0.2)  # two fractional totals on disjoint sets (one laminar family)
YNM = 4 * (1 * 1 / 5); L.total([0, 1, 2, 3], 0, 0, lo=YNM)  # between sign over all units: Y >= 0.8
tgt = L.vec([(i, 0, 0) for i in range(4)])
clo, _ = L.lp(tgt, 1); chi, _ = L.lp(tgt, -1)
latt = [k for k in range(0, 5) if clo - 1e-9 <= k <= chi + 1e-9]
ints = [k for k in range(0, 5) if L.ip(tgt, 1, extra_eq=(tgt, k))[0] is not None]
report("RC11 corrected C2(b): 'with integer boxes' is not enough; fractional totals also shrink the lattice", "CE",
       latt == [1, 2, 3, 4] and ints == [2, 3, 4],
       f"four units (N=5,V=1,yes=1), no boxes, totals y0+y1 >= 0.2 and y2+y3 >= 0.2, between sign Y >= {YNM:.1f}: "
       f"continuous range [{clo:.1f},{chi:.1f}] contains {latt}; integer values {ints}; condition should read "
       f"'all unit rows with integer right-hand sides'")

# RC5b [CE]: C2(a)/(b) need all sign rows on one unit set S and the target Y_rc(S): sign set != target set
L = Sys(a=[[2, 2], [2, 2]], b=[[2, 2], [2, 2]])
YNM0 = 2 * 2 / 4                                             # sign on unit set S = {0} only: Y_11({0}) >= 1
L.total([0], 0, 0, lo=YNM0)
tgt = L.vec([(0, 0, 0), (1, 0, 0)])                          # target Y_11 over all units
lo_, _ = L.rounded().lp(tgt, 1)
YNMall = 2 * (2 * 2 / 4)
report("RC5b C2: thresholds act on Y_rc(S) of the sign set, not on Y_rc of the target set", "CE",
       abs(lo_ - 1) < 1e-9 and YNMall == 2,
       f"two units (N=4,V=2,yes=2), between sign over S={{0}}: lower end of Y_11(all) = {lo_:.0f}, while reading C2(b) with "
       f"Y^NM computed over all units gives {YNMall:.0f}; C2 must fix one sign set S and the target Y_rc(S)")

# ===========================================================================
print("=== Lemma 2.1 / §9 item 3 ===")
# RC6 [V/CE for C2 contract]: the one-vs-rest gap D^r_c (R >= 3), written over the rest's cells, is a
# 'several same-sign cells of one column' target: the package's cellwise DD is outer; the exact range is
# the single-cell DD range mapped by D = kappa*Y_rc - Y_c/(N-A_r).
worst = 0.0; nwide = 0
for _ in range(200):
    I, R, C = int(rng.integers(2, 6)), int(rng.integers(3, 5)), 2
    tabs = [rng.integers(0, 8, size=(R, C)) + (rng.random((R, C)) < .2) for _ in range(I)]
    a = np.array([t.sum(1) for t in tabs], float); b = np.array([t.sum(0) for t in tabs], float); Ni = a.sum(1)
    if np.any(Ni == 0): continue
    A = a.sum(0); Nt = Ni.sum(); r, c = 0, 0
    lo = np.array([[np.sum(np.maximum(0, a[:, s] + b[:, k] - Ni)) / A[s] for k in range(C)] for s in range(R)])
    hi = np.array([[np.sum(np.minimum(a[:, s], b[:, k])) / A[s] for k in range(C)] for s in range(R)])
    w = np.zeros((R, C)); w[r, c] = 1; w[1:, c] = -A[1:] / (Nt - A[r])         # D = B_rc - sum_s A_s B_sc / (N - A_r)
    ia = (np.sum(np.where(w >= 0, w * lo, w * hi)), np.sum(np.where(w >= 0, w * hi, w * lo)))
    kap = 1 / A[r] + 1 / (Nt - A[r]); Yc = b[:, c].sum()
    ex = (kap * lo[r, c] * A[r] - Yc / (Nt - A[r]), kap * hi[r, c] * A[r] - Yc / (Nt - A[r]))
    wd = (ia[1] - ia[0]) - (ex[1] - ex[0]); worst = max(worst, wd); nwide += wd > 1e-9
    assert ia[0] <= ex[0] + 1e-12 and ia[1] >= ex[1] - 1e-12
report("RC6 one-vs-rest gap D^r_c (R >= 3): cellwise interval arithmetic is outer; exact = single-cell DD", "V",
       nwide > 0,
       f"200 random designs (R = 3, 4): cellwise strictly wider in {nwide}, max excess width {worst:.3f}; "
       f"6a §9 item 3 / contract C2 should name D^r_c (the M-SIGN default contrast) explicitly")

# ===========================================================================
print("=== Lemma 4.1 ===")
# RC7 [CE for 'S_phi is a property of phi']: the allowance depends on the representation of the target
# (adding multiples of the column-margin identity), and representations using outcome margins in
# phi_0 change E_phi and B~_phi. Both representations below are valid (each Var <= its S^2), with different S.
I = 400; R, C = 3, 2
Ni = rng.integers(50, 400, size=I); X = rng.dirichlet(np.ones(R) * 2, size=I); a = Ni[:, None] * X
A = a.sum(0); Nt = Ni.sum(); r = 0
kap = 1 / A[r] + 1 / (Nt - A[r])
rho1 = a[:, r] * kap                                                   # D = kappa*Y_rc - Y_c/(N-A_r)
rho2 = a[:, r] / A[r] + a[:, 1:].sum(1) / (Nt - A[r])                  # D = Y_rc/A_r - sum_{s!=r} Y_sc/(N-A_r)
S1, S2 = .5 * np.sqrt(np.sum(rho1 ** 2)), .5 * np.sqrt(np.sum(rho2 ** 2))
m = np.clip(0.3 + 0.4 * X, 0.05, 0.95)                                 # model means per unit and group (column 1)
v = m * (1 - m)                                                        # exact Var of a Bernoulli(m) unit rate (extreme law)
V1 = kap ** 2 * np.sum(a[:, r] ** 2 * v[:, r])
V2 = np.sum(a[:, r] ** 2 * v[:, r]) / A[r] ** 2 + np.sum(a[:, 1:] ** 2 * v[:, 1:]) / (Nt - A[r]) ** 2
beta = (rng.random((I, R)) < m).astype(float); e = beta - m          # one draw: the two E's are different numbers
E1 = kap * np.sum(a[:, r] * e[:, r]); E2 = np.sum(a[:, r] * e[:, r]) / A[r] - np.sum(a[:, 1:] * e[:, 1:]) / (Nt - A[r])
report("RC7 Lemma 4.1: S_phi depends on how the target is written (column-margin identity), E_phi too", "CE",
       abs(S1 - S2) / S2 > 0.05 and V1 <= S1 ** 2 and V2 <= S2 ** 2 and abs(E1 - E2) > 1e-4,
       f"R=3 one-vs-rest gap, 400 units: S_phi = {S1:.5f} (kappa*Y_rc - Y_c/(N-A_r)) vs {S2:.5f} (rest written over cells); "
       f"exact sd(E) under independent Bernoulli unit rates {np.sqrt(V1):.5f} vs {np.sqrt(V2):.5f} (each within its own bound); "
       f"in one draw E = {E1:.4f} vs {E2:.4f}: they differ by (Y_c - E[Y_c|design])/(N-A_r), so E_phi, B~_phi and S_phi "
       f"depend on the representation when phi_0 uses outcome margins")

# ===========================================================================
print("=== R27 ===")
# RC8 [CE for the probability reading of 'under H_C and H_M']: a count row that is a BELIEF about the
# realised table is a random event under H_M. Conditional on it holding, neither coverage nor size is
# controlled; unconditionally the bound is alpha + P(row fails).
I = 300; Ni = rng.integers(80, 200, size=I); x = rng.uniform(.1, .9, size=I); V = np.round(Ni * x).astype(int)
m1 = 0.3 + 0.4 * x                                                     # LCE mean for group 1 (pooled 2x2), m2 irrelevant
sdE = np.sqrt(np.sum((V / V.sum()) ** 2 * m1 * (1 - m1)))            # oracle sd of E for Bernoulli(m) unit rates
Btil = np.sum(V / V.sum() * m1); z = 1.959964
def rc8(kt, nrep=20000):
    t = Btil + kt * sdE                                                # belief row: Y_11/A_1 >= t (a stated floor on the realised count)
    cov_u = cov_c = rej_c = nrow = 0
    for _ in range(nrep):
        beta1 = (rng.random(I) < m1).astype(float)                     # each unit's group-1 rate is 0 or 1 (mean m1)
        B = np.sum(V * beta1) / V.sum()
        CI = (Btil - z * sdE, Btil + z * sdE)                         # exact-level oracle actual-target interval
        inter = (max(CI[0], t), CI[1]); empty = inter[0] > inter[1]   # count set C = [t, DD upper], DD upper >= B
        hit = (not empty) and inter[0] <= B <= inter[1]
        cov_u += hit
        if B >= t: nrow += 1; cov_c += hit; rej_c += empty
    return nrow / nrep, cov_c / max(nrow, 1), rej_c / max(nrow, 1), 1 - cov_u / nrep
pA, cA, rA, mA = rc8(1.0); pB, cB, rB, mB = rc8(2.4)
se = 3 * np.sqrt(.25 / 20000)
report("RC8 R27(i)-(ii): with belief rows, 'under H_C and H_M' must be read unconditionally with + P(H_C fails)", "CE",
       cA < 0.92 and rB > 0.9 and mA <= 0.05 + (1 - pA) + se and mB <= 0.05 + (1 - pB) + se,
       f"pooled 2x2 LCE true, oracle 95% interval for B; belief row 'B >= t'. t = B~ + 1.0 sd: P(row) = {pA:.3f}, coverage of "
       f"C∩CI given the row {cA:.3f} (< .95); t = B~ + 2.4 sd: P(row) = {pB:.4f}, rejection rate given the row {rB:.3f} (> .05). "
       f"Unconditional miss rates {mA:.3f}, {mB:.3f} obey alpha + P(row fails) = {0.05+1-pA:.3f}, {0.05+1-pB:.3f}")

# ===========================================================================
print("=== §1.4 rounding: numerical hazard (contract C4) ===")
# RC9 [B]: thresholds that are exact integers in rational arithmetic can round the wrong way in floating
# point (Y^ER involves several divisions); the contract should require exact rationals or a stated tolerance.
nhaz = ntry = 0
for _ in range(20000):
    I = int(rng.integers(2, 5)); Ns = rng.integers(2, 12, size=I); Vs = np.array([rng.integers(1, n) for n in Ns])
    Ys = np.array([rng.integers(0, n + 1) for n in Ns])
    N = int(Ns.sum()); p = [Fr(int(n), N) for n in Ns]; Xr = [Fr(int(v), int(n)) for v, n in zip(Vs, Ns)]
    Tc = [Fr(int(y), int(n)) for y, n in zip(Ys, Ns)]
    mu = sum(pi * xi for pi, xi in zip(p, Xr)); var = sum(pi * (xi - mu) ** 2 for pi, xi in zip(p, Xr))
    if var == 0: continue
    ER = sum(pi * (ti - sum(q * u for q, u in zip(p, Tc))) * (xi - mu) for pi, ti, xi in zip(p, Tc, Xr)) / var
    A = int(Vs.sum()); Yc = int(Ys.sum()); kap = Fr(1, A) + Fr(1, N - A)
    YER = (ER + Fr(Yc, N - A)) / kap
    pf = Ns / N; xf = Vs / Ns; tf = Ys / Ns; muf = np.sum(pf * xf); varf = np.sum(pf * (xf - muf) ** 2)
    ERf = np.sum(pf * (tf - np.sum(pf * tf)) * (xf - muf)) / varf
    YERf = (ERf + Yc / (N - A)) / (1 / A + 1 / (N - A))
    if YER.denominator == 1:
        ntry += 1; nhaz += math.floor(YERf) != int(YER)
report("RC9 floating-point floor of an integer-valued threshold Y^ER can be off by one", "B", ntry > 0,
       f"{ntry} random small populations with integer Y^ER (exact rationals): floor in double precision wrong in {nhaz}")

print(f"\n{sum(RESULTS)}/{len(RESULTS)} PASS")
