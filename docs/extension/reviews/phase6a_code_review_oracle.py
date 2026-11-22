#!/usr/bin/env python3
"""Independent oracle for the code review of eiCounts() (PHASE6a, M-COUNT).

Written without reading the package's row builders: every row is built here from its
definition (6a Sec. 1.3; 2c Definition 1 with realised rates beta = n / a), as an affine
function of the table whose coefficients are *extracted by evaluating the definition* at
0 and at unit vectors, in exact rational arithmetic (fractions.Fraction).

For each random tiny instance:
  * brute force: enumerate every integer table with the unit margins, keep those that
    satisfy every row exactly (CR = union of two systems), and record the exact integer
    min / max of each target (or integer infeasibility);
  * continuous endpoints: scipy.optimize.linprog (HiGHS) over the unrounded rows, per CR
    branch (an LP independent of lpSolve);
  * also the MILP endpoints with scipy.optimize.milp (cross-check of the brute force).

Usage:
  python3 phase6a_code_review_oracle.py gen   OUT.json [n_instances] [seed]
  python3 phase6a_code_review_oracle.py check INSTANCES.json RESULTS.json
The R side (phase6a_code_review_checks.R) reads OUT.json, runs eiCounts() on each
instance, and writes RESULTS.json (endpoints, statuses, witness tables); 'check'
compares and verifies each witness exactly against the rows built here.
"""
import itertools, json, math, random, sys
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds


# ---------------------------------------------------------------------------
# tables of one unit
# ---------------------------------------------------------------------------
def unit_tables(a, b):
    """all nonnegative integer R x C tables with row sums a and column sums b"""
    R, C = len(a), len(b)
    out = []
    def rec(r, colrem, acc):
        if r == R - 1:
            row = list(colrem)
            if sum(row) == a[r] and all(v >= 0 for v in row):
                out.append(acc + [row])
            return
        # compositions of a[r] into C parts bounded by colrem
        def comp(c, rem, cur):
            if c == C - 1:
                if rem <= colrem[c]:
                    yield cur + [rem]
                return
            for v in range(0, min(rem, colrem[c]) + 1):
                yield from comp(c + 1, rem - v, cur + [v])
        for row in comp(0, a[r], []):
            rec(r + 1, [colrem[c] - row[c] for c in range(C)], acc + [row])
    rec(0, list(b), [])
    return out


class Design:
    def __init__(self, A, B):
        self.A = [list(map(int, r)) for r in A]
        self.B = [list(map(int, r)) for r in B]
        self.p = len(A); self.R = len(A[0]); self.C = len(B[0])
        self.N = [sum(r) for r in self.A]
        self.nv = self.p * self.R * self.C

    def k(self, i, r, c):
        return (i * self.R + r) * self.C + c

    def beta(self, x, i, r, c):
        a = self.A[i][r]
        return F(x[self.k(i, r, c)], a) if a > 0 else F(0)


# affine extraction: f(x) = c0 + sum coef_k x_k, f evaluated exactly
def affine(f, nv):
    z = [F(0)] * nv
    c0 = f(z)
    coef = {}
    for k in range(nv):
        e = list(z); e[k] = F(1)
        v = f(e) - c0
        if v != 0:
            coef[k] = v
    return coef, c0


# ---------------------------------------------------------------------------
# definitions (2c Definition 1, realised rates), over a unit set S (0-based)
# ---------------------------------------------------------------------------
def pvec(D, S):
    Ntot = sum(D.N[i] for i in S)
    return {i: F(D.N[i], Ntot) for i in S}

def X(D, i, r):
    return F(D.A[i][r], D.N[i])

def wcov(w, f, g):
    sw = sum(w.values())
    if sw == 0:
        return F(0)
    mf = sum(w[i] * f[i] for i in w) / sw
    mg = sum(w[i] * g[i] for i in w) / sw
    return sum(w[i] * (f[i] - mf) * (g[i] - mg) for i in w) / sw

def delta_B_pair(D, S, r, s, c):
    p = pvec(D, S)
    return lambda x: sum(p[i] * X(D, i, r) * X(D, i, s) * (D.beta(x, i, r, c) - D.beta(x, i, s, c)) for i in S)

def delta_B_ovr(D, S, r, c):
    fs = [delta_B_pair(D, S, r, s, c) for s in range(D.R) if s != r]
    return lambda x: sum(f(x) for f in fs)

def delta_W_W1(D, S, r, c):
    p = pvec(D, S)
    mu = sum(p[i] * X(D, i, r) for i in S)
    def f(x):
        xr = {i: X(D, i, r) for i in S}
        br = {i: D.beta(x, i, r, c) for i in S}
        wr = {i: p[i] * xr[i] / mu for i in S}
        bn = {}
        for i in S:
            if xr[i] < 1:
                bn[i] = sum(X(D, i, s) * D.beta(x, i, s, c) for s in range(D.R) if s != r) / (1 - xr[i])
            else:
                bn[i] = F(0)
        wn = {i: p[i] * (1 - xr[i]) / (1 - mu) for i in S}
        return mu * wcov(wr, br, xr) + (1 - mu) * wcov(wn, bn, xr)
    return f

def delta_W_W2(D, S, r, c):
    p = pvec(D, S)
    def f(x):
        tot = F(0)
        xr = {i: X(D, i, r) for i in S}
        for s in range(D.R):
            mus = sum(p[i] * X(D, i, s) for i in S)
            if mus == 0:
                continue
            ws = {i: p[i] * X(D, i, s) / mus for i in S}
            bs = {i: D.beta(x, i, s, c) for i in S}
            tot += mus * wcov(ws, xr, bs)
        return tot
    return f

# slack scales: the package's probability-unit convention (deviation 6), over S
def scale_of(D, S, kind, r, other=None):
    p = pvec(D, S)
    xr = {i: X(D, i, r) for i in S}
    mu = {g: sum(p[i] * X(D, i, g) for i in S) for g in range(D.R)}
    def wvar(w):
        sw = sum(w.values())
        w = {i: w[i] / sw for i in w}
        return sum(w[i] * xr[i] ** 2 for i in w) - sum(w[i] * xr[i] for i in w) ** 2
    if kind == "dB_W1":
        return sum(p[i] * xr[i] * (1 - xr[i]) for i in S)
    if kind == "dB_pair":
        return sum(p[i] * xr[i] * X(D, i, other) for i in S)
    if kind == "dW_W1":
        return mu[r] * wvar({i: p[i] * xr[i] for i in S}) + (1 - mu[r]) * wvar({i: p[i] * (1 - xr[i]) for i in S})
    if kind == "dW_W2":
        return sum(mu[g] * wvar({i: p[i] * X(D, i, g) for i in S}) for g in range(D.R) if mu[g] > 0)
    raise ValueError(kind)


# ---------------------------------------------------------------------------
# restrictions -> alternatives (list of systems; each system a list of rows)
# a row: (coef dict, c0, lo, hi) meaning lo <= c0 + coef.x <= hi (None = open)
# a row may also be a predicate ("pred", fn) for brute force, with LP bounds given separately
# ---------------------------------------------------------------------------
def rows_of(D, spec):
    t = spec["type"]
    U = [u - 1 for u in spec.get("units", list(range(1, D.p + 1)))]
    if t in ("known", "zero", "box", "total", "order_unit", "sign", "cr"):
        r = spec["row"] - 1; c = spec["col"] - 1
    if t == "known":
        return [[({D.k(i, r, c): F(1)}, F(0), F(spec["value"]), F(spec["value"])) for i in U]]
    if t == "zero":
        return [[({D.k(i, r, c): F(1)}, F(0), F(0), F(0)) for i in U]]
    if t == "box":
        lo = None if spec.get("lower") is None else F(spec["lower"])
        hi = None if spec.get("upper") is None else F(spec["upper"])
        return [[({D.k(i, r, c): F(1)}, F(0), lo, hi) for i in U]]
    if t == "total":
        if spec.get("value") is not None:
            lo = hi = F(spec["value"])
        else:
            lo = None if spec.get("lower") is None else F(spec["lower"])
            hi = None if spec.get("upper") is None else F(spec["upper"])
        return [[({D.k(i, r, c): F(1) for i in U}, F(0), lo, hi)]]
    if t == "order_unit":
        sl = F(spec.get("slack", 0)); sg = spec.get("sign", 1)
        rows = []
        for i in U:
            if spec.get("over") is None:
                if D.A[i][r] == 0 or D.A[i][r] == D.N[i]:
                    continue
                # beta_r - beta_notr >= -slack   (sign +1)
                f = lambda x, i=i: (D.beta(x, i, r, c) - sum(F(x[D.k(i, s, c)]) for s in range(D.R) if s != r) / (D.N[i] - D.A[i][r]))
            else:
                q = spec["over"] - 1
                if D.A[i][r] == 0 or D.A[i][q] == 0:
                    continue
                f = lambda x, i=i, q=q: D.beta(x, i, r, c) - D.beta(x, i, q, c)
            co, c0 = affine(f, D.nv)
            if sg < 0:
                co = {k: -v for k, v in co.items()}; c0 = -c0
            rows.append((co, c0, -sl, None))
        return [rows]
    if t == "linear":
        co = {}
        for (u, rr, cc, v) in spec["cells"]:
            k = D.k(u - 1, rr - 1, cc - 1); co[k] = co.get(k, F(0)) + F(v)
        lo = None if spec.get("lower") is None else F(spec["lower"])
        hi = None if spec.get("upper") is None else F(spec["upper"])
        return [[(co, F(0), lo, hi)]]
    if t == "odds":
        rows = []
        psiL, psiU = F(spec["psi"][0]), F(spec["psi"][1])
        for i in U:
            a, b, N = D.A[i][0], D.B[i][0], D.N[i]
            if a in (0, N) or b in (0, N):
                continue
            rows.append(("odds", i, psiL, psiU))
        return [rows]
    if t in ("sign", "cr"):
        sl = F(spec.get("slack", 0))
        ver = spec.get("version", "W1")
        other = spec.get("other"); other = None if other is None else other - 1
        def mk(tp, v, sg):
            if tp == "dB":
                f = delta_B_ovr(D, U, r, c) if v == "W1" else delta_B_pair(D, U, r, other, c)
                sc = scale_of(D, U, "dB_W1" if v == "W1" else "dB_pair", r, other)
            else:
                f = delta_W_W1(D, U, r, c) if v == "W1" else delta_W_W2(D, U, r, c)
                sc = scale_of(D, U, "dW_W1" if v == "W1" else "dW_W2", r)
            co, c0 = affine(f, D.nv)
            if sg < 0:
                co = {k: -x for k, x in co.items()}; c0 = -c0
            return (co, c0, -sl * sc, None)
        if t == "sign":
            return [[mk(spec.get("stype", "dB"), ver, spec.get("sign", 1))]]
        vB = "pair" if ver == "pair" else "W1"; vW = "W1" if ver == "W1" else "W2"
        return [[mk("dB", vB, 1), mk("dW", vW, 1)], [mk("dB", vB, -1), mk("dW", vW, -1)]]
    raise ValueError(t)


def odds_ok(D, x, i, psiL, psiU):
    n11 = x[D.k(i, 0, 0)]; n12 = x[D.k(i, 0, 1)]; n21 = x[D.k(i, 1, 0)]; n22 = x[D.k(i, 1, 1)]
    num = n11 * n22; den = n12 * n21
    if den == 0:
        # ratio infinite (num > 0) or 0/0 at a corner; endpoint of the range: in the box iff psiU = inf
        if num == 0:
            return True  # both 0: cannot happen with all margins strictly inside; treat as unconstrained
        return math.isinf(float(psiU))
    if num == 0:
        return psiL == 0
    rt = F(num, den)
    return psiL <= rt <= psiU


def odds_cont_box(D, i, psiL, psiU):
    a, b, N = D.A[i][0], D.B[i][0], D.N[i]
    L, Uu = max(0, a + b - N), min(a, b)
    def g(k, psi):
        return k * (N - a - b + k) - psi * (a - k) * (b - k)
    def root(psi):
        lo, hi = float(L), float(Uu)
        for _ in range(200):
            m = (lo + hi) / 2
            if g(m, psi) < 0: lo = m
            else: hi = m
        return (lo + hi) / 2
    lo = root(float(psiL)) if psiL > 0 else None
    hi = root(float(psiU)) if not math.isinf(float(psiU)) else None
    return lo, hi


def row_ok(D, row, x):
    if row[0] == "odds":
        return odds_ok(D, x, row[1], row[2], row[3])
    co, c0, lo, hi = row
    v = c0 + sum(co[k] * x[k] for k in co)
    return (lo is None or v >= lo) and (hi is None or v <= hi)


# ---------------------------------------------------------------------------
# targets: coefficient vectors (exact) on n, plus constant
# ---------------------------------------------------------------------------
def target_vec(D, tg):
    t = tg["kind"]
    S = [u - 1 for u in tg.get("units", list(range(1, D.p + 1)))]
    co = {}
    def add(Su, rows, cols, w):
        for i in Su:
            for r in rows:
                for c in cols:
                    k = D.k(i, r, c); co[k] = co.get(k, F(0)) + w
    if t == "rate":
        rows = [r - 1 for r in tg["row"]]; cols = [c - 1 for c in tg["col"]]
        den = sum(D.A[i][r] for i in S for r in rows)
        add(S, rows, cols, F(1, den))
    elif t == "gap":
        S2 = [u - 1 for u in tg.get("units2", tg.get("units", list(range(1, D.p + 1))))]
        r1, r2, c = tg["row1"] - 1, tg["row2"] - 1, tg["col"] - 1
        add(S, [r1], [c], F(1, sum(D.A[i][r1] for i in S)))
        add(S2, [r2], [c], F(-1, sum(D.A[i][r2] for i in S2)))
    elif t == "ovr_gap":
        r, c = tg["row"] - 1, tg["col"] - 1
        Ar = sum(D.A[i][r] for i in S); N = sum(D.N[i] for i in S)
        add(S, [r], [c], F(1, Ar))
        add(S, [s for s in range(D.R) if s != r], [c], F(-1, N - Ar))
    elif t == "share":
        r, c = tg["row"] - 1, tg["col"] - 1
        add(S, [r], [c], F(1, sum(D.B[i][c] for i in S)))
    else:
        raise ValueError(t)
    return co


def tval(co, x):
    return sum(co[k] * x[k] for k in co)


# ---------------------------------------------------------------------------
# solve one instance
# ---------------------------------------------------------------------------
def systems(D, specs):
    alts = [rows_of(D, s) for s in specs]
    return [sum(choice, []) for choice in itertools.product(*alts)] if alts else [[]]


def lp_bounds(D, rows):
    """rows -> A_ub, b_ub, A_eq, b_eq for scipy (floats; odds as continuous boxes)"""
    Aub, bub, Aeq, beq = [], [], [], []
    for i in range(D.p):
        for r in range(D.R):
            v = np.zeros(D.nv); v[[D.k(i, r, c) for c in range(D.C)]] = 1; Aeq.append(v); beq.append(D.A[i][r])
        for c in range(D.C):
            v = np.zeros(D.nv); v[[D.k(i, r, c) for r in range(D.R)]] = 1; Aeq.append(v); beq.append(D.B[i][c])
    for row in rows:
        if row[0] == "odds":
            i = row[1]; lo, hi = odds_cont_box(D, i, row[2], row[3])
            v = np.zeros(D.nv); v[D.k(i, 0, 0)] = 1
            if lo is not None: Aub.append(-v); bub.append(-lo)
            if hi is not None: Aub.append(v); bub.append(hi)
            continue
        co, c0, lo, hi = row
        v = np.zeros(D.nv)
        for k, x in co.items(): v[k] = float(x)
        if lo is not None and hi is not None and lo == hi:
            Aeq.append(v); beq.append(float(lo - c0)); continue
        if lo is not None: Aub.append(-v); bub.append(-float(lo - c0))
        if hi is not None: Aub.append(v); bub.append(float(hi - c0))
    return Aub, bub, Aeq, beq


def lp_range(D, rows, co, integer=False):
    Aub, bub, Aeq, beq = lp_bounds(D, rows)
    c = np.zeros(D.nv)
    for k, x in co.items(): c[k] = float(x)
    res = []
    for sgn in (1, -1):
        if not integer:
            r = linprog(sgn * c, A_ub=np.array(Aub) if Aub else None, b_ub=bub or None,
                        A_eq=np.array(Aeq), b_eq=beq, bounds=[(0, None)] * D.nv, method="highs")
            if r.status != 0:
                return None
            res.append(sgn * r.fun)
        else:
            cons = [LinearConstraint(np.array(Aeq), beq, beq)]
            if Aub: cons.append(LinearConstraint(np.array(Aub), -np.inf, bub))
            r = milp(sgn * c, constraints=cons, integrality=np.ones(D.nv), bounds=Bounds(0, np.inf))
            if r.status != 0:
                return None
            res.append(sgn * r.fun)
    return (res[0], res[1])


def brute(D, sysl, cos):
    per = [unit_tables(D.A[i], D.B[i]) for i in range(D.p)]
    n_all = 1
    for t in per: n_all *= len(t)
    lo = [None] * len(cos); hi = [None] * len(cos); feas = 0
    for combo in itertools.product(*per):
        x = [v for tab in combo for row in tab for v in row]
        ok = any(all(row_ok(D, row, x) for row in rows) for rows in sysl)
        if not ok:
            continue
        feas += 1
        for j, co in enumerate(cos):
            v = tval(co, x)
            if lo[j] is None or v < lo[j]: lo[j] = v
            if hi[j] is None or v > hi[j]: hi[j] = v
    return lo, hi, feas, n_all


# ---------------------------------------------------------------------------
# random instances
# ---------------------------------------------------------------------------
def rand_design(rng, p, R, C, Nmax):
    A, B = [], []
    for _ in range(p):
        N = rng.randint(2, Nmax)
        a = rand_comp(rng, N, R); b = rand_comp(rng, N, C)
        A.append(a); B.append(b)
    return A, B

def rand_comp(rng, N, k):
    # mostly positive parts (degenerate designs make most rows trivially infeasible)
    if N >= k and rng.random() < 0.85:
        cuts = sorted(rng.sample(range(1, N), k - 1))
    else:
        cuts = sorted(rng.randint(0, N) for _ in range(k - 1))
    return [y - x for x, y in zip([0] + cuts, cuts + [N])]

def frac(rng, lo, hi):
    return round(rng.uniform(lo, hi), 2)

def laminar_pair(p):
    # family A (nested / disjoint) and family B crossing A
    A = [list(range(1, p + 1)), list(range(1, p // 2 + 1)) if p >= 2 else [1]]
    B = [list(range(2, p + 1))] if p >= 3 else []
    return A, B


def gen_instance(rng, scenario):
    """returns dict(design, specs, targets, scenario)"""
    if scenario in ("tu_i", "tu_i_cross3", "cr_w1", "sign_w1"):
        R, C = rng.choice([(2, 2), (2, 3), (3, 2), (3, 3)])
    elif scenario in ("tu_ii", "tu_ii_t", "mix", "pairwise", "w2", "pair_dB", "linear"):
        R, C = (3, 3) if scenario in ("mix",) else rng.choice([(3, 2), (3, 3), (2, 3)])
        if scenario in ("pairwise", "w2", "pair_dB"):
            R = 3; C = rng.choice([2, 3])
    elif scenario == "tu_iii":
        R, C = 2, rng.choice([2, 3])
    elif scenario == "tu_iii_t":
        R, C = rng.choice([2, 3]), 2
    elif scenario == "odds":
        R, C = 2, 2
    else:
        R, C = 2, 2
    p = 4 if scenario in ("tu_i", "tu_i_cross3", "tu_ii", "tu_iii", "tu_iii_t", "tu_ii_t") and R * C <= 4 else 3
    Nmax = 7 if R * C <= 4 else (5 if R * C <= 6 else 4)
    if p == 4 and R * C > 4:
        p = 3
    A, B = rand_design(rng, p, R, C, Nmax)
    D = Design(A, B)
    specs = []
    allu = list(range(1, p + 1))
    def cell():
        return rng.randint(1, R), rng.randint(1, C)
    def tot(r, c, U):
        sl = sum(1 for _ in U)
        cap = sum(min(D.A[i - 1][r - 1], D.B[i - 1][c - 1]) for i in U)
        u = rng.random()
        if u < 0.4:
            return dict(type="total", row=r, col=c, units=U, lower=round(rng.uniform(0, 0.6) * cap, 2))
        if u < 0.8:
            return dict(type="total", row=r, col=c, units=U, upper=round(rng.uniform(0.4, 1.0) * cap + 0.01, 2))
        return dict(type="total", row=r, col=c, units=U, value=rng.randint(0, cap))
    # common unit rows sometimes: fractional box, known cell
    if rng.random() < 0.6:
        r, c = cell(); i = rng.randint(1, p)
        cap = min(D.A[i - 1][r - 1], D.B[i - 1][c - 1])
        specs.append(dict(type="box", row=r, col=c, units=[i], lower=round(rng.uniform(0, 0.8) * cap + 0.1, 2)))
    if rng.random() < 0.3:
        r, c = cell(); i = rng.randint(1, p)
        # known value = some feasible value: pick within Frechet bounds
        a, b, N = D.A[i - 1][r - 1], D.B[i - 1][c - 1], D.N[i - 1]
        specs.append(dict(type="known", row=r, col=c, units=[i], value=rng.randint(max(0, a + b - N), min(a, b))))
    if scenario == "tu_i":
        r, c = cell()
        fams = [allu, allu[:2], allu[2:]] + ([allu[1:3]] if p >= 3 else [])
        for U in fams:
            if len(U) >= 2 and rng.random() < 0.7:
                specs.append(tot(r, c, U))
        if rng.random() < 0.5:
            specs.append(dict(type="order_unit", row=r, col=c, units=allu, sign=rng.choice([1, -1])))
    elif scenario == "tu_i_cross3":
        r, c = cell()
        for U in ([1, 2], [2, 3], [1, 3]):
            specs.append(tot(r, c, U))
    elif scenario == "tu_ii":
        c = rng.randint(1, C)
        for r in range(1, R + 1):
            if rng.random() < 0.7:
                specs.append(tot(r, c, rng.choice([allu, allu[:2], allu[1:]])))
    elif scenario == "tu_ii_t":
        r = rng.randint(1, R)
        for c in range(1, C + 1):
            if rng.random() < 0.7:
                specs.append(tot(r, c, rng.choice([allu, allu[:2], allu[1:]])))
    elif scenario == "tu_iii":
        for c in range(1, C + 1):
            fam = rng.choice([[allu], [allu[:2]], [allu, allu[:2]], [allu[1:]]])
            for U in fam:
                specs.append(tot(rng.randint(1, 2), c, U))
    elif scenario == "tu_iii_t":
        for r in range(1, R + 1):
            fam = rng.choice([[allu], [allu[:2]], [allu, allu[:2]], [allu[1:]]])
            for U in fam:
                specs.append(tot(r, rng.randint(1, 2), U))
    elif scenario == "mix":
        specs.append(tot(1, 1, [1, 2])); specs.append(tot(2, 2, [2, 3]))
        if rng.random() < 0.5: specs.append(tot(1, 2, allu))
    elif scenario == "pairwise":
        r, c = cell(); q = rng.choice([x for x in range(1, R + 1) if x != r])
        specs.append(dict(type="order_unit", row=r, col=c, over=q, units=allu, sign=rng.choice([1, -1])))
    elif scenario == "pair_dB":
        r, c = cell(); q = rng.choice([x for x in range(1, R + 1) if x != r])
        specs.append(dict(type="sign", stype="dB", version="pair", row=r, other=q, col=c, units=allu, sign=rng.choice([1, -1])))
    elif scenario == "w2":
        r, c = cell()
        specs.append(dict(type="sign", stype="dW", version="W2", row=r, col=c, units=allu, sign=rng.choice([1, -1])))
    elif scenario == "linear":
        cells = []
        for _ in range(3):
            cells.append([rng.randint(1, p), rng.randint(1, R), rng.randint(1, C), rng.choice([1, 2, -1, 0.5, 1.5])])
        tgt = frac(rng, -1, 2)
        specs.append(dict(type="linear", cells=cells, lower=tgt))
    elif scenario == "sign_w1":
        r, c = cell()
        stp = rng.choice(["dB", "dW"])
        specs.append(dict(type="sign", stype=stp, version="W1", row=r, col=c, units=allu, sign=rng.choice([1, -1]),
                          slack=rng.choice([0, 0, 0.1])))
    elif scenario == "cr_w1":
        r, c = cell()
        specs.append(dict(type="cr", version="W1", row=r, col=c, units=allu, slack=rng.choice([0, 0, 0.05])))
    elif scenario == "odds":
        lo = rng.choice([0, 0.5, 1, 2]); hi = lo + rng.choice([0, 0.5, 2, 10])
        specs.append(dict(type="odds", psi=[lo, hi], units=allu))
    # targets
    targets = []
    r, c = cell()
    targets.append(dict(kind="rate", row=[r], col=[c]))
    if p >= 2:
        targets.append(dict(kind="rate", row=[r], col=[c], units=allu[:2]))
    if R >= 2:
        r2 = rng.choice([x for x in range(1, R + 1) if x != r])
        targets.append(dict(kind="gap", row1=r, row2=r2, col=c))
        targets.append(dict(kind="ovr_gap", row=r, col=c))
    if R >= 3:
        targets.append(dict(kind="rate", row=[1, 2], col=[c]))
    if C >= 3:
        targets.append(dict(kind="rate", row=[r], col=[1, 2]))
    targets.append(dict(kind="share", row=r, col=c))
    if p >= 3:
        targets.append(dict(kind="gap", row1=r, row2=r, col=c, units=[1, 2], units2=[3]))
    # drop degenerate targets (zero denominators)
    good = []
    for t in targets:
        try:
            target_vec(D, t); good.append(t)
        except ZeroDivisionError:
            pass
    # validate specs (dW needs Var > 0 etc.)
    okspecs = []
    for s in specs:
        try:
            rows_of(D, s); okspecs.append(s)
        except ZeroDivisionError:
            pass
    return dict(A=A, B=B, specs=okspecs, targets=good, scenario=scenario)


def fstr(x):
    return None if x is None else str(x)


def solve_instance(inst):
    D = Design(inst["A"], inst["B"])
    sysl = systems(D, inst["specs"])
    cos = [target_vec(D, t) for t in inst["targets"]]
    lo, hi, feas, n_all = brute(D, sysl, cos)
    cont = []; milpres = []
    for co in cos:
        best = None; bi = None
        for rows in sysl:
            r = lp_range(D, rows, co)
            if r is not None:
                best = r if best is None else (min(best[0], r[0]), max(best[1], r[1]))
            m = lp_range(D, rows, co, integer=True) if not any(rw[0] == "odds" for rw in rows) else None
            if m is not None:
                bi = m if bi is None else (min(bi[0], m[0]), max(bi[1], m[1]))
        cont.append(best); milpres.append(bi)
    return dict(int_lo=[fstr(v) for v in lo], int_hi=[fstr(v) for v in hi], n_feasible=feas, n_tables=n_all,
                cont=cont, milp=milpres)


SCEN = ["tu_i", "tu_i_cross3", "tu_ii", "tu_ii_t", "tu_iii", "tu_iii_t", "mix", "pairwise", "pair_dB",
        "w2", "linear", "sign_w1", "cr_w1", "odds"]


def cmd_gen(out, n, seed):
    rng = random.Random(seed)
    import os
    scen = os.environ.get("ORACLE_SCEN", ",".join(SCEN)).split(",")
    insts = []
    tries = 0
    while len(insts) < n and tries < 20 * n:
        tries += 1
        sc = scen[len(insts) % len(scen)]
        inst = gen_instance(rng, sc)
        D = Design(inst["A"], inst["B"])
        tot = 1
        for i in range(D.p): tot *= len(unit_tables(D.A[i], D.B[i]))
        if tot > 60000 or not inst["targets"]:
            continue
        inst["oracle"] = solve_instance(inst)
        inst["id"] = len(insts) + 1
        insts.append(inst)
    json.dump(insts, open(out, "w"), indent=0)
    print("wrote", len(insts), "instances to", out)


def frac_of(s):
    return None if s is None else F(s)


def cmd_check(inst_file, res_file):
    insts = {d["id"]: d for d in json.load(open(inst_file))}
    res = json.load(open(res_file))
    nbad = 0; nint = 0; ncont = 0; nwit = 0; nwit_bad = 0
    by_scen = {}
    issues = []
    for r in res:
        inst = insts[r["id"]]; D = Design(inst["A"], inst["B"]); o = inst["oracle"]
        sysl = systems(D, inst["specs"])
        sc = inst["scenario"]; by_scen.setdefault(sc, [0, 0])
        if r.get("error"):
            issues.append((r["id"], sc, "R error: " + r["error"])); nbad += 1; continue
        for j, t in enumerate(inst["targets"]):
            by_scen[sc][0] += 1
            lo, hi = frac_of(o["int_lo"][j]), frac_of(o["int_hi"][j])
            Rlo, Rhi = r["lower"][j], r["upper"][j]
            # integer endpoints
            if lo is None:
                if Rlo is not None:
                    issues.append((r["id"], sc, f"target {j+1}: oracle integer-infeasible, R gives [{Rlo}, {Rhi}] status {r['status'][j]}"))
                    by_scen[sc][1] += 1; nbad += 1
            else:
                nint += 1
                tol = 1e-7 if sc.startswith("large") else 1e-9
                if Rlo is None or abs(Rlo - float(lo)) > tol or abs(Rhi - float(hi)) > tol:
                    issues.append((r["id"], sc, f"target {j+1}: integer [{float(lo):.6g}, {float(hi):.6g}] vs R [{Rlo}, {Rhi}] status {r['status'][j]} class {r['class'][j]}"))
                    by_scen[sc][1] += 1; nbad += 1
            # continuous endpoints
            ct = o["cont"][j]
            if ct is None:
                if r["lp_lower"][j] is not None:
                    issues.append((r["id"], sc, f"target {j+1}: oracle LP infeasible, R lp [{r['lp_lower'][j]}, {r['lp_upper'][j]}]"))
                    nbad += 1
            else:
                ncont += 1
                if r["lp_lower"][j] is None or abs(r["lp_lower"][j] - ct[0]) > 1e-7 or abs(r["lp_upper"][j] - ct[1]) > 1e-7:
                    issues.append((r["id"], sc, f"target {j+1}: LP [{ct[0]:.8g}, {ct[1]:.8g}] vs R lp [{r['lp_lower'][j]}, {r['lp_upper'][j]}]"))
                    nbad += 1
            # witnesses: exact check against the oracle's own rows
            for e, w in enumerate(r["witness"][j]):
                if w is None:
                    continue
                nwit += 1
                x = [F(v).limit_denominator(10**9) if not float(v).is_integer() else int(v) for v in w["x"]]
                integer = all(isinstance(v, int) for v in x)
                if not integer:
                    issues.append((r["id"], sc, f"target {j+1} end {e+1}: witness not integer")); nwit_bad += 1; continue
                ok_m = all(sum(x[D.k(i, rr, c)] for c in range(D.C)) == D.A[i][rr] for i in range(D.p) for rr in range(D.R)) and \
                       all(sum(x[D.k(i, rr, c)] for rr in range(D.R)) == D.B[i][c] for i in range(D.p) for c in range(D.C)) and \
                       all(v >= 0 for v in x)
                ok_r = any(all(row_ok(D, row, x) for row in rows) for rows in sysl)
                v = tval(target_vec(D, t), x)
                ok_v = abs(float(v) - (Rlo if e == 0 else Rhi)) < 1e-9
                if not (ok_m and ok_r and ok_v):
                    issues.append((r["id"], sc, f"target {j+1} end {e+1}: witness margins {ok_m} rows {ok_r} value {ok_v}"))
                    nwit_bad += 1
    print(f"instances: {len(res)}; integer endpoints compared: {nint}; continuous: {ncont}; witnesses: {nwit} (bad {nwit_bad}); mismatches: {nbad}")
    for sc, (n, b) in sorted(by_scen.items()):
        print(f"  {sc:12s} targets {n:4d}  integer mismatches {b}")
    for it in issues[:60]:
        print("  ISSUE", it)
    return nbad + nwit_bad


if __name__ == "__main__" and sys.argv[1] in ("gen", "check"):
    if sys.argv[1] == "gen":
        cmd_gen(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 280, int(sys.argv[4]) if len(sys.argv) > 4 else 20261122)
    else:
        sys.exit(1 if cmd_check(sys.argv[2], sys.argv[3]) else 0)


# ---------------------------------------------------------------------------
# larger instances (no enumeration): continuous endpoints by HiGHS LP, integer endpoints by
# HiGHS MILP, both over the oracle's own rows. Usage: ... genlarge OUT.json n seed
# ---------------------------------------------------------------------------
def gen_large(rng, kind):
    p = rng.randint(60, 150); R, C = 3, 3
    A, B = [], []
    for _ in range(p):
        N = rng.randint(20, 300)
        A.append(rand_comp(rng, N, R)); B.append(rand_comp(rng, N, C))
    D = Design(A, B); allu = list(range(1, p + 1))
    specs = []
    i = rng.randint(1, p)
    specs.append(dict(type="known", row=1, col=1, units=[i], value=rng.randint(max(0, A[i-1][0] + B[i-1][0] - D.N[i-1]), min(A[i-1][0], B[i-1][0]))))
    blocks = [allu[k:k + 10] for k in range(0, p, 10)]
    cap = lambda U, r, c: sum(min(A[u-1][r-1], B[u-1][c-1]) for u in U)
    if kind == "tu":
        for U in blocks[:5]:
            specs.append(dict(type="total", row=1, col=1, units=U, lower=round(0.3 * cap(U, 1, 1) + 0.37, 2)))
        specs.append(dict(type="sign", stype="dB", version="W1", row=1, col=1, units=allu, sign=1))
        specs.append(dict(type="order_unit", row=2, col=1, units=allu[:30], sign=1))
    elif kind == "mix":
        for U in blocks[:4]:
            specs.append(dict(type="total", row=1, col=1, units=U, upper=round(0.8 * cap(U, 1, 1) + 0.37, 2)))
            specs.append(dict(type="total", row=2, col=2, units=U[3:] + [u + 3 for u in U[:3] if u + 3 <= p], lower=round(0.2 * cap(U, 2, 2) + 0.21, 2)))
    else:  # pairwise orderings in 20 units
        specs.append(dict(type="order_unit", row=1, col=1, over=2, units=allu[:20], sign=1))
    targets = [dict(kind="rate", row=[1], col=[1]), dict(kind="gap", row1=1, row2=2, col=1), dict(kind="rate", row=[1, 2], col=[1]),
               dict(kind="ovr_gap", row=1, col=1), dict(kind="rate", row=[1], col=[1, 2], units=allu[:40])]
    return dict(A=A, B=B, specs=specs, targets=targets, scenario="large_" + kind)


def cmd_genlarge(out, n, seed):
    rng = random.Random(seed); insts = []
    for k in range(n):
        inst = gen_large(rng, ["tu", "mix", "pairwise"][k % 3])
        D = Design(inst["A"], inst["B"]); sysl = systems(D, inst["specs"])
        cos = [target_vec(D, t) for t in inst["targets"]]
        cont = [lp_range(D, sysl[0], co) for co in cos]
        mi = [lp_range(D, sysl[0], co, integer=True) for co in cos]
        inst["oracle"] = dict(int_lo=[None if m is None else repr(m[0]) for m in mi], int_hi=[None if m is None else repr(m[1]) for m in mi],
                              cont=cont, milp=mi, n_feasible=None, n_tables=None)
        inst["id"] = k + 1; insts.append(inst)
    json.dump(insts, open(out, "w"))
    print("wrote", len(insts))


if __name__ == "__main__" and sys.argv[1] == "genlarge":
    cmd_genlarge(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
