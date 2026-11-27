"""Checks for PHASE2e_did_general.md (general difference-in-differences with ecological data).

Run:  python3 -I verify_phase2e.py > verify_phase2e_output.txt

Tags (as in 2c, 5a, 6a): [F] fails on an earlier or tempting ("naive") statement and passes on the one
made in PHASE2e; [V] verifies a statement of PHASE2e; [I] illustrates; [P2d] a check of PHASE2d, run
unchanged. Every check prints PASS or FAIL. Check names start with the result they concern (E1, ...,
G1, ..., W for the worked example). Sections IR and PT were added after the independent review
(reviews/PHASE2e_independent_review.md); they port the cases of reviews/phase2e_independent_checks.py and
reviews/phase2e_independent_checks_pt.py.

Conventions: periods t = 0..T-1 in code (event time in the document: base period t0 is the last
pre-period); rows r (groups), columns c (outcomes); n[t, i, r, c] = latent count; a[i, r] = group
counts (fixed population: the same in every period); b[t, i, c] = outcome counts.
"""
import itertools
import os
import runpy
import contextlib
import io

import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

rng = np.random.default_rng(20261122)
results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""), flush=True)


# =====================================================================================
# P. PHASE2d's 15 checks, run unchanged (its own seed), counted here.
# =====================================================================================
here = os.path.dirname(os.path.abspath(__file__))
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ns2d = runpy.run_path(os.path.join(here, "verify_phase2d.py"), run_name="verify_phase2d")
for line in buf.getvalue().splitlines():
    if line.startswith("PASS") or line.startswith("FAIL"):
        check("[P2d] " + line[6:], line.startswith("PASS"))


# =====================================================================================
# Helpers: the T-period count LP (fixed population).
# =====================================================================================
class CountLP:
    """Variables n[t, i, r, c] >= 0 with margins sum_c n = a[i, r] (every t), sum_r n = b[t, i, c].
    Extra rows: add_ge(coef, rhs) means coef . n >= rhs; add_eq likewise. Extra continuous variables
    (for the causal layer) can be appended with add_var."""

    def __init__(self, a, b):
        self.a = np.asarray(a, float)
        self.b = np.asarray(b, float)
        self.T, self.p, self.C = self.b.shape
        self.R = self.a.shape[1]
        self.nn = self.T * self.p * self.R * self.C
        self.nextra = 0
        self.extra_bounds = []
        self.Aeq, self.beq, self.Aub, self.bub = [], [], [], []
        self.lb = np.zeros(self.nn)
        self.ub = np.full(self.nn, np.inf)
        for t in range(self.T):
            for i in range(self.p):
                for r in range(self.R):
                    v = np.zeros(self.nn)
                    for c in range(self.C):
                        v[self.ix(t, i, r, c)] = 1
                    self._eq(v, self.a[i, r])
                for c in range(self.C):
                    v = np.zeros(self.nn)
                    for r in range(self.R):
                        v[self.ix(t, i, r, c)] = 1
                    self._eq(v, self.b[t, i, c])

    def ix(self, t, i, r, c):
        return ((t * self.p + i) * self.R + r) * self.C + c

    def _eq(self, v, k):
        self.Aeq.append(v)
        self.beq.append(k)

    def nvar(self):
        return self.nn + self.nextra

    def add_var(self, lo=None, hi=None):
        self.nextra += 1
        self.extra_bounds.append((lo, hi))
        return self.nn + self.nextra - 1

    def pad(self, v):
        v = np.asarray(v, float)
        if len(v) < self.nvar():
            v = np.concatenate([v, np.zeros(self.nvar() - len(v))])
        return v

    def add_ge(self, coef, rhs):
        self.Aub.append(-np.asarray(coef, float))
        self.bub.append(-rhs)

    def add_le(self, coef, rhs):
        self.Aub.append(np.asarray(coef, float))
        self.bub.append(rhs)

    def add_eq(self, coef, rhs):
        self._eq(np.asarray(coef, float), rhs)

    def _mats(self):
        n = self.nvar()
        Aeq = np.array([self.pad(v) for v in self.Aeq])
        Aub = np.array([self.pad(v) for v in self.Aub]) if self.Aub else None
        bounds = list(zip(self.lb, [None if np.isinf(u) else u for u in self.ub])) + self.extra_bounds
        return Aeq, np.array(self.beq), Aub, (np.array(self.bub) if self.Aub else None), bounds

    def opt(self, coef, sense=1):
        Aeq, beq, Aub, bub, bounds = self._mats()
        r = linprog(sense * self.pad(coef), A_eq=Aeq, b_eq=beq, A_ub=Aub, b_ub=bub, bounds=bounds,
                    method="highs")
        if r.status != 0:
            return None, None
        return sense * r.fun, r.x

    def range(self, coef, const=0.0):
        lo, _ = self.opt(coef, 1)
        hi, _ = self.opt(coef, -1)
        if lo is None or hi is None:
            return None
        return lo + const, hi + const

    def feasible(self):
        return self.opt(np.zeros(self.nvar()))[0] is not None

    def irange(self, coef, const=0.0):
        """Integer (MILP) range over integer n (extra variables continuous)."""
        Aeq, beq, Aub, bub, bounds = self._mats()
        n = self.nvar()
        cons = [LinearConstraint(Aeq, beq, beq)]
        if Aub is not None:
            cons.append(LinearConstraint(Aub, -np.inf, bub))
        lo = np.array([-np.inf if b[0] is None else b[0] for b in bounds])
        hi = np.array([np.inf if b[1] is None else b[1] for b in bounds])
        integ = np.concatenate([np.ones(self.nn), np.zeros(self.nextra)])
        out = []
        for s in (1, -1):
            r = milp(s * self.pad(coef), constraints=cons, bounds=Bounds(lo, hi), integrality=integ)
            if r.status != 0:
                return None
            out.append(s * r.fun)
        return out[0] + const, out[1] + const


def rate_vec(L, r, c, t, units=None):
    """B_{rct}(S) = Y_{rct}(S) / A_r(S) as a coefficient vector."""
    units = range(L.p) if units is None else units
    A = sum(L.a[i, r] for i in units)
    v = np.zeros(L.nn)
    for i in units:
        v[L.ix(t, i, r, c)] = 1 / A
    return v


def rest_vec(L, r, c, t):
    """B_{not r, c, t}: rate of everyone outside row r (written over cells)."""
    A = L.a.sum() - L.a[:, r].sum()
    v = np.zeros(L.nn)
    for i in range(L.p):
        for s in range(L.R):
            if s != r:
                v[L.ix(t, i, s, c)] = 1 / A
    return v


def did_vec(L, r, s, c, t0, t1):
    """D = [B_{r t1} - B_{r t0}] - [B_{s t1} - B_{s t0}] (s = 'rest' gives the one-vs-rest DiD)."""
    if s == "rest":
        g = lambda t: rate_vec(L, r, c, t) - rest_vec(L, r, c, t)
    else:
        g = lambda t: rate_vec(L, r, c, t) - rate_vec(L, s, c, t)
    return g(t1) - g(t0)


def rand_margins(p, R, C, T, Nlo=20, Nhi=60, integer=True):
    """Random integer fixed-population margins: a[i, r], b[t, i, c] from random latent tables."""
    N = rng.integers(Nlo, Nhi + 1, p)
    a = np.array([rng.multinomial(N[i], rng.dirichlet(np.ones(R))) for i in range(p)])
    n = np.zeros((T, p, R, C), int)
    for t in range(T):
        for i in range(p):
            for r in range(R):
                n[t, i, r] = rng.multinomial(a[i, r], rng.dirichlet(np.ones(C)))
    b = n.sum(2)
    return a, b, n


def rand_polar(p, T, C=2, share=(0.02, 0.98), Nlo=200, Nhi=400, R=2):
    """Near-pure units (group-0 share near share[0] or share[1]): narrow margins-only ranges."""
    N = rng.integers(Nlo, Nhi + 1, p)
    a = np.zeros((p, R), int)
    for i in range(p):
        s0 = share[i % 2]
        a[i, 0] = int(round(s0 * N[i])); a[i, 1] = N[i] - a[i, 0]
    n = np.zeros((T, p, R, C), int)
    base = rng.uniform(0.2, 0.8, (R, C)); base /= base.sum(1, keepdims=True)
    for t in range(T):
        for i in range(p):
            for r in range(R):
                pr = np.clip(base[r] + rng.normal(0, 0.05, C) + 0.03 * t * (r == 0) * np.eye(C)[0], 0.01, None)
                n[t, i, r] = rng.multinomial(a[i, r], pr / pr.sum())
    return a, n.sum(2), n


def true_val(coef, n):
    return float(np.asarray(coef) @ n.reshape(-1))


# =====================================================================================
# E1. Stacking and separability (Lemma E1, Proposition E3).
# =====================================================================================
ok_sep, ok_dd, ok_shear, n_shear = True, True, True, 0
for rep in range(40):
    p, R, C, T = 4, 3, 3, 3
    a, b, n = rand_margins(p, R, C, T)
    L = CountLP(a, b)
    # separable target: sum over periods of per-period targets -> range = Minkowski sum of ranges
    tg = [rng.normal(size=L.nn) * (np.repeat(np.arange(T), p * R * C) == t) for t in range(T)]
    full = L.range(sum(tg))
    parts = [L.range(v) for v in tg]
    ok_sep &= abs(full[0] - sum(x[0] for x in parts)) < 1e-7 and abs(full[1] - sum(x[1] for x in parts)) < 1e-7
    # one-column two-row DiD range = difference of per-period cellwise DD gap ranges (6a Lemma 2.1(b))
    gaps = []
    for t in range(T):
        lo = sum(max(0, a[i, 0] + b[t, i, 1] - a[i].sum()) for i in range(p)) / a[:, 0].sum() - \
             sum(min(a[i, 1], b[t, i, 1]) for i in range(p)) / a[:, 1].sum()
        hi = sum(min(a[i, 0], b[t, i, 1]) for i in range(p)) / a[:, 0].sum() - \
             sum(max(0, a[i, 1] + b[t, i, 1] - a[i].sum()) for i in range(p)) / a[:, 1].sum()
        gaps.append((lo, hi))
    d = L.range(did_vec(L, 0, 1, 1, 0, 2))
    ok_dd &= abs(d[0] - (gaps[2][0] - gaps[0][1])) < 1e-9 and abs(d[1] - (gaps[2][1] - gaps[0][0])) < 1e-9
    # the event-study pair (D_{0->1}, D_{1->2}) shares period 1: the joint set is not the box of ranges
    D01, D12 = did_vec(L, 0, 1, 1, 0, 1), did_vec(L, 0, 1, 1, 1, 2)
    r01, r12 = L.range(D01), L.range(D12)
    if r01[1] - r01[0] > 1e-6 and r12[1] - r12[0] > 1e-6:
        # corner (max D01, max D12) needs gap_1 at its max and at its min at once
        L2 = CountLP(a, b)
        L2.add_ge(D01, r01[1] - 1e-9)
        hi12 = L2.range(D12)[1]
        n_shear += 1
        ok_shear &= hi12 < r12[1] - 1e-6
check("E1 [V] Lemma E1 / Prop. E3(a): T-period ranges of separable targets are Minkowski sums of per-period ranges", ok_sep)
check("E1 [V] Prop. E3(b): two-row one-column DiD range = [G-_t1 - G+_t0, G+_t1 - G-_t0] (any R, C, T)", ok_dd)
check("E1 [F] Prop. E3(c): the joint set of event-study DiDs is not the box of their ranges (shared base table)",
      ok_shear and n_shear > 20, f"{n_shear} designs, the box corner is infeasible in all")


# =====================================================================================
# E2. Panel (within-person) tables add nothing for period-margin targets (Lemma E2), integer too.
# =====================================================================================
ok_panel = True
for rep in range(200):
    p, R, C, T = 3, 2, 3, 3
    a, b, n = rand_margins(p, R, C, T, 5, 12)
    # integer panel table with the given (i, r, t) margins: list each person once per period
    for i in range(p):
        for r in range(R):
            labels = [np.repeat(np.arange(C), n[t, i, r]) for t in range(T)]
            rng.shuffle(labels[1])
            panel = np.zeros((C,) * T, int)
            for j in range(a[i, r]):
                panel[tuple(labels[t][j] for t in range(T))] += 1
            for t in range(T):
                ax = tuple(k for k in range(T) if k != t)
                ok_panel &= np.array_equal(panel.sum(axis=ax), n[t, i, r])
check("E2 [V] Lemma E2: every family of period tables extends to an integer within-person panel table", ok_panel)


# =====================================================================================
# E5/E6/E26. Brute force of small two-period integer tables, every count-layer rung (2x2 units).
#     Units are 2x2: n[t, i, 0, 0] (treated-group yes) determines the unit table.
# =====================================================================================
def brute_2x2(a, b, extra=None):
    """Enumerate integer two-period 2x2 tables. Returns list of dicts with the target pieces."""
    T, p, _ = b.shape
    ranges = []
    for t in range(T):
        for i in range(p):
            N = a[i].sum()
            lo, hi = max(0, a[i, 0] + b[t, i, 0] - N), min(a[i, 0], b[t, i, 0])
            ranges.append(range(lo, hi + 1))
    out = []
    for combo in itertools.product(*ranges):
        y = np.array(combo).reshape(T, p)          # treated-group 'yes' counts
        out.append(y)
    return out


def pieces_2x2(a, b, y, t0=0, t1=1):
    """One-vs-rest (= two-row) DiD, W1 change associations, ER/NM on changes, for the yes column."""
    N = a.sum(1).astype(float)
    pw = N / N.sum()
    x = a[:, 0] / N
    A1, A0 = a[:, 0].sum(), a[:, 1].sum()
    Y = b[:, :, 0]
    dY1 = y[t1] - y[t0]                         # treated change counts per unit
    dYt = Y[t1] - Y[t0]
    D = dY1.sum() / A1 - (dYt.sum() - dY1.sum()) / A0
    dT = dYt / N
    mu = np.sum(pw * x)
    varx = np.sum(pw * x * x) - mu ** 2
    ER = (np.sum(pw * x * dT) - mu * np.sum(pw * dT)) / varx
    NM = (np.sum(pw * x * dT) - mu * np.sum(pw * dT)) / (mu * (1 - mu))
    dW = varx * (ER - D)
    dB = mu * (1 - mu) * (D - NM)
    return D, ER, NM, dW, dB


ok_brute = {"margins": True, "dW>=0": True, "dB>=0": True, "CR": True, "contrast": True}
n_cr_int_empty, n_designs = 0, 0
for rep in range(60):
    p, T = 3, 2
    N = rng.integers(3, 7, p)
    a = np.column_stack([rng.integers(1, N), np.zeros(p, int)])
    a[:, 1] = N - a[:, 0]
    if len(set(np.round(a[:, 0] / N, 12))) == 1:      # Var_p(X) > 0 is assumed (Prop. E5)
        continue
    b = np.zeros((T, p, 2), int)
    for t in range(T):
        b[t, :, 0] = [rng.integers(0, N[i] + 1) for i in range(p)]
        b[t, :, 1] = N - b[t, :, 0]
    tabs = brute_2x2(a, b)
    vals = [pieces_2x2(a, b, y) for y in tabs]
    Ds = np.array([v[0] for v in vals])
    _, ER, NM, _, _ = vals[0]
    A1, A0 = a[:, 0].sum(), a[:, 1].sum()
    step = 1 / A1 + 1 / A0
    Dlo, Dhi = Ds.min(), Ds.max()
    # margins only: every lattice point between the extremes is attained (Prop. E26(a))
    lattice = np.round((Ds - Dlo) / step).astype(int)
    ok_brute["margins"] &= set(lattice) == set(range(int(round((Dhi - Dlo) / step)) + 1)) and \
        np.allclose(Dlo + lattice * step, Ds)
    # continuous range by LP equals brute-force extremes (integer endpoints attained)
    L = CountLP(a, b.astype(float))
    cr = L.range(did_vec(L, 0, 1, 0, 0, 1))
    ok_brute["margins"] &= abs(cr[0] - Dlo) < 1e-9 and abs(cr[1] - Dhi) < 1e-9
    # W1 signs on changes: integer set = lattice n continuous half-line (Prop. E26(b))
    for key, cond, cont in (("dW>=0", lambda v: v[3] >= -1e-12, (Dlo, min(Dhi, ER))),
                            ("dB>=0", lambda v: v[4] >= -1e-12, (max(Dlo, NM), Dhi))):
        S = sorted({round(v[0], 12) for v in vals if cond(v)})
        expect = sorted({round(d, 12) for d in Ds if cont[0] - 1e-12 <= d <= cont[1] + 1e-12})
        ok_brute[key] &= S == expect
    # CR on changes (union): integer set = lattice points between NM and ER; may be empty although the
    # continuous set is not (Remark M2.1 analogue)
    S = sorted({round(v[0], 12) for v in vals if v[3] * v[4] >= -1e-15})
    lo_c, hi_c = (min(NM, ER), max(NM, ER))
    expect = sorted({round(d, 12) for d in Ds if lo_c - 1e-12 <= d <= hi_c + 1e-12})
    ok_brute["CR"] &= S == expect
    n_cr_int_empty += (len(S) == 0)
    n_designs += 1
    # contrast on changes between bins L = {unit 0}, H = {units 1, 2} (Prop. E6), if t > 0
    xL = a[0, 0] / a[0].sum()
    xH = a[1:, 0].sum() / a[1:].sum()
    if abs(xH - xL) > 1e-9:
        sgn = 1 if xH > xL else -1
        Lb, Hb = ([0], [1, 2]) if sgn > 0 else ([1, 2], [0])
        def binpieces(y):
            r = {}
            for nm, U in (("L", Lb), ("H", Hb)):
                A1b, A0b, Nb = a[U, 0].sum(), a[U, 1].sum(), a[U].sum()
                d1 = (y[1, U] - y[0, U]).sum() / A1b
                d0 = ((b[1, U, 0] - y[1, U]) - (b[0, U, 0] - y[0, U])).sum() / A0b
                r[nm] = (A1b / Nb, d1, d0, (b[1, U, 0] - b[0, U, 0]).sum() / Nb)
            tt = r["H"][0] - r["L"][0]
            Dbar = 0.5 * ((r["H"][1] - r["H"][2]) + (r["L"][1] - r["L"][2]))
            mix = tt * Dbar
            xb = 0.5 * (r["H"][0] + r["L"][0])
            rate = xb * (r["H"][1] - r["L"][1]) + (1 - xb) * (r["H"][2] - r["L"][2])
            dtau = r["H"][3] - r["L"][3]
            return Dbar, mix, rate, dtau, tt
        bp = [binpieces(y) for y in tabs]
        ok_brute["contrast"] &= all(abs(m + rt - dt) < 1e-12 for (_, m, rt, dt, _) in bp)   # Kitagawa
        dtau, tt = bp[0][3], bp[0][4]
        S = sorted({round(q[0], 12) for q in bp if q[1] >= -1e-12 and q[2] >= -1e-12})
        allD = sorted({round(q[0], 12) for q in bp})
        expect = sorted({d for d in allD if -1e-12 <= d <= dtau / tt + 1e-12})
        ok_brute["contrast"] &= S == expect
for k, v in ok_brute.items():
    check(f"E5/E6/E26 [V] brute force, two-period integer 2x2 tables, rung '{k}': integer set as stated (Props E3, E5, E6, E26)",
          v, f"{n_designs} designs" + (f"; integer CR empty in {n_cr_int_empty} with nonempty continuous set"
                                       if k == "CR" else ""))


# =====================================================================================
# E4/E5. Rows on changes in R x C (Lemma E4) and Delta-4 beyond 2x2 (Prop. E5).
# =====================================================================================
def design_stats(a):
    N = a.sum(1).astype(float)
    pw = N / N.sum()
    X = a / N[:, None]
    return N, pw, X


def assoc_from_beta(a, beta, Tm, r):
    """Realised W1 associations of row r for unit rates beta[i, s] and unit outcome shares Tm[i]
    (one column), from the definitions of 2c Def. 1 (any bounded outcome, e.g. a change)."""
    N, pw, X = design_stats(a)
    xr = X[:, r]
    mu = np.sum(pw * xr)
    wr = pw * xr / mu
    wnr = pw * (1 - xr) / (1 - mu)
    bnr = (Tm - xr * beta[:, r]) / np.where(1 - xr > 0, 1 - xr, 1)
    cov = lambda f, g, w: np.sum(w * f * g) - np.sum(w * f) * np.sum(w * g)
    dW = mu * cov(beta[:, r], xr, wr) + (1 - mu) * cov(bnr, xr, wnr)
    dB = np.sum(pw * xr * (beta[:, r] - Tm))
    return dW, dB


ok_id, ok_cf, ok_nm, ok_one, n_cf = True, True, True, True, 0
for rep in range(60):
    p, R, C, T = 6, 3, 3, 2
    a, b, n = rand_margins(p, R, C, T)
    N, pw, X = design_stats(a)
    r, c = 0, 1
    xr = X[:, r]
    if np.var(xr) < 1e-6:
        continue
    beta = n[:, :, :, c] / np.where(a > 0, a, 1)[None]          # beta[t, i, s]
    dbeta = beta[1] - beta[0]
    Tm = b[:, :, c] / N[None]
    dT = Tm[1] - Tm[0]
    mu = np.sum(pw * xr)
    varx = np.sum(pw * xr * xr) - mu ** 2
    ER = (np.sum(pw * xr * dT) - mu * np.sum(pw * dT)) / varx
    NM = (np.sum(pw * xr * dT) - mu * np.sum(pw * dT)) / (mu * (1 - mu))
    L = CountLP(a, b)
    Dv = did_vec(L, r, "rest", c, 0, 1)
    Dtrue = true_val(Dv, n)
    dW, dB = assoc_from_beta(a, dbeta, dT, r)
    # Lemma E4 / collapsed identities on changes (2c M1(c) for the outcome Delta Y)
    ok_id &= abs(ER - Dtrue - dW / varx) < 1e-10 and abs(NM - Dtrue + dB / (mu * (1 - mu))) < 1e-10
    # closed form (Prop. E5) vs LP over the R x C two-period product with the W1 change rows
    Dlo, Dhi = L.range(Dv)
    for sgnW, sgnB in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)):
        LL = CountLP(a, b)
        if sgnW:
            LL.add_ge(sgnW * (-Dv), sgnW * (-ER))      # sgnW*dW >= 0  <=>  sgnW*(ER - D) >= 0
        if sgnB:
            LL.add_ge(sgnB * Dv, sgnB * NM)            # sgnB*dB >= 0  <=>  sgnB*(D - NM) >= 0
        lo, hi = Dlo, Dhi
        if sgnW == 1: hi = min(hi, ER)
        if sgnW == -1: lo = max(lo, ER)
        if sgnB == 1: lo = max(lo, NM)
        if sgnB == -1: hi = min(hi, NM)
        got = LL.range(Dv)
        if lo > hi + 1e-9:
            ok_cf &= got is None
        else:
            ok_cf &= got is not None and abs(got[0] - lo) < 1e-7 and abs(got[1] - hi) < 1e-7
        n_cf += 1
    ok_nm &= Dlo - 1e-9 <= NM <= Dhi + 1e-9
    gamma = varx / (mu * (1 - mu))
    ok_one &= abs(NM - gamma * ER) < 1e-12 and (gamma < 1 - 1e-12)
check("E4 [V] Lemma E4: change associations satisfy ER_D - D = dW_D / Var(X_r), NMc_D - D = -dB_D / (mu(1-mu)) in R x C", ok_id)
check("E5 [V] Prop. E5(a)-(b): one-vs-rest closed forms on changes = LP over the R x C two-period product", ok_cf, f"{n_cf} LPs")
check("E5 [V] Prop. E5(c)-(d): NM on changes always in the margins-only range; NM_D = gamma ER_D, gamma < 1", ok_nm and ok_one)

# pairwise between sign on changes (R = 3) and joint columns: valid rows, and joint tightening
ok_valid, n_tight_pair, n_tight_joint, n_rep = True, 0, 0, 0
for rep in range(40):
    p, R, C, T = 6, 3, 3, 2
    a, b, n = rand_margins(p, R, C, T)
    N, pw, X = design_stats(a)
    L = CountLP(a, b)
    Dv = did_vec(L, 0, 1, 0, 0, 1)            # pairwise DiD, rows 0 vs 1, column 0
    Dtrue = true_val(Dv, n)
    # realised pairwise between association on changes, delta_B(0,1) = E_p[X0 X1 (Dbeta_0 - Dbeta_1)]
    v = np.zeros(L.nn)
    for i in range(p):
        for t, sg in ((1, 1), (0, -1)):
            v[L.ix(t, i, 0, 0)] += sg * pw[i] * X[i, 1] / N[i]
            v[L.ix(t, i, 1, 0)] -= sg * pw[i] * X[i, 0] / N[i]
    s0 = np.sign(true_val(v, n)) or 1.0
    L1 = CountLP(a, b); L1.add_ge(s0 * v, 0)
    base, got = L.range(Dv), L1.range(Dv)
    ok_valid &= got[0] - 1e-9 <= Dtrue <= got[1] + 1e-9
    n_tight_pair += (got[1] - got[0]) < (base[1] - base[0]) - 1e-7
    # joint W1 within-signs on changes for columns 0 and 1 (row 0), each with the truth's sign
    LLs = []
    rowsW = []
    for c in (0, 1):
        Tm = b[:, :, c] / N[None]; dT = Tm[1] - Tm[0]; xr = X[:, 0]; mu = np.sum(pw * xr)
        varx = np.sum(pw * xr * xr) - mu ** 2
        ER = (np.sum(pw * xr * dT) - mu * np.sum(pw * dT)) / varx
        Dc = did_vec(L, 0, "rest", c, 0, 1)
        sg = np.sign(ER - true_val(Dc, n)) or 1.0
        rowsW.append((sg * (-Dc), sg * (-ER)))
    tgt = did_vec(L, 0, "rest", 0, 0, 1)
    Ls = CountLP(a, b); Ls.add_ge(*rowsW[0])
    Lj = CountLP(a, b); Lj.add_ge(*rowsW[0]); Lj.add_ge(*rowsW[1])
    rs, rj = Ls.range(tgt), Lj.range(tgt)
    tv = true_val(tgt, n)
    ok_valid &= rj[0] - 1e-9 <= tv <= rj[1] + 1e-9
    n_tight_joint += (rj[1] - rj[0]) < (rs[1] - rs[0]) - 1e-7
    n_rep += 1
check("E5 [V] Prop. E5(e): pairwise and joint-column change rows are valid LP rows (truth inside); they can tighten",
      ok_valid and n_tight_pair > 0 and n_tight_joint > 0,
      f"pairwise sign tightened {n_tight_pair}/{n_rep}; second column's sign tightened {n_tight_joint}/{n_rep}")

# T = 3: CR on changes in two period pairs = union of four branches; truth inside when CR holds in both
ok_union, n_u = True, 0
for rep in range(400):
    p, R, C, T = 5, 2, 2, 3
    a, b, n = rand_margins(p, R, C, T)
    N, pw, X = design_stats(a)
    xr = X[:, 0]; mu = np.sum(pw * xr); varx = np.sum(pw * xr * xr) - mu ** 2
    if varx < 1e-6:
        continue
    L = CountLP(a, b)
    pairs = [(0, 1), (1, 2)]
    info = []
    for (t0, t1) in pairs:
        dT = (b[t1, :, 0] - b[t0, :, 0]) / N
        ER = (np.sum(pw * xr * dT) - mu * np.sum(pw * dT)) / varx
        NM = (np.sum(pw * xr * dT) - mu * np.sum(pw * dT)) / (mu * (1 - mu))
        Dv = did_vec(L, 0, 1, 0, t0, t1)
        info.append((Dv, ER, NM, true_val(Dv, n)))
    holds = all((ER - Dt) * (Dt - NM) >= -1e-12 for (_, ER, NM, Dt) in info)
    if not holds:
        continue
    tgt = did_vec(L, 0, 1, 0, 0, 2)           # long difference, also covered
    union = []
    for sgns in itertools.product((1, -1), repeat=2):
        LL = CountLP(a, b)
        for sg, (Dv, ER, NM, _) in zip(sgns, info):
            LL.add_ge(sg * (-Dv), sg * (-ER)); LL.add_ge(sg * Dv, sg * NM)
        rg = LL.range(tgt)
        if rg is not None:
            union.append(rg)
    tv = true_val(tgt, n)
    ok_union &= any(lo - 1e-9 <= tv <= hi + 1e-9 for lo, hi in union)
    n_u += 1
check("E5 [V] Prop. E5(f): T = 3, CR on changes in two pairs: the union of the 2^2 branch LPs contains the truth", ok_union and n_u > 5, f"{n_u} designs with realised CR in both pairs")


# =====================================================================================
# E6. Finite contrasts on changes in the count layer, R x C rows (Prop. E6), LP = closed form.
# =====================================================================================
ok_c6, n_c6, n_empty6 = True, 0, 0
for rep in range(60):
    p, R, C, T = 6, 2, 3, 2
    a, b, n = rand_margins(p, R, C, T)
    N, pw, X = design_stats(a)
    order = np.argsort(X[:, 0])
    Lb, Hb = list(order[:3]), list(order[3:])
    xL, xH = a[Lb, 0].sum() / a[Lb].sum(), a[Hb, 0].sum() / a[Hb].sum()
    if xH - xL < 1e-6 or min(a[Lb, 0].sum(), a[Lb, 1].sum(), a[Hb, 0].sum(), a[Hb, 1].sum()) == 0:
        continue
    L = CountLP(a, b)
    c = 0
    def bin_gap(U):
        return (rate_vec(L, 0, c, 1, U) - rate_vec(L, 0, c, 0, U)) - (rate_vec(L, 1, c, 1, U) - rate_vec(L, 1, c, 0, U))
    Dbar = 0.5 * (bin_gap(Hb) + bin_gap(Lb))
    tt = xH - xL
    dtau = (b[1, Hb, c].sum() - b[0, Hb, c].sum()) / a[Hb].sum() - (b[1, Lb, c].sum() - b[0, Lb, c].sum()) / a[Lb].sum()
    LL = CountLP(a, b)
    LL.add_ge(tt * Dbar, 0)                    # mix term >= 0
    LL.add_le(tt * Dbar, dtau)                 # rate term = dtau - mix >= 0
    got = LL.range(Dbar)
    base = L.range(Dbar)
    lo, hi = max(base[0], 0), min(base[1], dtau / tt)
    if lo > hi + 1e-9:
        ok_c6 &= got is None
        n_empty6 += 1
    else:
        ok_c6 &= got is not None and abs(got[0] - lo) < 1e-7 and abs(got[1] - hi) < 1e-7
    n_c6 += 1
check("E6 [V] Prop. E6: count-layer contrast on changes = [0, (dtau_H - dtau_L)/t] n margins-only range (LP, 2x3)", ok_c6,
      f"{n_c6} designs; empty in {n_empty6} (the joint sign test rejects)")


# =====================================================================================
# E7. FD in R x C (Lemma E7): Delta Q identified from E[Delta T | X]; Delta K not; the DiD for a
#     column depends on Delta K only through that column's Delta K_c.
# =====================================================================================
def monomials(Xm):
    R = Xm.shape[1]
    cols, idx = [], []
    for r in range(R):
        for s in range(r, R):
            cols.append(Xm[:, r] * Xm[:, s] * (1 if r == s else 2)); idx.append((r, s))
    return np.column_stack(cols), idx


ok7 = True
for rep in range(30):
    p, R, C = 60, 3, 2
    Xm = rng.dirichlet(np.ones(R), p)
    w = rng.integers(50, 500, p).astype(float); pw = w / w.sum()
    lev = 0.3 + 0.4 * np.sin(3 * Xm @ rng.normal(size=(R, R)) + rng.normal(size=R)) ** 2   # non-linear levels, col 0
    dH = rng.uniform(-0.05, 0.05, (R, R))                                                  # column 0; column 1 = -dH
    dbeta0 = Xm @ dH.T                                     # M_dH(x)_r = (dH x)_r
    dT = np.einsum("ir,ir->i", Xm, dbeta0)
    Z, idx = monomials(Xm)
    q = np.linalg.lstsq(Z * np.sqrt(pw)[:, None], dT * np.sqrt(pw), rcond=None)[0]
    dQ = (dH + dH.T) / 2
    ok7 &= max(abs(q[k] - dQ[r, s]) for k, (r, s) in enumerate(idx)) < 1e-10
    # add an antisymmetric Delta K' to column 0 (and -Delta K' to column 1): same Delta T, different DiD
    Kp = rng.uniform(-0.03, 0.03, (R, R)); Kp = Kp - Kp.T
    dbeta0b = Xm @ (dH + Kp).T
    ok7 &= np.allclose(np.einsum("ir,ir->i", Xm, dbeta0b), dT)
    wr = lambda r: pw * Xm[:, r] / np.sum(pw * Xm[:, r])
    did = lambda db: np.sum(wr(0) * db[:, 0]) - np.sum(wr(1) * db[:, 1])
    ok7 &= abs(did(dbeta0b) - did(dbeta0)) > 1e-6
    # the change in the DiD is linear in K'_c with design coefficients: (w_0 . X K'^T)_0 - (w_1 . X K'^T)_1
    pred = np.sum(wr(0) * (Xm @ Kp.T)[:, 0]) - np.sum(wr(1) * (Xm @ Kp.T)[:, 1])
    ok7 &= abs(did(dbeta0b) - did(dbeta0) - pred) < 1e-12
check("E7 [V] Lemma E7: FD in R x C, Delta Q identified exactly with non-linear levels; Delta K_c moves the DiD, not E[Delta T|X]", ok7)


# =====================================================================================
# E8. Sharp identified set under FD (Prop. E8) versus the outer set FD4 (2d), 2x2, discrete design.
#     Sharp: Dm_v(x_j) in [E L(T1) - E U(T0), E U(T1) - E L(T0)] given x_j; outer: the same at tau.
# =====================================================================================
Lf = lambda t, x: np.maximum(0, (t - (1 - x)) / x)
Uf = lambda t, x: np.minimum(1, t / x)
ok8_in, ok8_sub, n8_strict, n8 = True, True, 0, 0
for rep in range(100):
    xs = np.array([0.2, 0.5, 0.8]); m = 30
    x = np.repeat(xs, m); p = len(x); pw = np.full(p, 1 / p)
    dQvv, dQvn, dQnn, dk = rng.uniform(-0.04, 0.04, 4)
    dmv = lambda xx, k: dQvv * xx + (dQvn + k) * (1 - xx)
    dmn = lambda xx, k: (dQvn - k) * xx + dQnn * (1 - xx)
    bv0, bn0 = rng.uniform(0.1, 0.9, p), rng.uniform(0.1, 0.9, p)
    ev, en = rng.uniform(-0.05, 0.05, p), rng.uniform(-0.05, 0.05, p)
    for j in range(3):                                   # FD holds exactly at the design
        sl = slice(j * m, (j + 1) * m); ev[sl] -= ev[sl].mean(); en[sl] -= en[sl].mean()
    bv1, bn1 = bv0 + dmv(x, dk) + ev, bn0 + dmn(x, dk) + en
    if (np.minimum(bv1, bn1) < 0).any() or (np.maximum(bv1, bn1) > 1).any():
        continue
    T0, T1 = x * bv0 + (1 - x) * bn0, x * bv1 + (1 - x) * bn1
    lo_s, hi_s, lo_o, hi_o = -np.inf, np.inf, -np.inf, np.inf
    for j, xj in enumerate(xs):
        sl = slice(j * m, (j + 1) * m)
        # Dm_v(xj, k) = const + (1 - xj) k ; interval for Dm_v -> interval for k
        const = dQvv * xj + dQvn * (1 - xj)
        sh = (np.mean(Lf(T1[sl], xj)) - np.mean(Uf(T0[sl], xj)), np.mean(Uf(T1[sl], xj)) - np.mean(Lf(T0[sl], xj)))
        t0m, t1m = T0[sl].mean(), T1[sl].mean()
        ou = (Lf(t1m, xj) - Uf(t0m, xj), Uf(t1m, xj) - Lf(t0m, xj))
        lo_s = max(lo_s, (sh[0] - const) / (1 - xj)); hi_s = min(hi_s, (sh[1] - const) / (1 - xj))
        lo_o = max(lo_o, (ou[0] - const) / (1 - xj)); hi_o = min(hi_o, (ou[1] - const) / (1 - xj))
    n8 += 1
    ok8_in &= lo_s - 1e-12 <= dk <= hi_s + 1e-12
    ok8_sub &= lo_o - 1e-12 <= lo_s and hi_s <= hi_o + 1e-12
    n8_strict += (hi_s - lo_s) < (hi_o - lo_o) - 1e-9
check("E8 [V] Prop. E8: the sharp FD set contains the truth and lies inside the FD4 outer set; strictly when T varies within x",
      ok8_in and ok8_sub and n8_strict > 0, f"{n8} designs; strictly narrower in {n8_strict}")
# noise-free (T_t degenerate given x): sharp = outer
ok8_eq = True
for rep in range(50):
    xs = np.array([0.2, 0.5, 0.8])
    dQvv, dQvn, dQnn = rng.uniform(-0.04, 0.04, 3)
    t0, t1 = rng.uniform(0.1, 0.9, 3), None
    t1 = np.clip(t0 + rng.uniform(-0.05, 0.05, 3), 0.01, 0.99)
    for j, xj in enumerate(xs):
        sh = (Lf(t1[j], xj) - Uf(t0[j], xj), Uf(t1[j], xj) - Lf(t0[j], xj))
        ok8_eq &= np.allclose(sh, (np.mean(Lf(np.array([t1[j]]), xj)) - np.mean(Uf(np.array([t0[j]]), xj)),
                                   np.mean(Uf(np.array([t1[j]]), xj)) - np.mean(Lf(np.array([t0[j]]), xj))))
check("E8 [I] Prop. E8(c): with T_t degenerate given X the sharp set equals the FD4 outer set", ok8_eq)


# =====================================================================================
# E9. DiD-S within FD (Prop. E9), 2x2: rows O_D >= 0 and L_D >= 0 on the hull <=> |Dk| <= min_hull Dq;
#     sharp set = sharp FD set n that interval (closed form sharp iff the interval lies in the FD set).
# =====================================================================================
ok9_cf, ok9_valid, n9_eq, n9_strict, n9 = True, True, 0, 0, 0
for rep in range(300):
    xs = np.array([0.2, 0.5, 0.8]); m = 20
    x = np.repeat(xs, m); p = len(x)
    dQvv, dQvn, dQnn = rng.uniform(-0.06, 0.06, 3)
    dq = lambda xx: dQvv * xx + dQvn * (1 - xx) - dQvn * xx - dQnn * (1 - xx)
    qmin = min(dq(xs[0]), dq(xs[-1]))
    if qmin < 0:
        continue                                          # DiD-S refuted by the data (test)
    dk = rng.uniform(-qmin, qmin)                         # truth satisfies DiD-S
    # rows by brute force on a fine grid of k: O = dq + k >= 0 and L = dq - k >= 0 on the hull
    grid = np.linspace(-0.2, 0.2, 40001)
    okg = np.all([(dq(xx) + grid >= -1e-12) & (dq(xx) - grid >= -1e-12) for xx in np.linspace(xs[0], xs[-1], 11)], axis=0)
    ok9_cf &= abs(grid[okg].min() + qmin) < 2e-5 and abs(grid[okg].max() - qmin) < 2e-5
    dmv = lambda xx, k: dQvv * xx + (dQvn + k) * (1 - xx)
    dmn = lambda xx, k: (dQvn - k) * xx + dQnn * (1 - xx)
    # levels concentrated near 0 and 1 (so that the FD bounds can bind), valid in both periods
    Dv_, Dn_ = dmv(x, dk), dmn(x, dk)
    bv0 = np.maximum(0, -Dv_) + rng.beta(0.2, 0.2, p) * (1 - np.abs(Dv_))
    bn0 = np.maximum(0, -Dn_) + rng.beta(0.2, 0.2, p) * (1 - np.abs(Dn_))
    if rep % 2 == 0:              # every unit at the middle composition starts at 0 (a floor that binds)
        sl = slice(m, 2 * m)
        if Dv_[sl].min() < 0 or Dn_[sl].min() < 0:
            continue
        bv0[sl], bn0[sl] = 0.0, 0.0
    bv1, bn1 = bv0 + Dv_, bn0 + Dn_
    T0, T1 = x * bv0 + (1 - x) * bn0, x * bv1 + (1 - x) * bn1
    lo_s, hi_s = -np.inf, np.inf
    for j, xj in enumerate(xs):
        sl = slice(j * m, (j + 1) * m)
        const = dQvv * xj + dQvn * (1 - xj)
        sh = (np.mean(Lf(T1[sl], xj)) - np.mean(Uf(T0[sl], xj)), np.mean(Uf(T1[sl], xj)) - np.mean(Lf(T0[sl], xj)))
        lo_s = max(lo_s, (sh[0] - const) / (1 - xj)); hi_s = min(hi_s, (sh[1] - const) / (1 - xj))
    lo, hi = max(lo_s, -qmin), min(hi_s, qmin)
    ok9_valid &= lo - 1e-12 <= dk <= hi + 1e-12
    n9 += 1
    if lo_s <= -qmin + 1e-12 and hi_s >= qmin - 1e-12:
        n9_eq += 1
    else:
        n9_strict += 1
# a constructed instance: dQ = (.03, 0, -.03), so Dq = .03 on the hull; at x = .5 every unit starts at 0,
# which forces Dm_v(.5) >= 0 and Dm_n(.5) >= 0, i.e. Delta kappa = -.03: sharp set {-.03} vs closed form [-.03, .03]
xj = 0.5; c_v, c_n = 0.03 * xj, -0.03 * (1 - xj)
k_lo, k_hi = -c_v / (1 - xj), -c_n / xj       # Dm_v = c_v + (1-x)k >= 0 ; Dm_n = c_n - x k >= 0
n9_strict += (abs(k_lo + 0.03) < 1e-12 and abs(k_hi + 0.03) < 1e-12)
check("E9 [V] Prop. E9(a): 2x2 DiD-S rows on the hull <=> |Delta kappa| <= min(Dq(x_L), Dq(x_U)) (closed form)", ok9_cf)
check("E9 [V] Prop. E9(b)-(c): DiD-S sharp set = sharp FD set n closed form; contains the truth; closed form sharp only sometimes",
      ok9_valid and n9_eq > 0 and n9_strict > 0, f"{n9} designs: closed form sharp in {n9_eq}, strictly wider in {n9_strict}")


# =====================================================================================
# E10. Contrasts on changes in the LCE layer (Prop. E10): with Dm linear, at compositions x_L, x_H with
#      x_H - x_L = t (e_r - e_s): mix = t O_D(xbar), rate = t L_D(xbar) (P-S objects at the midpoint).
# =====================================================================================
ok10 = True
for rep in range(200):
    R = 3
    dH = rng.normal(size=(R, R))
    xL = rng.dirichlet(np.ones(R)); r, s = 0, 2; t = rng.uniform(0, min(xL[s], 1))
    xH = xL + t * (np.eye(R)[r] - np.eye(R)[s])
    m = lambda xx: dH @ xx
    mix = (xH - xL) @ (m(xH) + m(xL)) / 2
    rate = ((xH + xL) / 2) @ (m(xH) - m(xL))
    xb = (xH + xL) / 2
    O = (dH @ xb)[r] - (dH @ xb)[s]
    Lc = (xb @ dH)[r] - (xb @ dH)[s]
    ok10 &= abs(mix - t * O) < 1e-12 and abs(rate - t * Lc) < 1e-12 and abs(mix + rate - (xH @ m(xH) - xL @ m(xL))) < 1e-12
check("E10 [V] Prop. E10: LCE-layer contrast on changes: mix = t O_D(xbar), rate = t L_D(xbar), sum = Dtau_H - Dtau_L", ok10)


# =====================================================================================
# E11. Per-period LCE with a model sign on changes (Prop. E11): delta*_W on changes is affine in
#      (kappa_0, kappa_1), here = identified - (a + b) Delta kappa (2x2 M0); sharp set = product n row.
# =====================================================================================
ok11, n11_tight = True, 0
for rep in range(100):
    p = 50
    x = rng.uniform(0.05, 0.9, p); w = rng.integers(50, 500, p).astype(float); pw = w / w.sum()
    mu = np.sum(pw * x)
    a_ = mu * (np.sum(pw * x * x * x) / mu - (np.sum(pw * x * x) / mu) ** 2)
    wn = pw * (1 - x); b_ = (1 - mu) * (np.sum(wn * x * x) / wn.sum() - (np.sum(wn * x) / wn.sum()) ** 2)
    def deltaW(Hm):     # Hm = [[Hvv, Hvn], [Hnv, Hnn]]; m_v = Hvv x + Hvn (1-x), m_n = Hnv x + Hnn (1-x)
        mv = Hm[0, 0] * x + Hm[0, 1] * (1 - x); mn = Hm[1, 0] * x + Hm[1, 1] * (1 - x)
        wv_ = pw * x / mu; wn_ = pw * (1 - x) / (1 - mu)
        cov = lambda f, g, ww: np.sum(ww * f * g) - np.sum(ww * f) * np.sum(ww * g)
        return mu * cov(mv, x, wv_) + (1 - mu) * cov(mn, x, wn_)
    Q0 = rng.uniform(0.2, 0.8, 3); Q1 = rng.uniform(0.2, 0.8, 3); k0, k1 = rng.uniform(-0.1, 0.1, 2)
    H = lambda Q, k: np.array([[Q[0], Q[1] + k], [Q[1] - k, Q[2]]])
    dW = deltaW(H(Q1, k1)) - deltaW(H(Q0, k0))
    ident = (a_ * ((Q1[0] - Q1[1]) - (Q0[0] - Q0[1])) + b_ * ((Q1[1] - Q1[2]) - (Q0[1] - Q0[2])))
    ok11 &= abs(dW - (ident - (a_ + b_) * (k1 - k0))) < 1e-12
check("E11 [V] Prop. E11: under per-period M0 (2x2) the model within-association on changes = identified - (a+b) Delta kappa", ok11)


# =====================================================================================
# E12. Two-period residual allowance (Lemma E12): Var of the actual-minus-model DiD <= S_D^2 with
#      rho_i = 2 (w_ir + w_is); attained by opposite +-1 changes.
# =====================================================================================
ok12, worst12 = True, 0.0
for rep in range(300):
    wr_, ws_ = rng.uniform(0, 1, 2)
    # random discrete joint law of (b_r0, b_r1, b_s0, b_s1) in [0,1]^4 with K atoms
    K = 6
    atoms = rng.uniform(0, 1, (K, 4)); pr = rng.dirichlet(np.ones(K) * 0.3)
    Z = wr_ * (atoms[:, 1] - atoms[:, 0]) - ws_ * (atoms[:, 3] - atoms[:, 2])
    var = np.sum(pr * Z ** 2) - np.sum(pr * Z) ** 2
    bound = (2 * (wr_ + ws_)) ** 2 / 4
    ok12 &= var <= bound + 1e-12
    worst12 = max(worst12, var / bound)
# extreme law: (b_r0, b_r1, b_s0, b_s1) = (0,1,1,0) or (1,0,0,1) w.p. 1/2 each
Zx = np.array([1 * 1 - 1 * (-1), 1 * (-1) - 1 * 1]); varx_ = np.mean(Zx ** 2) - np.mean(Zx) ** 2
ok12 &= abs(varx_ - (2 * 2) ** 2 / 4) < 1e-12
check("E12 [V] Lemma E12: actual-minus-model DiD, per-unit variance <= (2(w_r+w_s))^2/4; attained by opposite changes", ok12,
      f"max ratio over random laws {worst12:.2f}; extreme law 1.00")


# =====================================================================================
# E13. Cross-layer rule on the two-period product (Theorem E13): FD + R1 (model DiD point-identified),
#      actual-target interval = OLS interval + z S_D; count layer = margins-only range of the actual DiD.
#      Coverage of the intersection = coverage of the model interval; disjointness rate <= alpha.
# =====================================================================================
from scipy.stats import norm, chi2
reps, cov_m, cov_i, rej = 300, 0, 0, 0
p = 400
x = rng.uniform(0.05, 0.6, p); N = rng.integers(200, 1500, p).astype(float); pw = N / N.sum()
mv_, mn_ = np.sum(pw * x), np.sum(pw * (1 - x))
wv, wn = pw * x / mv_, pw * (1 - x) / mn_
dQ = np.array([0.06, 0.02, 0.0]); dk_R1 = dQ[1] - dQ[2]           # R1: Delta kappa = dQvn - dQnn
lev_v, lev_n = 0.3 + 0.3 * np.sin(6 * x) ** 2, 0.25 + 0.2 * np.cos(4 * x) ** 2
Dmv = dQ[0] * x + (dQ[1] + dk_R1) * (1 - x); Dmn = (dQ[1] - dk_R1) * x + dQ[2] * (1 - x)
Zq = np.column_stack([x * x, 2 * x * (1 - x), (1 - x) ** 2])
cvec = np.array([np.sum(wv * x), np.sum(wv * (1 - x)), 0]) - np.array([0, np.sum(wn * x), np.sum(wn * (1 - x))])
# model DiD under R1 as a linear functional of dQ: DiD = sum wv Dm_v - sum wn Dm_n with dk = dQvn - dQnn
def did_R1(q):
    k = q[1] - q[2]
    return np.sum(wv * (q[0] * x + (q[1] + k) * (1 - x))) - np.sum(wn * ((q[1] - k) * x + q[2] * (1 - x)))
G = np.array([did_R1(e) for e in np.eye(3)])                       # gradient (linear)
S_D = np.sqrt(np.sum((wv * 1 + wn * 1) ** 2 * 0 + (pw * x / mv_ + pw * (1 - x) / mn_) ** 2))
alpha = 0.05
for rep in range(reps):
    ev = rng.uniform(-0.12, 0.12, p); en = rng.uniform(-0.12, 0.12, p)       # unit residuals of changes
    bv0, bn0 = lev_v, lev_n
    bv1, bn1 = bv0 + Dmv + ev, bn0 + Dmn + en
    T0, T1 = x * bv0 + (1 - x) * bn0, x * bv1 + (1 - x) * bn1
    actual = np.sum(wv * (bv1 - bv0)) - np.sum(wn * (bn1 - bn0))
    model = did_R1(dQ)
    W = np.sqrt(pw)
    qh, *_ = np.linalg.lstsq(Zq * W[:, None], (T1 - T0) * W, rcond=None)
    res = (T1 - T0) - Zq @ qh
    Bread = np.linalg.inv((Zq * pw[:, None]).T @ Zq)
    Meat = (Zq * (pw * res)[:, None]).T @ (Zq * (pw * res)[:, None])
    V = Bread @ Meat @ Bread
    est, se = G @ qh, np.sqrt(G @ V @ G)
    ci = (est - norm.ppf(1 - alpha / 4) * se - norm.ppf(1 - alpha / 4) * S_D,
          est + norm.ppf(1 - alpha / 4) * se + norm.ppf(1 - alpha / 4) * S_D)
    # count layer: margins-only range of the actual DiD (two rows, one column; Prop. E3(b), closed form)
    gap = []
    for T in (T0, T1):
        lo = np.sum(pw * x * Lf(T, x)) / mv_ - np.sum(pw * (1 - x) * np.minimum(1, T / (1 - x))) / mn_
        hi = np.sum(pw * x * Uf(T, x)) / mv_ - np.sum(pw * (1 - x) * np.maximum(0, (T - x) / (1 - x))) / mn_
        gap.append((lo, hi))
    cl = (gap[1][0] - gap[0][1], gap[1][1] - gap[0][0])
    inter = (max(ci[0], cl[0]), min(ci[1], cl[1]))
    cov_m += ci[0] <= actual <= ci[1]
    cov_i += inter[0] <= actual <= inter[1]
    rej += inter[0] > inter[1]
check("E13 [V] Theorem E13 (R27 on the product): coverage of (count set n actual-target FD interval) = model interval's; rejection <= alpha",
      cov_i == cov_m and cov_m / reps >= 0.95 - 0.02 and rej / reps <= alpha,
      f"{reps} replications: model interval {cov_m / reps:.3f}, intersection {cov_i / reps:.3f}, disjoint {rej / reps:.3f}")


# =====================================================================================
# Causal layer helpers. Treated row 0, comparison row 1 (or "rest"). Extra variables: theta
# (treated group's untreated rate at the post period), S (spillover onto the comparison group's rate).
# ATT = B_{0,t1}(n) - theta ;  b = D(n) + S - ATT  (counterfactual trend gap, Lemma E14).
# =====================================================================================
def causal_lp(a, b, t0, t1, c=0, comp=1, Lam=0.0, int_rows=None):
    L = CountLP(a, b)
    if int_rows:
        for coef, rhs in int_rows:
            L.add_ge(coef, rhs)
    th = L.add_var(0.0, 1.0)
    S = L.add_var(-Lam, Lam)
    Bc1 = rate_vec(L, comp, c, t1) if comp != "rest" else rest_vec(L, 0, c, t1)
    # comparison group's untreated rate B_c1 - S in [0, 1]
    e = lambda j: np.eye(L.nvar())[j]
    L.add_ge(L.pad(Bc1) - e(S), 0.0)
    L.add_le(L.pad(Bc1) - e(S), 1.0)
    att = L.pad(rate_vec(L, 0, c, t1)) - e(th)
    D = L.pad(did_vec(L, 0, comp, c, t0, t1))
    bvec = D + e(S) - att                       # b = D + S - ATT
    return L, att, bvec, e, th, S


def union_range(ranges):
    ranges = [r for r in ranges if r is not None]
    if not ranges:
        return None
    return min(r[0] for r in ranges), max(r[1] for r in ranges)


# E14. Lemma E14: ATT = D - b + S with potential outcomes (any R, C), under no anticipation.
ok14 = True
for rep in range(200):
    A1, A0 = rng.integers(20, 200, 2)
    y = {k: rng.uniform(0, 1) for k in ["r_t0", "r_t1_obs", "r_t1_0", "s_t0", "s_t1_0"]}
    S_ = rng.uniform(-0.1, 0.1); s_t1_obs = y["s_t1_0"] + S_
    D_ = (y["r_t1_obs"] - y["r_t0"]) - (s_t1_obs - y["s_t0"])
    b_ = (y["r_t1_0"] - y["r_t0"]) - (y["s_t1_0"] - y["s_t0"])
    att_ = y["r_t1_obs"] - y["r_t1_0"]
    ok14 &= abs(att_ - (D_ - b_ + S_)) < 1e-14
check("E14 [V] Lemma E14: ATT_t = D_t - b_t + S_t (no anticipation; spillover S onto the comparison group)", ok14)

# E15. Prop. E15: with no trend assumption the ATT set is [B^L_{r,t1} - 1, B^U_{r,t1}] (sharp, LP).
ok15 = True
for rep in range(30):
    a, b, n = rand_margins(4, 2, 2, 2)
    L, att, bvec, e, th, S = causal_lp(a, b, 0, 1)
    got = L.range(att)
    Lr = CountLP(a, b); Br = Lr.range(rate_vec(Lr, 0, 0, 1))
    ok15 &= abs(got[0] - (Br[0] - 1)) < 1e-8 and abs(got[1] - Br[1]) < 1e-8
check("E15 [V] Prop. E15: no trend assumption => ATT set = [B^L_r1 - 1, B^U_r1]: the descriptive DiD alone says nothing", ok15)

# E16/E17. Trend dial (Cor. E17): LP with logical constraints; equals D-set (+) [-(M+Lam), M+Lam] when the
#          logical constraints are slack, strictly narrower when they bind; NS worked example arithmetic.
ok16, n16_bind, n16_eq = True, 0, 0
for rep in range(40):
    a, b, n = rand_polar(4, 2) if rep % 2 == 0 else rand_margins(4, 2, 2, 2)
    Lr = CountLP(a, b); Dr = Lr.range(did_vec(Lr, 0, 1, 0, 0, 1))
    for M, Lam in ((0.05, 0.0), (0.05, 0.03), (0.6, 0.0)):
        L, att, bvec, e, th, S = causal_lp(a, b, 0, 1, Lam=Lam)
        L.add_le(bvec, M); L.add_ge(bvec, -M)
        got = L.range(att)
        mink = (Dr[0] - M - Lam, Dr[1] + M + Lam)
        ok16 &= got[0] >= mink[0] - 1e-8 and got[1] <= mink[1] + 1e-8
        if abs(got[0] - mink[0]) < 1e-8 and abs(got[1] - mink[1]) < 1e-8:
            n16_eq += 1
        else:
            n16_bind += 1
check("E16/E17 [V] Theorem E16 / Cor. E17: trend-dial ATT set (joint LP) lies in D-set (+) [-(M+Lam), M+Lam]; equal unless the rate limits bind",
      ok16 and n16_eq > 0 and n16_bind > 0, f"equal in {n16_eq}, logical limits bind in {n16_bind}")
ws = (0.05, 0.25)
check("E17 [I] Cor. E17: NS worked example: D in [.05,.25]; M=.10 -> [-.05,.35]; M=.02 -> [.03,.27]; breakdown M* = .05",
      np.allclose((ws[0] - .10, ws[1] + .10), (-.05, .35)) and np.allclose((ws[0] - .02, ws[1] + .02), (.03, .27)))


# =====================================================================================
# E18/E19. Relative magnitudes with set-valued pre-trends (Props E18, E19); N21.
# =====================================================================================
def rm_range(a, b, Mbar, base, post, pre_pairs, Lam=0.0, logical=True, comp=1, extra=None):
    """RM (consecutive-difference form): |b_post| <= Mbar * max over pre pairs (s-1, s) of |D_{s-1 -> s}(n)|.
    One post period: b = delta_post - delta_base = counterfactual trend gap base -> post. Union over the
    argmax pair and its sign. Without 'logical', theta and S are not constrained to rates."""
    outs = []
    for (s0, s1) in pre_pairs:
        for sg in (1, -1):
            L, att, bvec, e, th, S = causal_lp(a, b, base, post, Lam=Lam, comp=comp, int_rows=extra)
            if not logical:
                L.extra_bounds[0] = (None, None)
                L.Aub = L.Aub[:0]; L.bub = L.bub[:0]
                if extra:
                    for coef, rhs in extra:
                        L.add_ge(coef, rhs)
            Dpre = L.pad(did_vec(L, 0, comp, 0, s0, s1))
            L.add_le(bvec - Mbar * sg * Dpre, 0.0)
            L.add_ge(bvec + Mbar * sg * Dpre, 0.0)
            outs.append(L.range(att))
    return union_range(outs)


# N21 (NS 3.2): two precincts (100 men; 20 and 60 veterans), elections t = -1, 0, 1 -> code t = 0, 1, 2.
a21 = np.array([[20, 80], [60, 40]])
Yv21 = {0: (30, 45), 1: (32, 48), 2: (30, 60)}
b21 = np.array([[[Yv21[t][i], 100 - Yv21[t][i]] for i in range(2)] for t in range(3)])
L21 = CountLP(a21, b21)
D0r = L21.range(did_vec(L21, 0, 1, 0, 0, 1)); D1r = L21.range(did_vec(L21, 0, 1, 0, 1, 2))
joint_nl = rm_range(a21, b21, 1.0, 1, 2, [(0, 1)], logical=False)
A_sup = max(abs(D0r[0]), abs(D0r[1]))
ew = (D1r[0] - A_sup, D1r[1] + A_sup)
joint_l = rm_range(a21, b21, 1.0, 1, 2, [(0, 1)], logical=True)
L15, att15, *_ = causal_lp(a21, b21, 1, 2)
noinfo = L15.range(att15)
ok21 = (np.allclose(np.round(D0r, 2), (-1.23, 1.27)) and np.allclose(np.round(D1r, 2), (-1.08, 1.42))
        and np.allclose(np.round(joint_nl, 2), (-2.35, 2.65)) and np.allclose(np.round(ew, 2), (-2.35, 2.69)))
check("E18 [V] N21 reproduced exactly by LP: D_pre [-1.23,1.27], D_post [-1.08,1.42]; RM(1) joint [-2.35,2.65], endpoint-wise [-2.35,2.69]",
      ok21, f"joint {joint_nl[0]:.4f},{joint_nl[1]:.4f}; endpoint-wise {ew[0]:.4f},{ew[1]:.4f}")
check("E18 [F] N21 corrected: with theta in [0,1] the RM set is [B^L-1, B^U] = the no-assumption set, not [-2.35,2.65]",
      np.allclose(joint_l, noinfo) and noinfo[0] > -1 and noinfo[1] <= 1,
      f"RM(1) with rate limits [{joint_l[0]:.3f},{joint_l[1]:.3f}] = no assumption [{noinfo[0]:.3f},{noinfo[1]:.3f}]")
# N21 narrow variant: pre-period veteran totals pinned to +/- 2 votes (as NS), then RM informative
def pin_rows(L, t, centre, half):
    v = np.zeros(L.nn)
    for i in range(L.p):
        v[L.ix(t, i, 0, 0)] = 1
    return [(v, centre - half), (-v, -(centre + half))]
Lp = CountLP(a21, b21)
box = lambda t: (sum(max(0, a21[i, 0] + b21[t, i, 0] - 100) for i in range(2)), sum(min(a21[i, 0], b21[t, i, 0]) for i in range(2)))
pins = []
for t in (0, 1):
    lo_, hi_ = box(t); pins += pin_rows(Lp, t, lo_ + (hi_ - lo_) / 2, 2)
for rr in pins:
    Lp.add_ge(*rr)
nD0 = Lp.range(did_vec(Lp, 0, 1, 0, 0, 1))
narrow_nl = rm_range(a21, b21, 1.0, 1, 2, [(0, 1)], logical=False, extra=pins)
narrow_l = rm_range(a21, b21, 1.0, 1, 2, [(0, 1)], logical=True, extra=pins)
check("E18 [I] N21 narrow pre-trends (pre totals pinned to +/-2): D_pre [-1/16, 5/48]; RM(1) [-29/48, 43/48], identical with and without rate limits",
      np.allclose(nD0, (-1 / 16, 5 / 48)) and np.allclose(narrow_nl, (-29 / 48, 43 / 48)) and np.allclose(narrow_l, narrow_nl),
      f"with rate limits [{narrow_l[0]:.3f},{narrow_l[1]:.3f}]")
# RM union-of-LPs equals a grid brute force (N21 coordinates s_t = veteran yes totals)
cf21 = 1 / 80 + 1 / 120
gfun = lambda t, s: s * cf21 - sum(Yv21[t]) / 120
S_ = [np.linspace(*box(t), 121) for t in range(3)]
g0, g1, g2 = np.meshgrid(*S_, indexing="ij")
Dp, Dq_ = gfun(1, g1) - gfun(0, g0), gfun(2, g2) - gfun(1, g1)
grid_rm = ((Dq_ - np.abs(Dp)).min(), (Dq_ + np.abs(Dp)).max())
check("E18 [V] Prop. E18(b): RM set = union of 2 x (pre pairs) LPs; equals grid brute force on N21",
      np.allclose(grid_rm, joint_nl, atol=1e-9), f"grid {grid_rm[0]:.4f},{grid_rm[1]:.4f}")

# Prop. E19(a): RM upper end = max_n (D + Mbar A(n)) <= max D + Mbar sup A, equality iff a common maximiser.
ok19a, n_eq, n_strict = True, 0, 0
for rep in range(40):
    a, b, n = rand_margins(3, 2, 2, 3, 10, 30)
    Mbar = rng.uniform(0.3, 1.5)
    rm = rm_range(a, b, Mbar, 1, 2, [(0, 1)], logical=False)
    L = CountLP(a, b)
    Dpost, Dpre = did_vec(L, 0, 1, 0, 1, 2), did_vec(L, 0, 1, 0, 0, 1)
    Dr, Pr = L.range(Dpost), L.range(Dpre)
    Asup = max(abs(Pr[0]), abs(Pr[1]))
    td = (Dr[0] - Mbar * Asup, Dr[1] + Mbar * Asup)
    ok19a &= rm[0] >= td[0] - 1e-8 and rm[1] <= td[1] + 1e-8
    # common maximiser test (upper end): exists n with D = max D and |Dpre| = Asup
    common = False
    for sg in (1, -1):
        Lc = CountLP(a, b); Lc.add_ge(Dpost, Dr[1] - 1e-9)
        val = Lc.opt(sg * Dpre, -1)[0]
        common |= val is not None and val >= Asup - 1e-7
    eq = abs(rm[1] - td[1]) < 1e-7
    ok19a &= (eq == common)
    n_eq += eq; n_strict += (not eq)
check("E19 [V] Prop. E19(a): RM set in trend-dial set with M = Mbar sup|D_pre|; upper ends equal iff a common maximiser exists",
      ok19a and n_eq > 0 and n_strict > 0, f"equal {n_eq}, strictly inside {n_strict}")

# Prop. E19(b): RM uninformative (= no-assumption set) iff the vacuity condition holds at maximisers of B_r1.
ok19b, n_vac, n_inf = True, 0, 0
for rep in range(40):
    a, b, n = rand_margins(3, 2, 2, 3, 10, 30)
    Mbar = rng.choice([0.2, 0.5, 1.0, 3.0])
    rm = rm_range(a, b, Mbar, 1, 2, [(0, 1)], logical=True)
    L0, att0, *_ = causal_lp(a, b, 1, 2)
    na = L0.range(att0)
    # criterion for the upper end: exists n with B_r1(n) = max and |thetaPT(n)| <= Mbar |Dpre(n)| (theta = 0 allowed),
    # thetaPT = B_r1 - D = B_r,base + B_s,post - B_s,base
    L = CountLP(a, b)
    Br1 = rate_vec(L, 0, 0, 2); Bmax = L.range(Br1)[1]
    thPT = rate_vec(L, 0, 0, 1) + rate_vec(L, 1, 0, 2) - rate_vec(L, 1, 0, 1)
    Dpre = did_vec(L, 0, 1, 0, 0, 1)
    crit = False
    for sg in (1, -1):
        Lc = CountLP(a, b); Lc.add_ge(Br1, Bmax - 1e-9)
        Lc.add_ge(Mbar * sg * Dpre - thPT, 0); Lc.add_ge(Mbar * sg * Dpre + thPT, 0)   # |thetaPT| <= Mbar A
        crit |= Lc.feasible()
    if rm is None:                 # RM (with the rate limits) refuted: theta_PT too far outside [0, 1]
        n_ref19 = globals().get("n_ref19", 0) + 1; globals()["n_ref19"] = n_ref19
        ok19b &= not crit
        continue
    eq = abs(rm[1] - na[1]) < 1e-7
    ok19b &= (eq == crit)
    n_vac += eq; n_inf += (not eq)
check("E19 [V] Prop. E19(b): RM upper end equals the no-assumption upper end iff the vacuity LP is feasible",
      ok19b and n_vac > 0 and n_inf > 0, f"vacuous {n_vac}, informative {n_inf}, refuted {globals().get('n_ref19', 0)}")

# Prop. E19(c): margins only, given the base gap g0 the allowance can reach Mbar max(|g0 - G-_pre|, |g0 - G+_pre|) >= Mbar w_pre/2;
# point-identified pre-trends reduce RM to the standard form D-set (+) [-Mbar|Dpre|, Mbar|Dpre|].
ok19c = True
for rep in range(30):
    a, b, n = rand_margins(3, 2, 2, 3, 10, 30)
    L = CountLP(a, b)
    gpre = rate_vec(L, 0, 0, 0) - rate_vec(L, 1, 0, 0); gb = rate_vec(L, 0, 0, 1) - rate_vec(L, 1, 0, 1)
    Gp, Gb = L.range(gpre), L.range(gb)
    for g0 in np.linspace(Gb[0], Gb[1], 5):
        Lg = CountLP(a, b); Lg.add_eq(gb, g0)
        r = Lg.range(gb - gpre)                  # = D_{pre->base}
        Amax = max(abs(r[0]), abs(r[1]))
        ok19c &= abs(Amax - max(abs(g0 - Gp[0]), abs(g0 - Gp[1]))) < 1e-8 and Amax >= (Gp[1] - Gp[0]) / 2 - 1e-9
    # pinned pre-trend: fix the pre-period table to the truth, base and post free
    pinsP = []
    for i in range(3):
        v = np.zeros(L.nn); v[L.ix(0, i, 0, 0)] = 1; pinsP += [(v, n[0, i, 0, 0]), (-v, -n[0, i, 0, 0])]
        v = np.zeros(L.nn); v[L.ix(1, i, 0, 0)] = 1; pinsP += [(v, n[1, i, 0, 0]), (-v, -n[1, i, 0, 0])]
    rmP = rm_range(a, b, 0.7, 1, 2, [(0, 1)], logical=False, extra=pinsP)
    LP_ = CountLP(a, b)
    for rr in pinsP:
        LP_.add_ge(*rr)
    Dpost = LP_.range(did_vec(LP_, 0, 1, 0, 1, 2)); Dp = true_val(did_vec(LP_, 0, 1, 0, 0, 1), n)
    ok19c &= np.allclose(rmP, (Dpost[0] - 0.7 * abs(Dp), Dpost[1] + 0.7 * abs(Dp)), atol=1e-8)
check("E19 [V] Prop. E19(c): margins only, allowance reaches Mbar max|g0 - G_pre| >= Mbar w_pre/2; pinned pre-trends give standard RR", ok19c)

# RR consecutive form vs NS levels form with two pre-periods: they differ ([I]).
a, b, n = rand_margins(3, 2, 2, 4, 10, 30)
L = CountLP(a, b)
rr_cons = rm_range(a, b, 1.0, 2, 3, [(0, 1), (1, 2)], logical=False)
# NS form: |b| <= Mbar max_{t <= base} |D_{t -> base}|, i.e. pairs (0,2), (1,2)
rr_lev = rm_range(a, b, 1.0, 2, 3, [(0, 2), (1, 2)], logical=False)
check("E19 [I] Def. RM: Rambachan-Roth consecutive differences and NS's 'levels relative to base' differ with two pre-periods",
      not np.allclose(rr_cons, rr_lev), f"consecutive [{rr_cons[0]:.3f},{rr_cons[1]:.3f}] vs levels [{rr_lev[0]:.3f},{rr_lev[1]:.3f}]")


# =====================================================================================
# E20. Smoothness (Prop. E20): one pre pair, one post: ATT in [G-_1 - 2G+_0 + G-_{-1} - M, G+_1 - 2G-_0 + G+_{-1} + M]
#      (margins only, no rate limits); width exceeds the PT width by w_0 + w_{-1} + 2M.
# =====================================================================================
ok20 = True
for rep in range(30):
    a, b, n = rand_margins(3, 2, 2, 3, 10, 30)
    M = rng.uniform(0, 0.1)
    L, att, bvec, e, th, S = causal_lp(a, b, 1, 2)
    L.extra_bounds[0] = (None, None); L.Aub = []; L.bub = []
    Dm1 = L.pad(did_vec(L, 0, 1, 0, 1, 0))          # delta_{-1} = gap_{-1} - gap_0
    # |(delta_1 - delta_0) - (delta_0 - delta_{-1})| = |b + delta_{-1}| <= M
    L.add_le(bvec + Dm1, M); L.add_ge(bvec + Dm1, -M)
    got = L.range(att)
    Lg = CountLP(a, b)
    G = [Lg.range(rate_vec(Lg, 0, 0, t) - rate_vec(Lg, 1, 0, t)) for t in range(3)]
    cf = (G[2][0] - 2 * G[1][1] + G[0][0] - M, G[2][1] - 2 * G[1][0] + G[0][1] + M)
    ok20 &= np.allclose(got, cf, atol=1e-8)
    w = [g[1] - g[0] for g in G]
    ok20 &= abs((cf[1] - cf[0]) - ((w[2] + w[1]) + (w[1] + w[0] + 2 * M))) < 1e-9
check("E20 [V] Prop. E20: smoothness with set-valued pre-trends: closed form = LP; extra width w_0 + w_-1 + 2M over PT", ok20)


# =====================================================================================
# E21. MTR / ND / PS (Prop. E21(c)): R x C, comparison = rest: U = B^U_r1 - max(B^L_r0, B_rest1(at B^U_r1) - Lam);
#      comparison = a named row (R = 3): U = B^U_r1 - max(B^L_r0, B^L_s1 - Lam), the ends being jointly attainable (6a L2.1(b)).
# =====================================================================================
ok21a, ok21b, n21_strict, n21 = True, True, 0, 0
for rep in range(40):
    a, b, n = rand_margins(4, 3, 3, 2)
    Lam = rng.choice([0.0, 0.03])
    for comp in ("rest", 1):
        L, att, bvec, e, th, S = causal_lp(a, b, 0, 1, comp=comp, Lam=Lam)
        Br1 = L.pad(rate_vec(L, 0, 0, 1)); Br0 = L.pad(rate_vec(L, 0, 0, 0))
        Bc1 = L.pad(rate_vec(L, comp, 0, 1) if comp != "rest" else rest_vec(L, 0, 0, 1))
        L.add_le(e(th) - Br1, 0)                 # MTR
        L.add_ge(e(th) - Br0, 0)                 # ND
        L.add_ge(e(th) - Bc1 + e(S), 0)          # PS on the untreated comparison rate
        got = L.range(att)
        Lc = CountLP(a, b)
        BU = Lc.range(rate_vec(Lc, 0, 0, 1))[1]; BL0 = Lc.range(rate_vec(Lc, 0, 0, 0))[0]
        if comp == "rest":
            Lx = CountLP(a, b); Lx.add_ge(rate_vec(Lx, 0, 0, 1), BU - 1e-10)
            brest = Lx.range(rest_vec(Lx, 0, 0, 1))
            U = BU - max(BL0, brest[0] - Lam)
            if U < -1e-9:
                ok21a &= got is None
            else:
                ok21a &= got is not None and abs(got[0]) < 1e-7 and abs(got[1] - U) < 1e-6
        else:
            BLs = Lc.range(rate_vec(Lc, 1, 0, 1))[0]
            U = BU - max(BL0, BLs - Lam)
            if U < -1e-9:
                ok21b &= got is None
            else:
                ok21b &= got is not None and abs(got[0]) < 1e-7 and abs(got[1] - U) < 1e-6
                n21 += 1
check("E21 [V] Prop. E21(c): MTR+ND+PS with spillover dial, R x C, comparison = rest: ATT set = [0, U] closed form = LP", ok21a)
check("E21 [V] Prop. E21(c): comparison = a named row (R = 3): U = B^U_r1 - max(B^L_r0, B^L_s1 - Lam) = LP (margins only: the two extremes are jointly attained)",
      ok21b and n21 > 0, f"{n21} non-empty cases")


# =====================================================================================
# E22. Spillovers (Prop. E22): pre-periods carry no information on S. Two worlds with identical observed
#      tables in every period and the same parallel trend (b = 0), different S: ATTs differ by the S gap.
#      Naive "RM calibrated on pre-periods covers ATT" fails when S != 0; adding [-Lam, Lam] restores it.
# =====================================================================================
a, b, n = rand_polar(4, 3, share=(0.0, 1.0))           # pure units: every rate identified
L = CountLP(a, b)
Dt = true_val(did_vec(L, 0, 1, 0, 1, 2), n)
Dpre_t = true_val(did_vec(L, 0, 1, 0, 0, 1), n)
att_A, att_B = Dt + 0.0, Dt + 0.08            # b = 0, S = 0 vs S = 0.08 (same observables)
ok22 = abs((att_B - att_A) - 0.08) < 1e-14
Mbar = 0.5
rm0 = rm_range(a, b, Mbar, 1, 2, [(0, 1)], logical=True, Lam=0.0)
rmL = rm_range(a, b, Mbar, 1, 2, [(0, 1)], logical=True, Lam=0.08)
# a truth at the RM limit with S = .08: b = -Mbar|Dpre|, ATT = D + S - b
att_true = Dt + 0.08 + Mbar * abs(Dpre_t)
Bobs = true_val(rate_vec(L, 0, 0, 2), n)
feasible_truth = 0 <= Bobs - att_true <= 1
ok22b = feasible_truth and att_true > rm0[1] + 1e-9 and rmL[0] - 1e-9 <= att_true <= rmL[1] + 1e-9
check("E22 [V] Prop. E22(a): identical observables in all periods, same b, different S -> ATTs differ by the S gap", ok22)
check("E22 [F] Prop. E22(b): RM from pre-periods alone can exclude the true ATT when S != 0; RM (+) [-Lam, Lam] covers it",
      ok22b, f"truth {att_true:.3f}; RM(Lam=0) [{rm0[0]:.3f},{rm0[1]:.3f}]; RM(Lam=.08) [{rmL[0]:.3f},{rmL[1]:.3f}]")


# =====================================================================================
# E24. Unions across branches and periods (Prop. E24): the union of branch confidence sets covers at
#      >= 1 - alpha; reporting only the plug-in branch's interval under-covers.
# =====================================================================================
reps, cov_u, cov_p = 20000, 0, 0
z = norm.ppf(0.975)
theta = 0.5                                # true branch: theta >= 0 (sigma = 1)
th = theta + rng.normal(size=reps)
lo_p, hi_p = np.maximum(0, th - z), th + z                 # '+' branch CI (empty if th + z < 0)
ne_p = th + z >= 0
lo_m, hi_m = th - z, np.minimum(0, th + z)                 # '-' branch CI (empty if th - z > 0)
ne_m = th - z <= 0
in_p = ne_p & (lo_p <= theta) & (theta <= hi_p)
in_m = ne_m & (lo_m <= theta) & (theta <= hi_m)
cov_u = np.mean(in_p | in_m)
cov_p = np.mean(np.where(th >= 0, in_p, in_m))
check("E24 [F] Prop. E24: union of branch intervals covers (>= .95); the plug-in branch's interval alone under-covers",
      cov_u >= 0.95 - 0.005 and cov_p < 0.9, f"union {cov_u:.3f}, plug-in branch {cov_p:.3f}")


# =====================================================================================
# E25. Joint confidence region for model-based pre-trends (Prop. E25 = E19(d)): one joint region covers
#      the RM lower end; the plug-in pre-trend with an interval for the post DiD alone under-covers.
#      g_t ~ N(g_t, s^2) independent estimates of three periods' gaps; D_pre = g0 - g_-1, D_post = g1 - g0.
# =====================================================================================
reps, s_ = 20000, 1.0
gtrue = np.array([0.0, 3.0, 4.0]); Mbar = 1.0
Dp, Dq = gtrue[1] - gtrue[0], gtrue[2] - gtrue[1]
att_true = Dq - Mbar * abs(Dp)               # ATT at the lower end of the RM set (b at its limit)
gh = gtrue + s_ * rng.normal(size=(reps, 3))
Dph, Dqh = gh[:, 1] - gh[:, 0], gh[:, 2] - gh[:, 1]
naive_lo = Dqh - Mbar * np.abs(Dph) - z * s_ * np.sqrt(2)
naive_hi = Dqh + Mbar * np.abs(Dph) + z * s_ * np.sqrt(2)
kq = np.sqrt(chi2.ppf(0.95, 3))
cs = [np.array([Mbar, -1 - Mbar, 1.0]), np.array([-Mbar, -1 + Mbar, 1.0])]   # D_post -/+ Mbar D_pre
j_lo = np.min([gh @ c - kq * s_ * np.linalg.norm(c) for c in cs], axis=0)
cu = [np.array([-Mbar, -1 + Mbar, 1.0]), np.array([Mbar, -1 - Mbar, 1.0])]   # D_post +/- Mbar D_pre
j_hi = np.max([gh @ c + kq * s_ * np.linalg.norm(c) for c in cu], axis=0)
cov_n = np.mean((naive_lo <= att_true) & (att_true <= naive_hi))
cov_j = np.mean((j_lo <= att_true) & (att_true <= j_hi))
check("E25 [F] Prop. E25 (= E19(d)): joint region over all periods covers the RM ATT set; plug-in pre-trend + post interval under-covers",
      cov_j >= 0.95 - 0.005 and cov_n < 0.93, f"joint {cov_j:.3f}, plug-in {cov_n:.3f}")


# =====================================================================================
# E27. Integer versions (Theorem E27): R = 2, change rows on one column and one period pair over a
#      laminar family of unit sets (rows rounded): the LP is integral (LP = MILP for random objectives).
#      Two overlapping period pairs are outside the theorem (searched; reported).
# =====================================================================================
def laminar_family(p):
    fam = [list(range(p))]
    cut = rng.integers(1, p)
    fam += [list(range(cut)), list(range(cut, p))]
    if cut > 1:
        fam.append([0])
    return fam


ok25, n25, gap_star, n_star = True, 0, 0, 0
for rep in range(60):
    p, R, C, T = 4, 2, 3, 4
    a, b, n = rand_margins(p, R, C, T, 4, 9)
    for mode in ("one", "path", "base2", "star3"):
        L = CountLP(a, b)
        pairs = {"one": [(0, 1)], "path": [(0, 1), (1, 2), (2, 3)], "base2": [(1, 0), (1, 2)],
                 "star3": [(0, 1), (0, 2), (0, 3)]}[mode]
        for (t0, t1) in pairs:
            fam = laminar_family(p)                          # each pair its own laminar family
            for U in fam:
                v = np.zeros(L.nn)
                for i in U:
                    v[L.ix(t1, i, 0, 0)] += 1; v[L.ix(t0, i, 0, 0)] -= 1
                tv = true_val(v, n)
                L.add_ge(v, tv - rng.integers(0, 2)); L.add_le(v, tv + rng.integers(0, 2))   # integer rhs
        for k in range(3):
            obj = rng.normal(size=L.nn)
            cr, ir = L.range(obj), L.irange(obj)
            same = abs(cr[0] - ir[0]) < 1e-6 and abs(cr[1] - ir[1]) < 1e-6
            if mode != "star3":
                ok25 &= same; n25 += 1
            else:
                gap_star += (not same); n_star += 1
check("E27 [V] Theorem E27: R=2, change rows on one column, period pairs forming paths, laminar unit sets: LP = MILP (integral)", ok25,
      f"{n25} objectives (one pair; consecutive path; two pairs sharing a base); three pairs sharing a base (outside the theorem): LP != MILP in {gap_star}/{n_star}")


# =====================================================================================
# G. Gap register: for each item, the naive statement fails when its assumption is dropped [F], and the
#    stated closing device holds [V] where one is given.
# =====================================================================================
def dd_rate(V, Y, N):
    """Duncan-Davis range of the group-0 rate over units (V group-0 counts, Y outcome counts, N sizes)."""
    lo = np.sum(np.maximum(0, V + Y - N)) / V.sum(); hi = np.sum(np.minimum(V, Y)) / V.sum()
    return lo, hi


# G1. Entry. Naive: the DiD of period-own group rates (and its Delta-0 range from period-own margins)
#     describes the stayers. Closing: entrant rows with known (or interval) entrant counts.
fail1, ok1, n1 = 0, True, 0
for rep in range(200):
    p = 5
    big = rng.integers(150, 300, p); small = rng.integers(1, 6, p)    # near-pure units: narrow ranges
    Ov = np.where(np.arange(p) % 2 == 0, big, small); On = np.where(np.arange(p) % 2 == 0, small, big)
    Ev = np.where(np.arange(p) % 2 == 0, rng.integers(40, 120, p), 0); En = np.where(np.arange(p) % 2 == 0, 0, rng.integers(40, 120, p))
    y0v = rng.binomial(Ov, 0.2); y0n = rng.binomial(On, 0.3)          # stayers' yes, t = 0
    y1v = rng.binomial(Ov, 0.5); y1n = rng.binomial(On, 0.35)         # stayers' yes, t = 1
    e1v = rng.binomial(Ev, 0.95); e1n = rng.binomial(En, 0.05)        # entrants vote differently
    did_stay = (y1v.sum() - y0v.sum()) / Ov.sum() - (y1n.sum() - y0n.sum()) / On.sum()
    # naive Delta-0 range on period-own populations
    g0 = dd_rate(Ov, y0v + y0n, Ov + On); g0n = dd_rate(On, y0v + y0n, Ov + On)
    g1 = dd_rate(Ov + Ev, y1v + y1n + e1v + e1n, Ov + On + Ev + En); g1n = dd_rate(On + En, y1v + y1n + e1v + e1n, Ov + On + Ev + En)
    naive = ((g1[0] - g1n[1]) - (g0[1] - g0n[0]), (g1[1] - g1n[0]) - (g0[0] - g0n[1]))
    fail1 += not (naive[0] - 1e-12 <= did_stay <= naive[1] + 1e-12)
    # closing device: stayer cells per period + entrant cells at t = 1, entrant counts known
    a_ = np.column_stack([Ov, On]); b0 = np.column_stack([y0v + y0n, Ov + On - y0v - y0n])
    # period-1 stayer table margins are unknown: the period-1 outcome total mixes entrants. LP variables:
    # stayer cells n[t,i,r,c] (fixed population, rows Ov/On) + entrant yes cells e_v, e_n per unit.
    nv = 2 * p * 4 + 2 * p
    ixs = lambda t, i, r, c: ((t * p + i) * 2 + r) * 2 + c
    ixe = lambda i, r: 2 * p * 4 + i * 2 + r
    Aeq, beq = [], []
    for i in range(p):
        for t in (0, 1):
            for r, tot in ((0, Ov[i]), (1, On[i])):
                v = np.zeros(nv); v[[ixs(t, i, r, 0), ixs(t, i, r, 1)]] = 1; Aeq.append(v); beq.append(tot)
        v = np.zeros(nv); v[[ixs(0, i, 0, 0), ixs(0, i, 1, 0)]] = 1; Aeq.append(v); beq.append(y0v[i] + y0n[i])
        v = np.zeros(nv); v[[ixs(1, i, 0, 0), ixs(1, i, 1, 0), ixe(i, 0), ixe(i, 1)]] = 1
        Aeq.append(v); beq.append(y1v[i] + y1n[i] + e1v[i] + e1n[i])
    bnds = [(0, None)] * (2 * p * 4) + sum([[(0, Ev[i]), (0, En[i])] for i in range(p)], [])
    obj = np.zeros(nv)
    for i in range(p):
        obj[ixs(1, i, 0, 0)] += 1 / Ov.sum(); obj[ixs(0, i, 0, 0)] -= 1 / Ov.sum()
        obj[ixs(1, i, 1, 0)] -= 1 / On.sum(); obj[ixs(0, i, 1, 0)] += 1 / On.sum()
    lo = linprog(obj, A_eq=np.array(Aeq), b_eq=beq, bounds=bnds, method="highs").fun
    hi = -linprog(-obj, A_eq=np.array(Aeq), b_eq=beq, bounds=bnds, method="highs").fun
    ok1 &= lo - 1e-9 <= did_stay <= hi + 1e-9
    n1 += 1
check("G1 [F] entry: period-own Delta-0 range can miss the stayers' DiD; the entrant-row polytope contains it",
      fail1 > 0 and ok1, f"naive misses in {fail1}/{n1}; polytope covers in all")

# G2. Exit. Naive: the base-period pool's rate equals the stayers' rate. Dial: B_pool - B_O = (1 - pi)(B_E - B_O).
ok2, fail2 = True, 0
for rep in range(300):
    O, E = rng.integers(10, 200, 2)
    yO, yE = rng.integers(0, O + 1), rng.integers(0, E + 1)
    BO, BE, Bpool, pi = yO / O, yE / E, (yO + yE) / (O + E), O / (O + E)
    ok2 &= abs((Bpool - BO) - (1 - pi) * (BE - BO)) < 1e-14
    fail2 += abs(Bpool - BO) > 1e-9
check("G2 [F] exit: pool rate != stayers' rate unless leavers vote like stayers; the dial identity holds exactly",
      ok2 and fail2 > 0, f"naive equality fails in {fail2}/300")

# G3. Aging-in (2d C2). Naive: the DiD on a common post-period base equals the stayers' attitude DiD.
fail3 = 0
for rep in range(300):
    O = rng.integers(20, 500, 2).astype(float); A = rng.integers(0, 200, 2).astype(float) * np.array([1.5, 1.0])
    BO0, BO1 = rng.uniform(0.01, 0.2, 2), rng.uniform(0, 0.8, 2)
    om = O / (O + A)
    did_orig = (BO1[0] - om[0] * BO0[0]) - (BO1[1] - om[1] * BO0[1])     # aged-in vote like stayers at t = 1
    did_O = (BO1[0] - BO0[0]) - (BO1[1] - BO0[1])
    fail3 += abs(did_orig - did_O - ((1 - om[0]) * BO0[0] - (1 - om[1]) * BO0[1])) < 1e-12 and abs(did_orig - did_O) > 1e-6
check("G3 [F] aging-in: common-base DiD != stayers' DiD; gap = (1-w_v)B^O_v0 - (1-w_n)B^O_n0 (2d C2)", fail3 > 250, f"{fail3}/300")

# G4. Stayers vs leavers (selection on staying). Naive: parallel trends in the base population implies
#     parallel trends among stayers. Counterexample by arithmetic: untreated rates for stayers (S) and leavers (L).
pi_v, pi_n = 0.6, 0.9                                 # staying shares differ by group
# untreated changes: stayers v +.10, leavers v -.05 ; stayers n +.04, leavers n +.04 (pool trends equal?)
chg = {"vS": 0.10, "vL": -0.05, "nS": 0.04, "nL": 0.04}
pool_v = pi_v * chg["vS"] + (1 - pi_v) * chg["vL"]; pool_n = pi_n * chg["nS"] + (1 - pi_n) * chg["nL"]
b_pool, b_stay = pool_v - pool_n, chg["vS"] - chg["nS"]
check("G4 [F] stayers vs leavers: PT in the pool (b = 0) with b != 0 among stayers when staying differs by group",
      abs(b_pool) < 1e-12 and abs(b_stay) > 0.05, f"b_pool {b_pool:.3f}, b_stayers {b_stay:.3f}")

# G5. Interval bases. Naive: a census-projected base N_hat treated as exact. Closing: N in [N_L, N_U], N >= Y, N >= V.
fail5, ok5, n5 = 0, True, 0
for rep in range(300):
    p = 6
    V = rng.integers(20, 80, p); M = rng.integers(20, 80, p); N = V + M
    yv, ym = rng.binomial(V, 0.2), rng.binomial(M, 0.8); Y = yv + ym
    Nhat = np.maximum(np.round(N * rng.uniform(0.6, 0.9, p)).astype(int), np.maximum(V, Y))     # census undercount
    truth = yv.sum() / V.sum()
    lo, hi = dd_rate(V, Y, Nhat)
    fail5 += not (lo - 1e-12 <= truth <= hi + 1e-12)
    NL, NU = np.maximum(np.minimum(Nhat, N), np.maximum(V, Y)), np.maximum(Nhat, N) + 5   # an interval containing N
    # LP over N_i in [NL, NU]: the veteran-yes count y_i in [max(0, V+Y-N), min(V, Y)]; the lower end decreases in N
    lo_i, hi_i = np.sum(np.maximum(0, V + Y - NU)) / V.sum(), np.sum(np.minimum(V, Y)) / V.sum()
    ok5 &= lo_i - 1e-12 <= truth <= hi_i + 1e-12
    n5 += 1
check("G5 [F] interval bases: DD with a census-projected base taken as exact can exclude the truth; the interval base covers it",
      fail5 > 0 and ok5, f"naive misses in {fail5}/{n5}")

# G6. Changing unit boundaries. Blocks k = 1, 2, 3 (fixed people); period-0 units {1,2},{3}; period-1 units {1},{2,3}.
#     Naive: allocate period-0 counts to period-1 units by fixed shares (here 1/2 of unit {1,2} to each side),
#     then use the fixed-population Delta-0 bound. Closing: coarsen to the common unit {1,2,3} (valid, wider).
fail6, ok6, n6 = 0, True, 0
for rep in range(300):
    V = np.array([95, 5, 50]) + rng.integers(0, 5, 3); M = np.array([5, 95, 50]) + rng.integers(0, 5, 3)   # near-pure blocks
    pol = lambda: rng.choice([0.05, 0.95], 3)                                                         # polarised rates
    y0v, y0m = rng.binomial(V, pol()), rng.binomial(M, pol())
    y1v, y1m = rng.binomial(V, pol()), rng.binomial(M, pol())
    truth = (y1v.sum() - y0v.sum()) / V.sum() - (y1m.sum() - y0m.sum()) / M.sum()
    # coarse: one unit
    a_c = np.array([[V.sum(), M.sum()]]); b_c = np.array([[[ (y0v + y0m).sum(), (V + M).sum() - (y0v + y0m).sum()]],
                                                         [[ (y1v + y1m).sum(), (V + M).sum() - (y1v + y1m).sum()]]])
    Lc = CountLP(a_c, b_c); rc = Lc.range(did_vec(Lc, 0, 1, 0, 0, 1))
    ok6 &= rc[0] - 1e-9 <= truth <= rc[1] + 1e-9
    # naive: period-1 units U1 = {1}, U2 = {2,3}; period-0 counts of {1,2} split 50/50 between U1 and U2
    Y12 = y0v[0] + y0v[1] + y0m[0] + y0m[1]; Y3 = y0v[2] + y0m[2]
    aU = np.array([[V[0], M[0]], [V[1] + V[2], M[1] + M[2]]])
    Y0U = np.array([0.5 * Y12, 0.5 * Y12 + Y3])
    Y0U = np.minimum(Y0U, aU.sum(1))
    Y1U = np.array([y1v[0] + y1m[0], y1v[1] + y1v[2] + y1m[1] + y1m[2]])
    bU = np.array([np.column_stack([Y0U, aU.sum(1) - Y0U]), np.column_stack([Y1U, aU.sum(1) - Y1U])])
    Ln = CountLP(aU, bU); rn = Ln.range(did_vec(Ln, 0, 1, 0, 0, 1))
    n6 += 1
    if rn is None or not (rn[0] - 1e-9 <= truth <= rn[1] + 1e-9):
        fail6 += 1
check("G6 [F] changing boundaries: share-allocated margins can exclude the truth; coarsening to common units covers it",
      fail6 > 0 and ok6, f"naive misses or infeasible in {fail6}/{n6}")

# G7. Uncertain margins. (a) Naive: expected (matched) group counts taken as exact can exclude the truth; the union
#     over a margin set containing the true counts covers (R25(a)). (b) A composition common to both periods
#     enters once: the union over V (same in both periods) is inside the union over (V_0, V_1) separately.
#     (c) R26 on changes: E[Delta T | q] = qbar' dQ qbar + Var(X|q) (dQvv - 2 dQvn + dQnn) (Berkson), exact.
fail7, ok7a, n7, n7_strict = 0, True, 0, 0
for rep in range(200):
    p = 3
    Ncert = rng.integers(10, 30, p); Nunc = rng.integers(5, 20, p); M = rng.integers(20, 60, p)
    q = rng.uniform(0.2, 0.8, p)
    Vtrue = Ncert + rng.binomial(Nunc, q); Vexp = Ncert + Nunc * q
    N = Ncert + Nunc + M
    yv0, yv1 = rng.binomial(Vtrue, 0.02), rng.binomial(Vtrue, 0.98)       # polarised rates: narrow ranges
    ym0, ym1 = rng.binomial(N - Vtrue, 0.02), rng.binomial(N - Vtrue, 0.03)
    Y0, Y1 = yv0 + ym0, yv1 + ym1
    truth = (yv1.sum() - yv0.sum()) / Vtrue.sum() - ((Y1 - yv1).sum() - (Y0 - yv0).sum()) / (N - Vtrue).sum()
    def did_range(V0, V1):
        # fixed population needs V0 == V1; for the 'separate' union allow them to differ (two populations)
        a0 = np.column_stack([V0, N - V0]); a1 = np.column_stack([V1, N - V1])
        g = []
        for (aa, Y) in ((a0, Y0), (a1, Y1)):
            lo = np.sum(np.maximum(0, aa[:, 0] + Y - N)) / aa[:, 0].sum() - np.sum(np.minimum(aa[:, 1], Y)) / aa[:, 1].sum()
            hi = np.sum(np.minimum(aa[:, 0], Y)) / aa[:, 0].sum() - np.sum(np.maximum(0, aa[:, 1] + Y - N)) / aa[:, 1].sum()
            g.append((lo, hi))
        return g[1][0] - g[0][1], g[1][1] - g[0][0]
    naive = did_range(Vexp, Vexp)
    fail7 += not (naive[0] - 1e-12 <= truth <= naive[1] + 1e-12)
    # margin set: every V with Ncert <= V <= Ncert + Nunc (probability 1 here) -> union over integer V vectors
    grid = [range(Ncert[i], Ncert[i] + Nunc[i] + 1, max(1, Nunc[i] // 4)) for i in range(p)]
    grid = [sorted(set(list(gr) + [Vtrue[i]])) for i, gr in enumerate(grid)]
    common = [did_range(np.array(v), np.array(v)) for v in itertools.product(*grid)]
    uc = (min(c[0] for c in common), max(c[1] for c in common))
    ok7a &= uc[0] - 1e-12 <= truth <= uc[1] + 1e-12
    if rep < 20:
        sep = [did_range(np.array(v0), np.array(v1)) for v0 in itertools.product(*grid) for v1 in itertools.product(*grid)]
        us = (min(c[0] for c in sep), max(c[1] for c in sep))
        ok7a &= us[0] <= uc[0] + 1e-12 and uc[1] <= us[1] + 1e-12
        n7_strict += (us[1] - us[0]) > (uc[1] - uc[0]) + 1e-9
    n7 += 1
check("G7 [F] uncertain margins: expected counts taken as exact can miss the truth; the union over the margin set covers it",
      fail7 > 0 and ok7a, f"naive misses in {fail7}/{n7}")
check("G7 [V] uncertain margins: a composition common to both periods enters once (union over V inside the union over (V0, V1))",
      ok7a and n7_strict > 0, f"strictly narrower in {n7_strict}/20")
ok7c, bias7 = True, 0.0
for rep in range(50):
    Nm = 8
    qj = rng.uniform(0, 1, Nm); dQ_ = rng.uniform(-0.1, 0.1, 3)
    # exact law of X = (sum of independent Bernoulli(qj)) / Nm by enumeration of the Poisson-binomial
    pmf = np.array([1.0])
    for qq in qj:
        pmf = np.convolve(pmf, [1 - qq, qq])
    xs_ = np.arange(Nm + 1) / Nm
    f = lambda X: dQ_[0] * X ** 2 + 2 * dQ_[1] * X * (1 - X) + dQ_[2] * (1 - X) ** 2
    exact = np.sum(pmf * f(xs_))
    qbar = qj.mean(); varX = np.sum(qj * (1 - qj)) / Nm ** 2
    pred = f(qbar) + varX * (dQ_[0] - 2 * dQ_[1] + dQ_[2])
    ok7c &= abs(exact - pred) < 1e-13
    bias7 = max(bias7, abs(exact - f(qbar)))
check("G7 [F] R26 on changes: E[Delta T|q] = qbar'dQ qbar + Var(X|q)(dQvv - 2dQvn + dQnn) exactly; ignoring it is biased", ok7c and bias7 > 1e-5,
      f"max bias of the uncorrected fit {bias7:.2e}")

# G8. Covariates. Naive: unconditional DiD = ATT when PT holds only within covariate strata. The stratum-weighted
#     conditional DiD (weights = treated-group shares by stratum) recovers the ATT.
wz_v, wz_n = np.array([0.8, 0.2]), np.array([0.3, 0.7])            # group composition across strata
trend0 = np.array([0.02, 0.10])                                      # untreated trends by stratum (same for both groups)
att_z = np.array([0.05, 0.05])
D_unc = np.sum(wz_v * (trend0 + att_z)) - np.sum(wz_n * trend0)
D_cond = np.sum(wz_v * ((trend0 + att_z) - trend0))
att_ = np.sum(wz_v * att_z)
check("G8 [F] covariates: conditional PT within strata, unconditional DiD != ATT; treated-weighted stratum DiD = ATT",
      abs(D_unc - att_) > 0.01 and abs(D_cond - att_) < 1e-14, f"unconditional {D_unc:.3f}, conditional {D_cond:.3f}, ATT {att_:.3f}")

# G9. Time-varying group membership. Naive: the DiD of period-specific group rates equals the within-person
#     DiD of those in the group in both periods. Fails with switchers.
# exact shares (no sampling): v-stayers .34, v->n switchers .06, n-stayers .51, n->v switchers .09;
# switchers out of n had a low t = 0 rate (.20); n's t = 0 mean is .50; at t = 1 everyone takes the new group's rate
sh = {"vv": .34, "vn": .06, "nn": .51, "nv": .09}
r0 = {"vv": .30, "vn": .30, "nv": .20}; r0["nn"] = (.50 * .60 - sh["nv"] * r0["nv"]) / sh["nn"]
r1 = {"vv": .70, "nv": .70, "nn": .45, "vn": .45}
B = lambda keys, r: sum(sh[k] * r[k] for k in keys) / sum(sh[k] for k in keys)
did_groups = (B(["vv", "nv"], r1) - B(["vv", "vn"], r0)) - (B(["nn", "vn"], r1) - B(["nn", "nv"], r0))
did_within = (r1["vv"] - r0["vv"]) - (r1["nn"] - r0["nn"])
check("G9 [F] time-varying groups: DiD of period-specific group rates != within-person DiD of group stayers",
      abs(did_groups - did_within) > 0.01, f"period groups {did_groups:.3f}, stayers {did_within:.3f}")


# =====================================================================================
# W. Worked example of the plain-language section (§9): four precincts of 200 men with 10, 50, 150, 190
#    veterans; three elections (code t = 0, 1, 2 = document t = -1, 0, 1). Yes counts below.
# =====================================================================================
aW = np.array([[10, 190], [50, 150], [150, 50], [190, 10]])
yesW = np.array([[60, 62, 61, 56], [62, 66, 63, 61], [74, 86, 113, 120]])
bW = np.stack([np.column_stack([yesW[t], aW.sum(1) - yesW[t]]) for t in range(3)])
LW = CountLP(aW, bW)
from decimal import Decimal, ROUND_HALF_UP
r4 = lambda z: float(Decimal(repr(float(z) + (1e-12 if z >= 0 else -1e-12))).quantize(Decimal('0.0001'), ROUND_HALF_UP))
rW = lambda v, L=LW: tuple(r4(z) for z in L.range(v))
NW = aW.sum(1); pwW = NW / NW.sum(); xW = aW[:, 0] / NW; muW = np.sum(pwW * xW); vxW = np.sum(pwW * xW * xW) - muW ** 2
def er_nm(t0, t1):
    dT = (bW[t1, :, 0] - bW[t0, :, 0]) / NW
    cv = np.sum(pwW * xW * dT) - muW * np.sum(pwW * dT)
    return cv / vxW, cv / (muW * (1 - muW))
def cr_rowsW(t0, t1):
    ER, NM = er_nm(t0, t1); Dv = did_vec(LW, 0, 1, 0, t0, t1); sg = 1 if ER >= 0 else -1
    return [(sg * (-Dv), sg * (-ER)), (sg * Dv, sg * NM)]
def attW(M=None, Lam=0.0, rows=None, cl=False):
    Lc, attv, bvec, e, th, S = causal_lp(aW, bW, 1, 2, Lam=Lam, int_rows=rows)
    if M is not None:
        Lc.add_le(bvec, M); Lc.add_ge(bvec, -M)
    if cl:
        Br1 = Lc.pad(rate_vec(Lc, 0, 0, 2)); Br0 = Lc.pad(rate_vec(Lc, 0, 0, 1)); Bc1 = Lc.pad(rate_vec(Lc, 1, 0, 2))
        Lc.add_le(e(th) - Br1, 0); Lc.add_ge(e(th) - Br0, 0); Lc.add_ge(e(th) - Bc1 + e(S), 0)
    x = Lc.range(attv)
    return tuple(r4(z) for z in x)
crP, crQ = cr_rowsW(1, 2), cr_rowsW(0, 1)
LWc = CountLP(aW, bW)
for rr in crP:
    LWc.add_ge(*rr)
LWq = CountLP(aW, bW)
for rr in crQ:
    LWq.add_ge(*rr)
W = {
    "B_v,post": rW(rate_vec(LW, 0, 0, 2)),
    "D_post (margins)": rW(did_vec(LW, 0, 1, 0, 1, 2)),
    "D_pre (margins)": rW(did_vec(LW, 0, 1, 0, 0, 1)),
    "ER, NM on post changes": tuple(r4(z) for z in er_nm(1, 2)),
    "ER, NM on pre changes": tuple(r4(z) for z in er_nm(0, 1)),
    "D_post (CR on changes)": rW(did_vec(LWc, 0, 1, 0, 1, 2), LWc),
    "D_pre (CR on pre changes)": rW(did_vec(LWq, 0, 1, 0, 0, 1), LWq),
    "ATT, no assumption": attW(),
    "ATT, PT, margins only": attW(0.0),
    "ATT, PT + CR": attW(0.0, rows=crP),
    "ATT, dial M=.05 + CR": attW(0.05, rows=crP),
    "ATT, dial M=.05, spillover .02 + CR": attW(0.05, 0.02, rows=crP),
    "ATT, RM(1), margins only": tuple(r4(z) for z in rm_range(aW, bW, 1.0, 1, 2, [(0, 1)])),
    "ATT, RM(1), CR post only": tuple(r4(z) for z in rm_range(aW, bW, 1.0, 1, 2, [(0, 1)], extra=crP)),
    "ATT, RM(1), CR post and pre": tuple(r4(z) for z in rm_range(aW, bW, 1.0, 1, 2, [(0, 1)], extra=crP + crQ)),
    "ATT, MTR+ND+PS, margins only": attW(cl=True),
    "ATT, MTR+ND+PS + CR": attW(cl=True, rows=crP),
}
for k, v in W.items():
    print(f"      W  {k}: [{v[0]:.4f}, {v[1]:.4f}]")
check("W [I] worked example (§9): numbers as printed above, each rung narrows or reports no information as stated",
      W["ATT, RM(1), margins only"] == W["ATT, no assumption"] and W["ATT, PT + CR"] == W["D_post (CR on changes)"]
      and W["ATT, RM(1), CR post and pre"][0] > 0)


# =====================================================================================
# R. Checks added after the proof review (reviews/PHASE2e_proof_review.md, issues 1, 2, 12, 14).
# =====================================================================================
# Issue 1 (Cor. E17(b)): one 2x2 unit, a = (100, 100), 100 yes votes at both elections, one coupling row
# 3 y0 - 2 y1 = 140 (y_t = treated-yes count). The earlier sufficient condition holds; the sharp upper end under
# parallel trends (M = Lam = 0) is below max D; the exact criterion of E17(b) predicts this.
aR = np.array([[100, 100]]); bR = np.array([[[100, 100]], [[100, 100]]])
LR = CountLP(aR, bR)
row = np.zeros(LR.nn); row[LR.ix(0, 0, 0, 0)] = 3; row[LR.ix(1, 0, 0, 0)] = -2
Dv = did_vec(LR, 0, 1, 0, 0, 1)
LRr = CountLP(aR, bR); LRr.add_eq(row, 140)
Dmax = LRr.range(Dv)[1]
thPT = rate_vec(LRr, 0, 0, 0) + rate_vec(LRr, 1, 0, 1) - rate_vec(LRr, 1, 0, 0)
Bs1 = rate_vec(LRr, 1, 0, 1)
Lp_ = CountLP(aR, bR); Lp_.add_eq(row, 140); Lp_.add_ge(Dv, Dmax - 1e-9)
thp = Lp_.range(thPT); bsp = Lp_.range(Bs1)
Dmin = LRr.range(Dv)[0]
Lm_ = CountLP(aR, bR); Lm_.add_eq(row, 140); Lm_.add_le(Dv, Dmin + 1e-9)
thm = Lm_.range(thPT); bsm = Lm_.range(Bs1)
old_cond = thp[1] >= 0 and bsp[1] >= 0 and thm[0] <= 1 and bsm[0] <= 1
LcR, attR, bvecR, eR, thR, SR = causal_lp(aR, bR, 0, 1, int_rows=[(row, 140), (-row, -140)])
LcR.add_le(bvecR, 0.0); LcR.add_ge(bvecR, 0.0)
sharp = LcR.range(attR)
new_crit = thp[0] <= 1 + 1e-12 and thp[1] >= -1e-12          # some argmax n with thetaPT in [0, 1]
check("E17 [F] review issue 1: E17(b) as first stated (one-sided rate limits) fails with a coupling row; the exact criterion holds",
      old_cond and sharp[1] < Dmax - 1e-6 and not new_crit,
      f"max D {Dmax:.3f}, thetaPT at argmax {thp[0]:.3f}..{thp[1]:.3f}, sharp PT set [{sharp[0]:.3f},{sharp[1]:.3f}]")
# margins only: the inner conditions at n+ and n- imply the outer ones (n+ minimises, n- maximises thetaPT)
ok_m, n_m = True, 0
for rep in range(60):
    a, b, n = rand_margins(4, 2, 2, 2) if rep % 2 else rand_polar(4, 2)
    for comp in (1, "rest"):
        L = CountLP(a, b)
        Bc = (lambda t: rate_vec(L, comp, 0, t)) if comp != "rest" else (lambda t: rest_vec(L, 0, 0, t))
        D = rate_vec(L, 0, 0, 1) - rate_vec(L, 0, 0, 0) - Bc(1) + Bc(0)
        th = rate_vec(L, 0, 0, 0) + Bc(1) - Bc(0)
        Dr, thr = L.range(D), L.range(th)
        Lp = CountLP(a, b); Lp.add_ge(D, Dr[1] - 1e-9); Lm = CountLP(a, b); Lm.add_le(D, Dr[0] + 1e-9)
        ok_m &= abs(Lp.range(th)[0] - thr[0]) < 1e-7 and abs(Lm.range(th)[1] - thr[1]) < 1e-7
        n_m += 1
check("E17 [V] review issue 1: margins only, a maximiser of D minimises thetaPT and a minimiser maximises it (named row and rest)",
      ok_m, f"{n_m} cases")

# Issue 2 (Prop. E6(b)): four pure units; bins L = {group-1 unit of 20, group-2 unit of 80}, H = {group-1 unit of 80,
# group-2 unit of 20}; yes votes (10, 40, 40, 4) then (8, 40, 72, 16): Dbar+ < 0 <= cap, cap >= Dbar-, LP infeasible.
aE = np.array([[20, 0], [0, 80], [80, 0], [0, 20]])
yE = np.array([[10, 40, 40, 4], [8, 40, 72, 16]])
bE = np.stack([np.column_stack([yE[t], aE.sum(1) - yE[t]]) for t in range(2)])
LE = CountLP(aE, bE)
Lb_, Hb_ = [0, 1], [2, 3]
bg = lambda U: (rate_vec(LE, 0, 0, 1, U) - rate_vec(LE, 0, 0, 0, U)) - (rate_vec(LE, 1, 0, 1, U) - rate_vec(LE, 1, 0, 0, U))
Dbar = 0.5 * (bg(Hb_) + bg(Lb_))
lamE = aE[Hb_, 0].sum() / aE[Hb_].sum() - aE[Lb_, 0].sum() / aE[Lb_].sum()
dtauE = (yE[1, Hb_].sum() - yE[0, Hb_].sum()) / aE[Hb_].sum() - (yE[1, Lb_].sum() - yE[0, Lb_].sum()) / aE[Lb_].sum()
capE = dtauE / lamE
rgE = LE.range(Dbar)
LEb = CountLP(aE, bE); LEb.add_ge(lamE * Dbar, 0); LEb.add_le(lamE * Dbar, dtauE)
infeas = LEb.range(Dbar) is None
check("E6 [F] review issue 2: the emptiness criterion needs 'or Dbar+ < 0': here dtau_H > dtau_L and cap >= Dbar-, yet the set is empty",
      infeas and dtauE > 0 and capE >= rgE[0] and rgE[1] < 0 and max(0, rgE[0]) > min(rgE[1], capE),
      f"Dbar in [{rgE[0]:.3f},{rgE[1]:.3f}], lambda {lamE:.2f}, dtau_H - dtau_L {dtauE:.2f}, cap {capE:.3f}")

# Issue 12 (§9.2): RM(1) with CR on the pre-war change only vs on both changes.
pre_only = tuple(r4(z) for z in rm_range(aW, bW, 1.0, 1, 2, [(0, 1)], extra=crQ))
check("W [I] review issue 12: RM(1) with CR on the pre-war change alone gives [-0.4155, 0.7325]; on both changes [0.1352, 0.2783]",
      pre_only == (-0.4155, 0.7325) and W["ATT, RM(1), CR post and pre"] == (0.1352, 0.2783), f"pre only {pre_only}")

# Issue 14: the counts quoted in the note's gap-register table match this run.
doc = open(os.path.join(here, "PHASE2e_did_general.md")).read()
check("G [V] review issue 14: the gap-register counts quoted in PHASE2e (G1, G5, G6, G7) match this run",
      f"{fail1}/{n1}" in doc and f"{fail5}/{n5}" in doc and f"{fail6}/{n6}" in doc and f"{fail7}/{n7}" in doc,
      f"G1 {fail1}/{n1}, G5 {fail5}/{n5}, G6 {fail6}/{n6}, G7 {fail7}/{n7}")


# =====================================================================================
# IR. Checks added after the independent review (reviews/PHASE2e_independent_review.md, findings 1-12).
#     Cases are ported from reviews/phase2e_independent_checks.py (reproduction targets).
# =====================================================================================
from fractions import Fraction
from math import comb


def merged_union(ranges, tol=1e-9):
    """Union of branch intervals, overlaps merged; NOT the envelope (finding 12)."""
    out = []
    for lo, hi in sorted(r for r in ranges if r is not None):
        if out and lo <= out[-1][1] + tol:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out


def rm_branches(a, b, Mbar, base, post, pre_pairs, Lam=0.0, comp=1, extra=None):
    """As rm_range (with rate limits) but returns the list of branch intervals."""
    outs = []
    for (s0, s1) in pre_pairs:
        for sg in (1, -1):
            L, att, bvec, e, th, S = causal_lp(a, b, base, post, Lam=Lam, comp=comp, int_rows=extra)
            Dpre = L.pad(did_vec(L, 0, comp, 0, s0, s1))
            L.add_le(bvec - Mbar * sg * Dpre, 0.0)
            L.add_ge(bvec + Mbar * sg * Dpre, 0.0)
            outs.append(L.range(att))
    return outs


# Finding 1 (Prop. E18(b)): later RM horizons accumulate. Pure units, rates at elections -1, 0, 1, 2.
r_, s_, r0_ = np.array([.35, .5, .5, .5]), np.full(4, .5), np.array([.35, .5, .65, .80])
b_path = (r0_ - s_) - (r0_ - s_)[1]
A_ = abs(b_path[1] - b_path[0]); D2 = (r_[3] - r_[1]) - (s_[3] - s_[1]); att2 = r_[3] - r0_[3]
rm_ok = all(abs(b_path[t] - b_path[t - 1]) <= A_ + 1e-12 for t in (2, 3))
check("E18 [F] finding 1: at horizon 2 RM(1) holds but ATT_2 = -.30 lies outside the first-post endpoint-wise set D_2 +/- A; inside D_2 +/- 2A",
      rm_ok and np.isclose(att2, -.30) and not (D2 - A_ - 1e-12 <= att2 <= D2 + A_ + 1e-12) and (D2 - 2 * A_ - 1e-12 <= att2 <= D2 + 2 * A_ + 1e-12),
      f"A = {A_:.2f}, D_2 = {D2:.2f}")

# Finding 2 (Prop. E23(b)): finite-sample FD regression; P3.8's single-election range under-covers, the doubled range covers.
n_ = 20
rad = np.sqrt(np.log(1 / .025) / (2 * n_))
vals = 2 * np.arange(n_ + 1) / n_ - 1
wts = np.array([comb(n_, i) / 2 ** n_ for i in range(n_ + 1)])
cov1, cov2 = wts[np.abs(vals) <= rad].sum(), wts[np.abs(vals) <= 2 * rad].sum()
check("E23 [F] finding 2: fitted change = mean of 20 signed ones; P3.8 radius with one election's range covers .884682 < .95; doubled range .997423",
      round(cov1, 6) == 0.884682 and round(cov2, 6) == 0.997423, f"radius {rad:.4f}")

# Finding 3 (Prop. E25(a), zero-box clause): RM on actual quantities does not transfer to means.
atts, oks = [], True
for pre_r in (.40, .60):                       # treated pre rate; comparison .50 throughout; base .50, post .50
    Dpre = (.5 - pre_r) - 0.0                  # gap change pre -> base
    b1 = (.60 - .50) - 0.0                     # untreated treated post .60
    oks &= abs(b1) <= abs(Dpre) + 1e-12        # RM(1) at this realisation
    atts.append(.50 - .60)
mean_pre = np.mean([.5 - .40, .5 - .60])
check("E25 [F] finding 3: RM(1) holds at every realisation, mean programme effect -.10; RM on the mean pre-DiD (= 0) forces b = 0 and effect 0",
      oks and np.allclose(atts, -.10) and mean_pre == 0)

# Finding 4 (Thm E16(b)): integer witness needs (B_s - S) A_s integer as well as A_r theta.
Ar_, As_, Br_, Bs_, th_, S_4 = 1, 1, 1.0, 1.0, 0.0, 0.5
old_ok = (A_r_th := Ar_ * th_) == int(A_r_th) and 0 <= th_ <= Br_ and 0 <= Bs_ - S_4 <= 1
check("E16 [F] finding 4: theta A_r integer and continuous rate limits hold, but one comparison person cannot have untreated rate 1/2",
      old_ok and Fraction(Bs_ - S_4) * As_ != int(Fraction(Bs_ - S_4) * As_))

# Finding 5 (Remark 7.1): margins-only PT, m treated, m+1 comparison, m yes at both elections: integer {0}, continuous [-1, 1].
ok5 = True
for m in (2, 3, 5, 10):
    poss = set()
    for y0, y1 in itertools.product(range(m + 1), repeat=2):
        th = Fraction(y0, m) - Fraction(y1 - y0, m + 1)
        if 0 <= th <= 1 and (m * th).denominator == 1:
            poss.add(Fraction(y1, m) - th)
    row = [-(1 / m + 1 / (m + 1)), 1 / (m + 1), 1]; obj = np.array([0, 1 / m, -1])
    ends = [s * linprog(s * obj, A_eq=[row], b_eq=[0], bounds=[(0, m), (0, m), (0, 1)]).fun for s in (1, -1)]
    ok5 &= poss == {0} and np.allclose(ends, (-1, 1))
check("Remark 7.1 [F] finding 5: no one-step integer guarantee: integer ATT set {0}, continuous [-1, 1] (m = 2, 3, 5, 10)", ok5)

# Finding 6 (Thm E27): three pairs sharing a base: strict LP-integer gap (3/2 vs 1).
P6 = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
lp6 = -linprog(-np.ones(3), A_ub=P6, b_ub=np.ones(3), bounds=(0, 1)).fun
ip6 = max(sum(v) for v in itertools.product(range(2), repeat=3) if np.all(P6 @ v <= 1))
check("E27 [F] finding 6: degree-3 star with one change row per edge: LP maximum 3/2, integer maximum 1 (strict gap)",
      np.isclose(lp6, 1.5) and ip6 == 1)
# cycles: a 4-cycle of period pairs, each with its own laminar family, is TU (LP = MILP)
ok_cyc, n_cyc = True, 0
for rep in range(15):
    a, b, n = rand_margins(4, 2, 3, 4, 4, 9)
    L = CountLP(a, b)
    for (t0, t1) in [(0, 1), (1, 2), (2, 3), (0, 3)]:
        for U in laminar_family(4):
            v = np.zeros(L.nn)
            for i in U:
                v[L.ix(t1, i, 0, 0)] += 1; v[L.ix(t0, i, 0, 0)] -= 1
            tv = true_val(v, n)
            L.add_ge(v, tv - rng.integers(0, 2)); L.add_le(v, tv + rng.integers(0, 2))
    for k in range(2):
        obj = rng.normal(size=L.nn)
        cr, ir = L.range(obj), L.irange(obj)
        ok_cyc &= abs(cr[0] - ir[0]) < 1e-6 and abs(cr[1] - ir[1]) < 1e-6; n_cyc += 1
check("E27 [V] finding 6: a 4-cycle of period pairs (chronological orientation mixed) with laminar families: LP = MILP", ok_cyc, f"{n_cyc} objectives")

# Finding 7 (Lemma E7(d)): zero common-composition contrasts, electorate DiD loading .0075.
X7 = np.array([*np.eye(3)] + [(np.eye(3)[i] + np.eye(3)[j]) / 2 for i, j in [(0, 1), (0, 2), (1, 2), (0, 2)]])
K7 = np.array([[0, 0, 1], [0, 0, 1], [-1, -1, 0]]) / 10
mc7 = X7 @ K7.T
own7 = (X7 * mc7).sum(0) / X7.sum(0)
quad7 = np.column_stack([X7[:, i] * X7[:, j] for i in range(3) for j in range(i, 3)])
check("E7 [F] finding 7: every common-composition contrast for pair (1,2) is zero, yet the electorate DiD loads .0075 (full-rank quadratic design)",
      np.allclose(mc7[:, 0], mc7[:, 1]) and np.isclose(own7[0] - own7[1], .0075) and np.linalg.matrix_rank(quad7) == 6)

# Finding 8 (Prop. E8): pairwise FD sets do not stack across a shared election.
def path_feasible(x, incs):
    T = len(incs) + 1
    rows = np.kron(np.eye(T), np.array([x, 1 - x]))
    chg = np.kron(np.diff(np.eye(T), axis=0), np.eye(2))
    rhs = np.r_[np.full(T, .5), np.tile([1 - x, -x], len(incs)) * np.repeat(incs, 2)]
    return linprog(np.zeros(2 * T), A_eq=np.vstack([rows, chg]), b_eq=rhs, bounds=(0, 1)).success
check("E8 [F] finding 8: each of two successive change pairs is feasible at x = .4, .5, .6, the three-election path is not",
      all(path_feasible(x, [1.25]) and not path_feasible(x, [1.25, 1.25]) for x in (.4, .5, .6)))

# Finding 9 (Thm E16(d)): PT refuted jointly across two post-periods though feasible at each separately.
fe = [linprog(np.zeros(2), A_eq=[[1, 1]], b_eq=[1.2], A_ub=np.eye(2)[rs], b_ub=[.5] * len(rs), bounds=(.2, 1)).success
      for rs in ([0], [1], [0, 1])]
check("E16 [F] finding 9: comparison post rates s1 + s2 = 1.2; PT feasible at each post-period separately, refuted jointly", fe == [True, True, False])
# Cor. E17(c): the sharp breakdown value (LP in M) is at least the Minkowski value min D - Lam (worked example, PT + CR).
Lb, attb, bvb, eb, thb, Sb = causal_lp(aW, bW, 1, 2, int_rows=crP)
Mv = Lb.add_var(0.0, None)
eb = lambda j: np.eye(Lb.nvar())[j]
Lb.add_le(Lb.pad(bvb) - eb(Mv), 0.0); Lb.add_ge(Lb.pad(bvb) + eb(Mv), 0.0); Lb.add_le(Lb.pad(attb), 0.0)
Msharp = Lb.opt(eb(Mv), 1)[0]
Mmink = LWc.range(did_vec(LWc, 0, 1, 0, 1, 2))[0]
check("E17 [V] finding 9: sharp breakdown of 'ATT > 0' = min M s.t. joint feasibility and ATT <= 0 (one LP) >= Minkowski value",
      Msharp is not None and Msharp >= Mmink - 1e-9, f"sharp {Msharp:.4f}, Minkowski {Mmink:.4f}")

# Finding 10 (Prop. E19(c)): RM informative with the counts alone; two mixed precincts (99,1), (1,99).
aM = np.array([[99, 1], [1, 99]]); yM = np.array([[50, 50], [50, 50], [70, 50]])
bM = np.stack([np.column_stack([yM[t], aM.sum(1) - yM[t]]) for t in range(3)])
mixed = merged_union(rm_branches(aM, bM, 1.0, 1, 2, [(0, 1)]))
LMc, attM, *_ = causal_lp(aM, bM, 1, 2)
noaM = LMc.range(attM)
check("E19 [F] finding 10: counts alone, no pure precinct, RM(1): ATT in [.12, .28] (no-assumption [-.31, .71]); no pre-period model needed",
      len(mixed) == 1 and np.allclose(mixed[0], (.12, .28)) and np.allclose(noaM, (-.31, .71)),
      f"RM [{mixed[0][0]:.4f},{mixed[0][1]:.4f}]")

# Finding 11 (G6): clean witness. Precincts (90,10), (10,90), baseline yes 100 and 0; a 50/50 split of the merged
# baseline keeps the total but gives baseline gap [-.2, .2], excluding the true .8.
a11 = np.array([[90, 10], [10, 90]])
L11 = CountLP(a11, np.array([[[50, 50], [50, 50]]]))
g11 = L11.range(rate_vec(L11, 0, 0, 0) - rate_vec(L11, 1, 0, 0))
check("G6 [F] finding 11: a 50/50 split of a merged baseline gives gap [-.2, .2]; the true gap .8 is excluded",
      np.allclose(g11, (-.2, .2)) and not (g11[0] <= .8 <= g11[1]))

# Finding 12: an RM set with a hole although the count set is one polytope; the envelope would include -.50.
obj12, br12 = np.array([2., -1.]), []
for sg in (-1, 1):
    rows = [[-sg, 1], [-sg, -1]]
    f = [linprog(k * obj12, A_ub=rows, b_ub=[1.1, -1.1], bounds=[(-.2, .2), (0, 1)]) for k in (1, -1)]
    br12.append((.5 + f[0].fun, .5 - f[1].fun))
mu12 = merged_union(br12); env12 = union_range(br12)
check("E18 [F] finding 12: sharp RM projection [-.90,-.70] U [-.30,0] with one polytope of tables; the envelope wrongly contains -.50",
      len(mu12) == 2 and np.allclose(mu12, [(-.9, -.7), (-.3, 0)]) and env12[0] <= -.5 <= env12[1])
check("E18 [V] finding 12: N21 and the worked-example RM sets are single intervals after merging branch intervals (no hole)",
      all(len(merged_union(rm_branches(*args))) == 1 for args in
          [(a21, b21, 1.0, 1, 2, [(0, 1)]), (aW, bW, 1.0, 1, 2, [(0, 1)], 0.0, 1, crP + crQ)]))
# E3(c): consecutive differences with only election 0 uncertain and elections 1, 2 point-identified: a box.
G3 = [(-.2, .3), (.1, .1), (.25, .25)]
corners = [(g1 - g0, g2 - g1) for g0 in G3[0] for g1 in G3[1] for g2 in G3[2]]
check("E3 [V] finding 12: consecutive DiDs are a box when no interior election is uncertain; base-relative (base 0) is a segment",
      len({c[1] for c in corners}) == 1 and len({(g1 - g0, g2 - g0) for g0 in G3[0] for g1 in G3[1] for g2 in G3[2]}) == 2)


# =====================================================================================
# PT. Within-precinct parallel trends in expectation (Props E28-E30); ported from
#     reviews/phase2e_independent_checks_pt.py. Expected rates are conditional means given the design;
#     a precinct's aggregate law is a finite list of equally likely totals, and the feasible set of its
#     expected group rates is the Aumann expectation (per-state Duncan-Davis ends averaged).
# =====================================================================================
def mean_pt_lp(size, vets, pre, post, comparison, lce=("change",), eps=0.0, eps_var=False, Lam=0.0):
    vets = np.asarray(vets, float); size = np.broadcast_to(np.asarray(size, float), vets.shape)
    ctl = size - vets; x = vets / size; p = len(x)
    pre = np.asarray(pre, float).reshape(p, -1); post = np.asarray(post, float).reshape(p, -1)
    sup = ctl > 0
    # variables: v0, c0, v1, c1 (p each), theta_U (p), affine params 2 per (group, period) x 6, eps
    nv = 5 * p + 12 + 1 + p                                       # + spillover S_i (WPT sensitivity)
    V0, C0, V1, C1, TH = (np.arange(p) + k * p for k in range(5))
    par = lambda g, per: 5 * p + 2 * (2 * per + g)                 # g 0 = V, 1 = C; per 0 pre, 1 post, 2 change
    EPS = 5 * p + 12
    SS = 5 * p + 13 + np.arange(p)
    Aeq, beq, Aub, bub = [], [], [], []
    bnds = [(0, 1)] * (4 * p) + [(0, 1) if not sup[i] else (0, 0) for i in range(p)] + [(None, None)] * 12
    bnds += [(0, None) if eps_var else (eps, eps)]
    bnds += [(-Lam, Lam) if (sup[i] and comparison == "within") else (0, 0) for i in range(p)]
    for (tot, Vt, Ct) in ((pre, V0, C0), (post, V1, C1)):
        lo = np.maximum(0, tot - ctl[:, None]).mean(1) / vets
        hi = np.minimum(vets[:, None], tot).mean(1) / vets
        for i in range(p):
            bnds[Vt[i]] = (lo[i], hi[i])
            if not sup[i]:
                bnds[Ct[i]] = (0, 0)
            v = np.zeros(nv); v[Vt[i]] = x[i]; v[Ct[i]] = 1 - x[i]; Aeq.append(v); beq.append(tot[i].mean() / size[i])
    for per in lce:
        k = {"pre": 0, "post": 1, "change": 2}[per]
        for g, (G0, G1) in enumerate(((V0, V1), (C0, C1))):
            for i in range(p):
                if g == 1 and not sup[i]:
                    continue
                v = np.zeros(nv)
                if k == 0: v[G0[i]] = 1
                elif k == 1: v[G1[i]] = 1
                else: v[G1[i]] = 1; v[G0[i]] = -1
                v[par(g, k)] = -1; v[par(g, k) + 1] = -x[i]
                if k == 0 and (eps > 0 or eps_var):            # baseline curvature allowance |dev| <= eps
                    w = v.copy(); w[EPS] = -1; Aub.append(w); bub.append(0)
                    w = -v.copy(); w[EPS] = -1; Aub.append(w); bub.append(0)
                else:
                    Aeq.append(v); beq.append(0)
    wv = vets / vets.sum(); wc = ctl / ctl.sum()
    obj = np.zeros(nv)
    if comparison == "within":
        for i in range(p):
            if sup[i]:
                obj[V1[i]] += wv[i]; obj[V0[i]] -= wv[i]; obj[C1[i]] -= wv[i]; obj[C0[i]] += wv[i]; obj[SS[i]] += wv[i]
                v = np.zeros(nv); v[V0[i]] = 1; v[C1[i]] = 1; v[C0[i]] = -1; v[SS[i]] = -1   # 0 <= m_V0 + d_N - S_i <= 1
                Aub += [v, -v]; bub += [1, 0]
                v = np.zeros(nv); v[C1[i]] = 1; v[SS[i]] = -1                         # 0 <= m_Nt - S_i <= 1
                Aub += [v, -v]; bub += [1, 0]
            else:
                obj[V1[i]] += wv[i]; obj[TH[i]] -= wv[i]
    else:
        obj[V1] += wv; obj[V0] -= wv; obj[C1] -= wc; obj[C0] += wc
        v = np.zeros(nv); v[V0] = wv; v[C1] = wc; v[C0] -= wc
        Aub += [v, -v]; bub += [1, 0]
    return dict(obj=obj, Aeq=np.array(Aeq), beq=np.array(beq), Aub=np.array(Aub) if Aub else None,
                bub=np.array(bub) if bub else None, bnds=bnds, EPS=EPS, nv=nv, sup=sup, wv=wv, V1=V1, TH=TH)


def lp_range(P, obj=None):
    obj = P["obj"] if obj is None else obj
    out = []
    for s in (1, -1):
        r = linprog(s * obj, A_eq=P["Aeq"], b_eq=P["beq"], A_ub=P["Aub"], b_ub=P["bub"], bounds=P["bnds"], method="highs")
        if r.status != 0:
            return None
        out.append(s * r.fun)
    return tuple(100 * np.array(out))


rare = (20000, [5000, 8000, 11000, 14000, 17000], [40, 12000, 12000, 12000, 12000],
        [[50, 30], [9964, 9956], [8188, 8192], [6682, 6698], [5446, 5474]])
ordinary = (1000, [100, 300, 500, 700, 900], [500] * 5, [[682, 650], [606, 590], [530, 530], [454, 470], [378, 410]])
PT5 = {}
for nm, data in (("rare", rare), ("ordinary", ordinary)):
    for cmp in ("within", "pooled"):
        PT5[(nm, cmp)] = lp_range(mean_pt_lp(*data, comparison=cmp))
for k, v in PT5.items():
    print(f"      PT  five precincts, {k[0]}, {k[1]} mean PT: [{v[0]:+.2f}, {v[1]:+.2f}]")
exp5 = {("rare", "within"): (10.39, 12.25), ("rare", "pooled"): (-3.96, -2.22),
        ("ordinary", "within"): (-67.33, 43.78), ("ordinary", "pooled"): (-80.56, 52.13)}
check("E28 [V] five-precinct table reproduced (FD-LCE on changes): rare [+10.39,+12.25] within vs [-3.96,-2.22] pooled; ordinary [-67.33,+43.78] vs [-80.56,+52.13]",
      all(np.allclose(PT5[k], exp5[k], atol=0.0051) for k in exp5))

three = (800, [200, 400, 600], [4, 164, 324], [[204], [124], [4]])
exp3 = {(): (-58.00, 27.33), ("change",): (-49.33, 7.33), ("pre", "change"): (2.67, 7.33), ("pre", "post"): (2.67, 7.33),
        ("post", "change"): (2.67, 7.33), ("pre", "post", "change"): (2.67, 7.33)}
got3 = {k: lp_range(mean_pt_lp(*three, comparison="within", lce=k)) for k in exp3}
for k, v in got3.items():
    print(f"      PT  three precincts, LCE on {' + '.join(k) or 'nothing'}: [{v[0]:+.2f}, {v[1]:+.2f}]")
check("E30 [V] three-precinct check: no LCE [-58.00,27.33]; change only [-49.33,7.33]; any two of pre/post/change, or all three, [2.67,7.33]",
      all(np.allclose(got3[k], exp3[k], atol=0.0051) for k in exp3))
# baseline-curvature allowance eps: monotone, eps = 0 recovers two-election LCE, large eps recovers change-only;
# breakdown eps* = min eps s.t. ATT <= 0 is one LP (eps a variable).
epsr = [lp_range(mean_pt_lp(*three, comparison="within", lce=("pre", "change"), eps=e)) for e in (0.0, .001, 1 / 600, .005, .02, .05, .5, 1.0)]
mono = all(epsr[j + 1][0] <= epsr[j][0] + 1e-7 and epsr[j + 1][1] >= epsr[j][1] - 1e-7 for j in range(len(epsr) - 1))
Pb = mean_pt_lp(*three, comparison="within", lce=("pre", "change"), eps_var=True)
cE = np.zeros(Pb["nv"]); cE[Pb["EPS"]] = 1
rb = linprog(cE, A_eq=Pb["Aeq"], b_eq=Pb["beq"], A_ub=np.vstack([Pb["Aub"], Pb["obj"]]), b_ub=np.r_[Pb["bub"], 0.0],
             bounds=Pb["bnds"], method="highs")
eps_star = rb.fun
at_star = lp_range(mean_pt_lp(*three, comparison="within", lce=("pre", "change"), eps=eps_star))
check("E30 [I] baseline-curvature allowance: monotone in eps; positive for eps < 1/600, zero feasible at eps* = 1/600 (one LP); eps = 1/2 already gives change-only",
      mono and np.allclose(epsr[0], exp3[("pre", "change")], atol=.0051) and np.allclose(epsr[-1], exp3[("change",)], atol=.0051) and np.allclose(epsr[-2], epsr[-1]) and epsr[1][0] > 1e-6
      and abs(at_star[0]) < 1e-6, f"eps* = {100 * eps_star:.4f} points; at eps* [{at_star[0]:.2f},{at_star[1]:.2f}]")
# unsupported precincts (no comparison voters): full veteran target, no-assumption bound on that part.
four = (800, [200, 400, 600, 800], [4, 164, 324, 720], [[204], [124], [4], [80]])   # precinct 4: veterans only
P4 = mean_pt_lp(*four, comparison="within", lce=("change",))
full4 = lp_range(P4)
objS = P4["obj"].copy(); objS[P4["V1"][3]] = 0; objS[P4["TH"][3]] = 0
supS = lp_range(P4, objS)
piU = P4["wv"][3]; cU = 100 * piU * 80 / 800
check("E28 [V] unsupported precinct (all veterans): ATT interval = supported part + w_V B_V1 - [pi_U, 0], pi_U = share of veterans affected",
      np.allclose(full4, (supS[0] + cU - 100 * piU, supS[1] + cU), atol=1e-6), f"pi_U = {piU:.3f}; [{full4[0]:.2f},{full4[1]:.2f}]")

# The four-precinct example: realised descriptive targets D_W (within, veteran weights) and the pooled DiD, margins only.
LD = CountLP(aW, bW)
wvW = aW[:, 0] / aW[:, 0].sum(); wnW = aW[:, 1] / aW[:, 1].sum()
DW = np.zeros(LD.nn)
for i in range(4):
    for r, sgn in ((0, 1), (1, -1)):
        DW[LD.ix(2, i, r, 0)] += sgn * wvW[i] / aW[i, r]; DW[LD.ix(1, i, r, 0)] -= sgn * wvW[i] / aW[i, r]
DWr = 100 * np.array(LD.range(DW)); DPr = 100 * np.array(LD.range(did_vec(LD, 0, 1, 0, 1, 2)))
check("E28 [I] four-precinct example, margins only: D_W in [-93.87, +144.72], pooled DiD in [-40.75, +79.25]; precinct 4 weight 47.5% vs 2.5%",
      np.allclose(DWr, (-93.87, 144.72), atol=.0051) and np.allclose(DPr, (-40.75, 79.25), atol=.0051)
      and np.isclose(wvW[3], .475) and np.isclose(wnW[3], .025),
      f"D_W [{DWr[0]:.2f},{DWr[1]:.2f}]")

# Prop. E29: [L - r, U + r] covers the expected ATT; the realised-table value alone does not.
pS = 30
NS_ = rng.integers(60, 160, pS); xS = rng.uniform(.2, .8, pS); AV = np.round(NS_ * xS).astype(int); AN = NS_ - AV
pV0 = rng.uniform(.2, .6, pS); pN0 = rng.uniform(.2, .6, pS); trend = rng.uniform(-.1, .1, pS); tau = .05
pV1, pN1 = pV0 + trend + tau, pN0 + trend
wS = AV / AV.sum()
att_mean = np.sum(wS * ((pV1 - pV0) - (pN1 - pN0)))
rS = np.sqrt(.5 * np.log(2 / .05) * np.sum((4 * wS) ** 2))
neffS = 1 / np.sum(wS ** 2)                        # r = sqrt(8 log(2/alpha) / n_eff); 0.992 for 30 equal weights
rP = np.sqrt(.5 * np.log(2 / .05) * np.sum(4 * wS ** 2 * (1 / AV + 1 / AN)))   # persons independent: range 2 w/A per person
cov_pt, cov_int, cov_p = 0, 0, 0
for rep in range(400):
    yV = [rng.binomial(AV, q) for q in (pV0, pV1)]; yN = [rng.binomial(AN, q) for q in (pN0, pN1)]
    dw = np.sum(wS * ((yV[1] - yV[0]) / AV - (yN[1] - yN[0]) / AN))       # all cells audited: L = U = D_W
    cov_pt += abs(dw - att_mean) < 1e-12
    cov_int += dw - rS <= att_mean <= dw + rS
    cov_p += dw - rP <= att_mean <= dw + rP
check("E29 [F] realised-table D_W (cells known, L = U) is not a confidence set for the expected ATT; [L - r, U + r] covers (Hoeffding: rho_i = 4 w_Vi, or person-level)",
      cov_pt / 400 < .05 and cov_int / 400 >= .95 and cov_p / 400 >= .95
      and np.isclose(rS, np.sqrt(8 * np.log(2 / .05) / neffS)) and round(np.sqrt(8 * np.log(40) / 30), 3) == 0.992,
      f"point {cov_pt / 400:.3f}; n_eff = {neffS:.1f} (unequal weights); any within-precinct dependence r = {rS:.3f} (equal weights: 0.992), cover {cov_int / 400:.3f}; independent persons r = {rP:.3f}, cover {cov_p / 400:.3f}")


# =====================================================================================
# R2. Checks added after the independent review's round 2 on sections 5.8-5.9
#     (reviews/PHASE2e_independent_review_round2_WPT.md).
# =====================================================================================
# A1 (E28(c), PPT witness): veteran baselines (.1, .9), comparison changes (-.3, +.3), equal weights.
m0r, dNr, wr = np.array([.1, .9]), np.array([-.3, .3]), np.array([.5, .5])
local = m0r + dNr; p0 = wr @ m0r + wr @ dNr
check("E28 [F] round 2: under PPT the local imputation m_V0 + d_N = (-.2, 1.2) is not a rate; the pooled witness p0 = .5 in every precinct is valid",
      np.allclose(local, (-.2, 1.2)) and np.isclose(p0, .5) and 0 <= p0 <= 1 and not np.all((local >= 0) & (local <= 1)))
# A2 (E29(d), RM on increments): pre maximum .10, post gaps .10, .20 (relative to the base).
bpost = np.array([0, .10, .20]); preA = .10
inc_ok = all(abs(bpost[t] - bpost[t - 1]) <= preA + 1e-12 for t in (1, 2))
lvl_ok = all(abs(bpost[t]) <= preA + 1e-12 for t in (1, 2))
check("E29 [F] round 2: RM on increments |b_t - b_(t-1)| <= Mbar A holds (Mbar = 1), the cap on levels |b_t| <= Mbar A fails at the second post-election",
      inc_ok and not lvl_ok)
# A4 (spillover sensitivity under WPT): joint LP with S_i vs dilation of the zero-spillover interval by pi_S Lam.
for LamS in (.01, .05):
    jt = lp_range(mean_pt_lp(*three, comparison="within", lce=("change",), Lam=LamS))
    z0 = lp_range(mean_pt_lp(*three, comparison="within", lce=("change",)))
    dil = (z0[0] - 100 * LamS, z0[1] + 100 * LamS)          # pi_S = 1 here
    print(f"      R2  spillover {LamS}: joint LP [{jt[0]:.4f},{jt[1]:.4f}] vs dilation [{dil[0]:.4f},{dil[1]:.4f}]")
spill_rows = [(lp_range(mean_pt_lp(*three, comparison="within", lce=("change",), Lam=L_)),
               lp_range(mean_pt_lp(*three, comparison="within", lce=("change",)))) for L_ in (.01, .05)]
check("E28 [I] round 2: WPT spillover sensitivity by one joint LP in (fit, S_i) with local rate limits; reported beside the pi_S*Lam dilation (no sharpness claimed for the dilation)",
      all(j is not None and j[0] <= z[0] + 1e-9 and j[1] >= z[1] - 1e-9 for j, z in spill_rows))

# B10 (Prop. E34): one pre transition, two-group FD-LCE: pre increment q(x_i) = alpha + beta x_i (affine);
# exact local RM |db_i| <= Mbar |q(x_i)|. Split-and-orientation enumeration (<= 2(|S|+1) LPs) equals
# exhaustive enumeration of 2^|S| sign patterns on both extrema.
def rm_affine(xs, lo, hi, base, w, Mbar, patterns):
    pS = len(xs); nv_ = 2 + pS                      # alpha, beta, db_1..db_p
    out = []
    for sg in patterns:
        A, b = [], []
        for i in range(pS):
            q = np.zeros(nv_); q[0] = 1; q[1] = xs[i]
            A += [q, -q]; b += [hi[i], -lo[i]]                     # data: q(x_i) in [lo_i, hi_i]
            A.append(-sg[i] * q); b.append(0.0)                    # sign pattern sg_i q(x_i) >= 0
            e = np.zeros(nv_); e[2 + i] = 1
            A += [e - Mbar * sg[i] * q, -e - Mbar * sg[i] * q]; b += [0.0, 0.0]   # |db_i| <= Mbar sg_i q_i
            A += [e, -e]; b += [1 - base[i], base[i]]              # local rate limit: base_i + db_i in [0, 1]
        obj = np.zeros(nv_); obj[2:] = -w                          # ATT part: -sum w_i db_i
        r1 = linprog(obj, A_ub=A, b_ub=b, bounds=[(None, None)] * nv_, method="highs")
        r2 = linprog(-obj, A_ub=A, b_ub=b, bounds=[(None, None)] * nv_, method="highs")
        if r1.status == 0 and r2.status == 0:
            out.append((r1.fun, -r2.fun))
    return (min(o[0] for o in out), max(o[1] for o in out)) if out else None
ok34, n34 = True, 0
for rep in range(10):
    pS = 6
    xs = np.sort(rng.uniform(.1, .9, pS)); a_, b_ = rng.uniform(-.2, .2), rng.uniform(-.4, .4)
    tru = a_ + b_ * xs; sl = rng.uniform(.01, .08, pS)
    lo, hi = tru - sl, tru + sl
    base = rng.uniform(.2, .8, pS); w = rng.dirichlet(np.ones(pS)); Mbar = rng.uniform(.5, 2)
    split = [tuple(s_ * (1 if i < k else -1) for i in range(pS)) for k in range(pS + 1) for s_ in (1, -1)]
    allp = list(itertools.product((1, -1), repeat=pS))
    e1, e2 = rm_affine(xs, lo, hi, base, w, Mbar, split), rm_affine(xs, lo, hi, base, w, Mbar, allp)
    ok34 &= e1 is not None and np.allclose(e1, e2, atol=1e-8); n34 += 1
check("E34 [V] round 2: exact within-precinct RM with an affine pre-increment profile: split-and-orientation enumeration (<= 2(|S|+1) LPs) = exhaustive 2^|S| sign enumeration",
      ok34, f"{n34} six-precinct cases, both extrema")


# =====================================================================================
# E35. Realised within-precinct parallel trends with tolerance tau_i = M sqrt(1/a_iV + 1/a_iN)
#      (owner's default reading of WPT). Worked example: the four precincts of section 9.2, base -> post.
# =====================================================================================
def rwpt_lp(a, b, t0, t1, M=0.0, Mvar=False, Lam=0.0, rows=None, amin=1, tau=None):
    """E35: realised within-precinct PT with tolerance tau_i = M sqrt(1/a_iV + 1/a_iN) (or given tau).
    Variables: n, theta_i (realised untreated veteran rate at t1, every precinct with veterans),
    S_i (realised spillover on local non-veterans), M (optional variable). Treated row 0, comparison row 1."""
    L = CountLP(a, b)
    for coef, rhs in (rows or []):
        L.add_ge(coef, rhs)
    p = L.p; A = L.a
    wV = A[:, 0] / A[:, 0].sum()
    sup = A[:, 1] >= amin
    th = [L.add_var(0.0, 1.0) for i in range(p)]
    Sv = [L.add_var(-Lam, Lam) if sup[i] else L.add_var(0.0, 0.0) for i in range(p)]
    Mv = L.add_var(0.0, None) if Mvar else None
    e = lambda j: np.eye(L.nvar())[j]
    att = np.zeros(L.nvar())
    for i in range(p):
        if A[i, 0] == 0:
            continue
        bV = lambda t: L.pad(np.eye(L.nn)[L.ix(t, i, 0, 0)] / A[i, 0])
        att += wV[i] * (bV(t1) - e(th[i]))
        if sup[i]:
            bN = lambda t: L.pad(np.eye(L.nn)[L.ix(t, i, 1, 0)] / A[i, 1])
            gap = e(th[i]) - bV(t0) - (bN(t1) - e(Sv[i]) - bN(t0))     # untreated V change - untreated N change
            c = np.sqrt(1 / A[i, 0] + 1 / A[i, 1]) if tau is None else tau[i]
            if Mvar:
                L.add_le(gap - c * e(Mv), 0.0); L.add_ge(gap + c * e(Mv), 0.0)
            else:
                L.add_le(gap, c * M); L.add_ge(gap, -c * M)
            L.add_ge(bN(t1) - e(Sv[i]), 0.0); L.add_le(bN(t1) - e(Sv[i]), 1.0)
    return L, att, Mv, e


def audit_rows(L, cells):
    out = []
    for (t, i, r, val) in cells:
        v = np.zeros(L.nn); v[L.ix(t, i, r, 0)] = 1
        out += [(v, val), (-v, -val)]
    return out

L35 = CountLP(aW, bW)
aud35 = audit_rows(L35, [(1, 2, 1, 20), (2, 2, 1, 22), (1, 3, 1, 4), (2, 3, 1, 5)])   # hypothetical audited non-veteran cells
Lm, attm, _, _ = rwpt_lp(aW, bW, 1, 2, 0.0); mo35 = Lm.range(attm)
r35 = {}
for M in (0.0, 0.1, 0.2):
    Lx, ax, _, _ = rwpt_lp(aW, bW, 1, 2, M, rows=aud35); r35[M] = Lx.range(ax)
Lb35, ab35, Mv35, eb35 = rwpt_lp(aW, bW, 1, 2, Mvar=True, rows=aud35); Lb35.add_le(ab35, 0.0)
Mstar35 = Lb35.opt(eb35(Mv35), 1)[0]
Ls35, as35, _, _ = rwpt_lp(aW, bW, 1, 2, Mstar35, rows=aud35); at_star35 = Ls35.range(as35)
for k, v in r35.items():
    print(f"      E35  audited cells, M = {k}: [{v[0]:.4f}, {v[1]:.4f}]")
print(f"      E35  margins only, M = 0: [{mo35[0]:.4f}, {mo35[1]:.4f}];  breakdown M* = {Mstar35:.4f}")
check("E35 [V] realised WPT: intervals nested in M; breakdown M* of 'ATT > 0' by one LP (lower end 0 at M*, positive below)",
      r35[0.0][0] >= r35[0.1][0] >= r35[0.2][0] > 0 and r35[0.0][1] <= r35[0.1][1] <= r35[0.2][1]
      and abs(at_star35[0]) < 1e-7 and Mstar35 > 0.2, f"M* = {Mstar35:.4f}")
# M = 0, no spillover, full support: realised ATT = D_W, so the causal set lies in D_W's count set (equal if rate limits slack)
LD35 = CountLP(aW, bW)
for rr in aud35:
    LD35.add_ge(*rr)
DW35 = np.zeros(LD35.nn)
for i in range(4):
    for r, sg in ((0, 1), (1, -1)):
        DW35[LD35.ix(2, i, r, 0)] += sg * wvW[i] / aW[i, r]; DW35[LD35.ix(1, i, r, 0)] -= sg * wvW[i] / aW[i, r]
dw35 = LD35.range(DW35)
check("E35 [V] at M = 0 (exact realised WPT, no spillover, all precincts supported) the causal set lies in the count set of D_W; equal when the rate limits are slack",
      r35[0.0][0] >= dw35[0] - 1e-9 and r35[0.0][1] <= dw35[1] + 1e-9, f"ATT [{r35[0.0][0]:.4f},{r35[0.0][1]:.4f}], D_W [{dw35[0]:.4f},{dw35[1]:.4f}]")
check("E35 [I] margins only, realised WPT does not sign the effect even at M = 0 (each precinct's non-veteran change is free): close to the no-assumption set [-0.5675, 0.7325]",
      mo35[0] < 0 < mo35[1] and abs(mo35[0] + 0.5675) < 1e-6, f"[{mo35[0]:.4f},{mo35[1]:.4f}]")

# ---------------------------------------------------------------- E36
print("\n[E36] baseline-curvature allowance on the realised route (E35 + E13)")


def add_lce_eps(L, x, wV, bhat, h, base="eps", eps=0.0, epsv=None, link=True):
    """E36 model layer, two groups (row 0 = V, 1 = N): change LCE exact (FD), base level means m_g0(x_i) within
    eps of a_g + b_g x_i; beta = (gamma: base-level quadratic, delta: change quadratic) in the box bhat +- h
    (a polyhedral projection region). base: 'eps' (rows with eps or variable epsv) or 'none' (no baseline rows).
    link: the base-level regression identity gamma = P (x m_V0 + (1 - x) m_N0). Returns a function giving the
    model-implied D_W coefficient vector."""
    K = len(x); Z = np.column_stack([np.ones(K), x, x ** 2]); P = np.linalg.pinv(Z)
    m = [[L.add_var(0.0, 1.0) for i in range(K)] for g in range(2)]
    ab = [(L.add_var(), L.add_var()) for g in range(2)]
    ce = [(L.add_var(), L.add_var()) for g in range(2)]
    gam = [L.add_var(bhat[j] - h, bhat[j] + h) for j in range(3)]
    dl = [L.add_var(bhat[3 + j] - h, bhat[3 + j] + h) for j in range(3)]

    def v(d):
        out = np.zeros(L.nvar())
        for k, c in d.items():
            out[k] += c
        return out
    for g in range(2):
        for i in range(K):
            dev = {m[g][i]: 1.0, ab[g][0]: -1.0, ab[g][1]: -x[i]}
            if base == "eps":
                if epsv is None:
                    L.add_le(v(dev), eps); L.add_ge(v(dev), -eps)
                else:
                    L.add_le(v({**dev, epsv: -1.0}), 0.0); L.add_ge(v({**dev, epsv: 1.0}), 0.0)
            L.add_ge(v({m[g][i]: 1.0, ce[g][0]: 1.0, ce[g][1]: x[i]}), 0.0)
            L.add_le(v({m[g][i]: 1.0, ce[g][0]: 1.0, ce[g][1]: x[i]}), 1.0)
    if link:
        for j in range(3):
            d = {gam[j]: -1.0}
            for i in range(K):
                d[m[0][i]] = d.get(m[0][i], 0) + P[j, i] * x[i]; d[m[1][i]] = d.get(m[1][i], 0) + P[j, i] * (1 - x[i])
            L.add_eq(v(d), 0.0)
    (cV, eV), (cN, eN) = ce
    L.add_eq(v({dl[0]: 1.0, cN: -1.0}), 0.0)
    L.add_eq(v({dl[1]: 1.0, cV: -1.0, cN: 1.0, eN: -1.0}), 0.0)
    L.add_eq(v({dl[2]: 1.0, eV: -1.0, eN: 1.0}), 0.0)
    tgt = {cV: float(np.sum(wV)), eV: float(wV @ x), cN: -float(np.sum(wV)), eN: -float(wV @ x)}
    return lambda: L.pad(v(tgt))


x36 = np.array([0.25, 0.5, 0.75]); wV36 = x36 / x36.sum()
Z36 = np.column_stack([np.ones(3), x36, x36 ** 2])
T0_36 = np.array([0.005, 0.205, 0.405]); T1_36 = np.array([0.255, 0.155, 0.005])
bh36 = np.concatenate([np.linalg.solve(Z36, T0_36), np.linalg.solve(Z36, T1_36 - T0_36)])


def ci36(eps, h=0.0, base="eps", link=True):
    L = CountLP(np.array([[1, 0]]), np.array([[[1, 0]]]))
    tg = add_lce_eps(L, x36, wV36, bh36, h, base, eps, link=link)
    return L.range(tg())


a36 = np.array([[50, 150], [100, 100], [150, 50]])
yes36 = [np.array([1, 41, 81]), np.array([51, 31, 1])]
b36 = np.stack([np.column_stack([yes36[t], a36.sum(1) - yes36[t]]) for t in range(2)])


def dw_vec(L, a, t0, t1):
    w = a[:, 0] / a[:, 0].sum(); D = np.zeros(L.nn)
    for i in range(len(a)):
        for r, sg in ((0, 1), (1, -1)):
            D[L.ix(t1, i, r, 0)] += sg * w[i] / a[i, r]; D[L.ix(t0, i, r, 0)] -= sg * w[i] / a[i, r]
    return D


def comb36(eps, h=0.0, M=0.0):
    ci = ci36(eps, h)
    L, att, _, _ = rwpt_lp(a36, b36, 0, 1, M)
    D = dw_vec(L, a36, 0, 1); L.add_ge(D, ci[0]); L.add_le(D, ci[1])
    return ci, L.range(att)


def estar36_joint(h=0.0, M=0.0):
    """One LP: minimise eps s.t. E35 rows, two model witnesses (lower, upper end of CI_LCE(eps)) sharing eps,
    D_W(n) between them, ATT <= 0. Valid with a polyhedral projection region (here a box)."""
    L, att, _, _ = rwpt_lp(a36, b36, 0, 1, M)
    ev = L.add_var(0.0, None)
    t_lo = add_lce_eps(L, x36, wV36, bh36, h, epsv=ev); t_hi = add_lce_eps(L, x36, wV36, bh36, h, epsv=ev)
    D = L.pad(dw_vec(L, a36, 0, 1))
    L.add_le(t_lo() - D, 0.0); L.add_ge(t_hi() - D, 0.0); L.add_le(L.pad(att), 0.0)
    return L.opt(np.eye(L.nvar())[ev], 1)[0]


def estar36_bisect(h=0.0, M=0.0, lo=0.0, hi=0.5, it=60):
    for _ in range(it):
        mid = (lo + hi) / 2
        if comb36(mid, h, M)[1][0] > 1e-12:
            lo = mid
        else:
            hi = mid
    return hi


res36 = {e: comb36(e) for e in (0.0, 0.001, 0.01, 0.5)}
for e, (ci, r) in res36.items():
    print(f"      E36  eps = {e}: CI_LCE [{ci[0]:.4f}, {ci[1]:.4f}];  realised ATT (E35, M = 0) [{r[0]:.4f}, {r[1]:.4f}]")
ci_none = ci36(0.5, base="none"); ci_fd = ci36(0.5, base="none", link=False)
es_j = estar36_joint(); es_b = estar36_bisect()
es_jh = estar36_joint(h=0.0002); es_bh = estar36_bisect(h=0.0002)
print(f"      E36  change-only with base-level link {ci_none};  without link (fd) {ci_fd}")
print(f"      E36  eps* joint LP {es_j:.6f}, bisection {es_b:.6f};  box h = 0.0002: {es_jh:.6f}, {es_bh:.6f}")
check("E36 [V] nested in eps; eps = 1/2 equals change-only with the base-level link and lies inside fd() without it",
      all(res36[e1][0][0] >= res36[e2][0][0] - 1e-9 and res36[e1][0][1] <= res36[e2][0][1] + 1e-9
          and res36[e1][1][0] >= res36[e2][1][0] - 1e-9 and res36[e1][1][1] <= res36[e2][1][1] + 1e-9
          for e1, e2 in ((0.0, 0.001), (0.001, 0.01), (0.01, 0.5)))
      and np.allclose(res36[0.5][0], ci_none, atol=1e-8) and ci_fd[0] <= ci_none[0] + 1e-9 and ci_fd[1] >= ci_none[1] - 1e-9,
      f"eps=1/2 {res36[0.5][0]}, change-only+link {ci_none}, fd {ci_fd}")
check("E36 [V] breakdown eps* of 'ATT > 0' on the realised route: one LP (polyhedral projection) equals monotone bisection",
      res36[0.0][1][0] > 0 and abs(es_j - es_b) < 1e-7 and abs(es_jh - es_bh) < 1e-7 and es_jh < es_j,
      f"eps* {es_j:.6f} vs {es_b:.6f}; h=0.0002: {es_jh:.6f} vs {es_bh:.6f}")

print(f"\n{sum(results)}/{len(results)} PASS")
