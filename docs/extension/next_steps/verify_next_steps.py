"""Checks for NEXT_STEPS_RECOMMENDATION.md (next-steps planning, 20 Nov 2026).

Each check prints PASS/FAIL with the numbers it relies on. Labels N1..N12 are cited
in the recommendation. Run:  python3 -I verify_next_steps.py
"""
import itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

rng = np.random.default_rng(20261120)
TOL = 1e-8
results = []


def report(label, ok, msg):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'} {label}: {msg}")


def lp_max(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None):
    r = linprog(-np.asarray(c, float), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                bounds=bounds, method="highs")
    assert r.status == 0, r.message
    return -r.fun, r.x


def proper_subsets(C):
    for k in range(1, C):
        for A in itertools.combinations(range(C), k):
            yield A


# ---------------------------------------------------------------------------
# Unit sets.  2 x C: y_c = group-1 mass in column c; S(x,t) = {0<=y<=t, sum y = x}.
#             R x 2: y_r = mass of row r in the positive column; S = {0<=y<=x, sum y = t}.
# ---------------------------------------------------------------------------

def h_unit_2xC(u, x, t):
    """support function of S(x,t) for 2 x C (greedy: fill best columns first)"""
    order = np.argsort(-u)
    rem, val = x, 0.0
    for c in order:
        a = min(rem, t[c]); val += u[c] * a; rem -= a
    return val


def h_unit_Rx2(u, x, t):
    order = np.argsort(-u)
    rem, val = t, 0.0
    for r in order:
        a = min(rem, x[r]); val += u[r] * a; rem -= a
    return val


def h_unit_lp_2xC(u, x, t):
    C = len(t)
    v, _ = lp_max(u, A_eq=np.ones((1, C)), b_eq=[x], bounds=[(0, tc) for tc in t])
    return v


def h_unit_lp_Rx2(u, x, t):
    R = len(x)
    v, _ = lp_max(u, A_eq=np.ones((1, R)), b_eq=[t], bounds=[(0, xr) for xr in x])
    return v


def h_subset_polytope(u, rhs_fn, total, n, subsets):
    """max u.y s.t. y(A) <= rhs(A) for the given subsets, y(all) = total, y >= 0"""
    A_ub, b_ub = [], []
    for A in subsets:
        row = np.zeros(n); row[list(A)] = 1; A_ub.append(row); b_ub.append(rhs_fn(A))
    v, _ = lp_max(u, A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.ones((1, n)), b_eq=[total],
                  bounds=[(0, None)] * n)
    return v


# N1: subset-inequality identity, 2 x C, against the lifted LP (one table per atom)
maxdiff, ncase = 0.0, 0
for C in [2, 3, 4, 5]:
    for _ in range(40):
        J = rng.integers(2, 6)
        x = rng.uniform(0.1, 0.9)
        ts = rng.dirichlet(np.ones(C) * 0.7, size=J)
        p = rng.dirichlet(np.ones(J))
        rhs = lambda A: sum(p[j] * min(x, ts[j][list(A)].sum()) for j in range(J))
        for _ in range(5):
            u = rng.normal(size=C)
            # exact: Aumann expectation = sum_j p_j h_{S(x,t_j)}; also via the lifted LP
            h_exact = sum(p[j] * h_unit_lp_2xC(u, x, ts[j]) for j in range(J))
            h_greedy = sum(p[j] * h_unit_2xC(u, x, ts[j]) for j in range(J))
            h_sub = h_subset_polytope(u, rhs, x, C, list(proper_subsets(C)))
            maxdiff = max(maxdiff, abs(h_exact - h_sub), abs(h_exact - h_greedy))
            ncase += 1
report("N1 subset identity 2xC", maxdiff < 1e-7,
       f"{ncase} support-function comparisons (C=2..5, 2-5 atoms), max |diff| = {maxdiff:.1e}")

# N1b: lifted LP (explicit joint variables per atom) equals the sum of unit supports
maxdiff = 0.0
for _ in range(60):
    C = rng.integers(2, 5); J = rng.integers(2, 5); x = rng.uniform(.1, .9)
    ts = rng.dirichlet(np.ones(C), size=J); p = rng.dirichlet(np.ones(J)); u = rng.normal(size=C)
    # variables b_j (J*C); objective u . sum_j p_j b_j
    c = np.concatenate([p[j] * u for j in range(J)])
    A_eq = np.zeros((J, J * C))
    for j in range(J): A_eq[j, j * C:(j + 1) * C] = 1
    bnds = [(0, ts[j][cc]) for j in range(J) for cc in range(C)]
    v, _ = lp_max(c, A_eq=A_eq, b_eq=[x] * J, bounds=bnds)
    maxdiff = max(maxdiff, abs(v - sum(p[j] * h_unit_2xC(u, x, ts[j]) for j in range(J))))
report("N1b lifted LP = Minkowski sum", maxdiff < 1e-8, f"60 cases, max |diff| = {maxdiff:.1e}")

# N2: R x 2 orientation
maxdiff, ncase = 0.0, 0
for R in [2, 3, 4, 5]:
    for _ in range(40):
        J = rng.integers(2, 6)
        x = rng.dirichlet(np.ones(R))
        ts = rng.uniform(0, 1, size=J)
        p = rng.dirichlet(np.ones(J))
        rhs = lambda A: sum(p[j] * min(x[list(A)].sum(), ts[j]) for j in range(J))
        for _ in range(5):
            u = rng.normal(size=R)
            h_exact = sum(p[j] * h_unit_lp_Rx2(u, x, ts[j]) for j in range(J))
            h_sub = h_subset_polytope(u, rhs, float(p @ ts), R, list(proper_subsets(R)))
            maxdiff = max(maxdiff, abs(h_exact - h_sub)); ncase += 1
report("N2 subset identity Rx2", maxdiff < 1e-7,
       f"{ncase} comparisons (R=2..5), max |diff| = {maxdiff:.1e}")

# N3: coordinate-only (singletons + complements, i.e. averaged cell limits) is exact for
#     C <= 3, and can be strictly looser for C = 4 (random search + the review's example)
def coord_subsets(C):
    s = [(c,) for c in range(C)] + [tuple(sorted(set(range(C)) - {c})) for c in range(C)]
    return list(dict.fromkeys(s))

worst = {3: 0.0, 4: 0.0, 5: 0.0}
nstrict = {3: 0, 4: 0, 5: 0}
for C in [3, 4, 5]:
    for _ in range(150):
        J = rng.integers(2, 5); x = rng.uniform(.2, .8)
        ts = rng.dirichlet(np.ones(C) * 0.5, size=J); p = rng.dirichlet(np.ones(J))
        rhs = lambda A: sum(p[j] * min(x, ts[j][list(A)].sum()) for j in range(J))
        u = rng.normal(size=C)
        h_full = h_subset_polytope(u, rhs, x, C, list(proper_subsets(C)))
        h_coord = h_subset_polytope(u, rhs, x, C, coord_subsets(C))
        d = h_coord - h_full
        worst[C] = max(worst[C], d); nstrict[C] += d > 1e-7
report("N3a coordinate-only exact at C=3", worst[3] < 1e-7, f"150 cases, max excess = {worst[3]:.1e}")
report("N3b coordinate-only looser at C>=4", nstrict[4] > 0 and nstrict[5] > 0,
       f"strictly looser in {nstrict[4]}/150 (C=4), {nstrict[5]}/150 (C=5); "
       f"max excess {worst[4]:.3f}, {worst[5]:.3f}")

# the review's 2x4 example: x=1/2, T=(1/2,1/2,0,0) or (0,0,1/2,1/2) equally likely
x = 0.5; ts = np.array([[.5, .5, 0, 0], [0, 0, .5, .5]]); p = np.array([.5, .5])
rhs = lambda A: sum(p[j] * min(x, ts[j][list(A)].sum()) for j in range(2))
u = np.array([1, 1, 0, 0.])
hi_coord = h_subset_polytope(u, rhs, x, 4, coord_subsets(4))
lo_coord = -h_subset_polytope(-u, rhs, x, 4, coord_subsets(4))
hi_full = h_subset_polytope(u, rhs, x, 4, list(proper_subsets(4)))
lo_full = -h_subset_polytope(-u, rhs, x, 4, list(proper_subsets(4)))
cand = np.array([.25, .25, 0, 0])
cand_ok_coord = all(cand[list(A)].sum() <= rhs(A) + TOL for A in coord_subsets(4))
cand_ok_full = all(cand[list(A)].sum() <= rhs(A) + TOL for A in proper_subsets(4))
report("N3c review 2x4 example", cand_ok_coord and not cand_ok_full and abs(hi_full - .25) < TOL
       and abs(lo_full - .25) < TOL,
       f"group-1 mass to columns 1+2: coordinate-only [{lo_coord:.3f}, {hi_coord:.3f}], "
       f"exact [{lo_full:.3f}, {hi_full:.3f}] (as a share of group 1: [{lo_coord/x:.2f},{hi_coord/x:.2f}] vs {hi_full/x:.2f}); "
       f"candidate (1/4,1/4,0,0) admitted by coordinates: {cand_ok_coord}, by subsets: {cand_ok_full}")

# N4: size conditioning example (review Part I.3): at x=1/2, small precincts T in {0,1},
#     large precincts T=1/2. Group-1 conditional mean rate.
def dd_rate_2x2(x, t):  # group-1 rate interval in a 2x2 unit
    return max(0.0, (t - (1 - x)) / x), min(1.0, t / x)
small = [dd_rate_2x2(.5, 0.), dd_rate_2x2(.5, 1.)]
large = [dd_rate_2x2(.5, .5)]
s_lo, s_hi = np.mean([a for a, b in small]), np.mean([b for a, b in small])
l_lo, l_hi = large[0]
pooled = [(0.25, dd_rate_2x2(.5, 0.)), (0.5, dd_rate_2x2(.5, .5)), (0.25, dd_rate_2x2(.5, 1.))]
p_lo = sum(w * a for w, (a, b) in pooled); p_hi = sum(w * b for w, (a, b) in pooled)
report("N4 size conditioning", (s_lo, s_hi) == (.5, .5) and (p_lo, p_hi) == (.25, .75),
       f"small: [{s_lo},{s_hi}], large: [{l_lo},{l_hi}], intersection (size invariance) = 1/2; "
       f"pooled over size: [{p_lo},{p_hi}]")

# N5: plain-language M-SHARP example: x=1/2, half the precincts T=.1, half T=.9
lo = np.mean([dd_rate_2x2(.5, .1)[0], dd_rate_2x2(.5, .9)[0]])
hi = np.mean([dd_rate_2x2(.5, .1)[1], dd_rate_2x2(.5, .9)[1]])
report("N5 M-SHARP worked example", abs(lo - .4) < TOL and abs(hi - .6) < TOL,
       f"regression-only: group-1 mean at x=1/2 in [0,1]; with the precinct law: [{lo:.2f}, {hi:.2f}]")

# ---------------------------------------------------------------------------
# Count layer: worked examples (two precincts, veterans V / others O, outcome yes)
# P1: 100 men, 20 V, 30 yes.  P2: 100 men, 60 V, 50 yes.
# ---------------------------------------------------------------------------
N = np.array([100, 100]); V = np.array([20, 60]); Y = np.array([30, 50])

def dd_count(N, V, Y):
    return max(0, V + Y - N), min(V, Y)

d1, d2 = dd_count(N[0], V[0], Y[0]), dd_count(N[1], V[1], Y[1])
lo, hi = (d1[0] + d2[0]) / V.sum(), (d1[1] + d2[1]) / V.sum()
report("N6 DD worked example", (lo, hi) == (0.125, 0.875) and d1 == (0, 20) and d2 == (10, 50),
       f"P1 veteran-yes in {tuple(map(int,d1))}, P2 in {tuple(map(int,d2))}; veteran yes-rate in [{lo}, {hi}]")

# realised ordering in each precinct: y/V >= (Y-y)/(N-V)  <=>  y >= V*Y/N
ordlo = [int(np.ceil(V[i] * Y[i] / N[i])) for i in range(2)]
lo_o = (max(ordlo[0], d1[0]) + max(ordlo[1], d2[0])) / V.sum()
report("N7 realised ordering", ordlo == [6, 30] and lo_o == .45,
       f"lower limits become {ordlo}; veteran yes-rate in [{lo_o}, {hi}]")

# N8: audit additivity on margins-only (random), and its failure with a coupling row
maxdiff = 0.0
for _ in range(200):
    I = rng.integers(2, 7)
    Ni = rng.integers(20, 300, size=I); Vi = np.array([rng.integers(1, n) for n in Ni])
    Yi = np.array([rng.integers(0, n + 1) for n in Ni])
    # target: total veteran-yes count (rate up to the known denominator)
    # LP over y_i with DD bounds
    b = [dd_count(Ni[i], Vi[i], Yi[i]) for i in range(I)]
    width = lp_max(np.ones(I), bounds=b)[0] + lp_max(-np.ones(I), bounds=b)[0]
    di = [bb[1] - bb[0] for bb in b]
    for i in range(I):
        truth = rng.integers(b[i][0], b[i][1] + 1)
        b2 = list(b); b2[i] = (truth, truth)
        w2 = lp_max(np.ones(I), bounds=b2)[0] + lp_max(-np.ones(I), bounds=b2)[0]
        maxdiff = max(maxdiff, abs((width - w2) - di[i]))
report("N8a audit gain = d_i (separable)", maxdiff < 1e-7,
       f"200 random designs, every precinct and a random true table: max |gain - d_i| = {maxdiff:.1e}")

# coupling: a published total for precincts 1+2 (example): y1 + y2 = 45 (veteran yes total)
A_eq = np.array([[1, 1.]]); b = [d1, d2]
y2lo = -lp_max([0, -1], A_eq=A_eq, b_eq=[45], bounds=b)[0]; y2hi = lp_max([0, 1], A_eq=A_eq, b_eq=[45], bounds=b)[0]
report("N8b published total (E5) worked example", (y2lo, y2hi) == (25, 45),
       f"with the town total 45: P2 veteran-yes in [{y2lo:.0f},{y2hi:.0f}] "
       f"(rate [{y2lo/60:.3f},{y2hi/60:.3f}]) vs DD [10,50] (rate [{10/60:.3f},{50/60:.3f}])")

# with a coupling row, gains are no longer the separable d_i: target = P1 veteran-yes, a
# published total y1 + y2 = 45. Without the total, auditing P2 removes nothing from a P1 target
# (P2 contributes d = 0 to it); with the total, auditing P2 pins P1 completely.
def p1_width(bounds, coup):
    kw = dict(A_eq=np.array([[1, 1.]]), b_eq=[45]) if coup else {}
    c = np.array([1, 0.])
    return lp_max(c, bounds=bounds, **kw)[0] + lp_max(-c, bounds=bounds, **kw)[0]
bb = [d1, d2]
out = []
for coup in (False, True):
    w0 = p1_width(bb, coup)
    g = [w0 - p1_width([d1, (t, t)], coup) for t in (25, 35, 45)]
    out.append((w0, g))
report("N8c audit gain not separable with a coupling row",
       out[0][1] == [0, 0, 0] and out[1][0] == 20 and out[1][1] == [20, 20, 20],
       f"target P1: no total: width {out[0][0]:.0f}, auditing P2 removes {out[0][1]}; "
       f"with the total 45: width {out[1][0]:.0f}, auditing P2 removes {out[1][1]} (P2 truth 25/35/45)")

# N8d: with coupling rows, the remaining width after auditing unit j is concave in unit j's
# revealed count (max-LP value is concave, min-LP value convex in the right-hand side), so the
# worst-case remaining width is a one-dimensional concave maximisation per candidate unit.
viol = 0; nchk = 0
for _ in range(200):
    I = 4
    Ns = rng.integers(20, 120, size=I); Vs = np.array([rng.integers(1, n) for n in Ns])
    Ys = np.array([rng.integers(0, n + 1) for n in Ns]); b = [dd_count(Ns[i], Vs[i], Ys[i]) for i in range(I)]
    G = rng.normal(size=(2, I)); y0 = np.array([rng.integers(bb[0], bb[1] + 1) for bb in b])
    h = G @ y0 + rng.uniform(0, 3, size=2)        # two coupling rows, feasible at y0
    c = rng.uniform(0, 1, size=I); c[1] = 0        # target excludes unit 2 (index 1)
    ws = []
    for tval in range(b[1][0], b[1][1] + 1):
        b2 = list(b); b2[1] = (tval, tval)
        r1 = linprog(-c, A_ub=G, b_ub=h, bounds=b2, method="highs"); r2 = linprog(c, A_ub=G, b_ub=h, bounds=b2, method="highs")
        ws.append(-r1.fun - r2.fun if r1.status == 0 and r2.status == 0 else None)
    w = [v for v in ws if v is not None]
    # feasible t form an interval (projection of a polytope); check discrete concavity on it
    for k in range(1, len(w) - 1):
        nchk += 1; viol += w[k] < (w[k - 1] + w[k + 1]) / 2 - 1e-7
report("N8d remaining width concave in the audited count", viol == 0,
       f"200 random 4-unit designs with two coupling rows: {viol} concavity violations in {nchk} checks")

# N9: stratified margins (E5a) worked example: P2 split by age band
# young: 30 men, 25 V, 10 yes ; old: 70 men, 35 V, 40 yes
yb = dd_count(30, 25, 10); ob = dd_count(70, 35, 40)
report("N9 stratified margins", (yb[0] + ob[0], yb[1] + ob[1]) == (10, 45),
       f"young {tuple(map(int,yb))}, old {tuple(map(int,ob))}: P2 veteran-yes in [{yb[0]+ob[0]},{yb[1]+ob[1]}] vs [10,50] unstratified")

# N10: integrality. margins-only endpoints are integer; a realised cross-cell row can make
#      the LP endpoint fractional (so an integer check is needed).
nonint_margins = 0
for _ in range(300):
    R, C = rng.integers(2, 4), rng.integers(2, 5)
    n = rng.integers(0, 30, size=(R, C)); rs, cs = n.sum(1), n.sum(0)
    # LP over table cells with margins; random integer target
    A_eq = []; b_eq = []
    for r in range(R):
        row = np.zeros(R * C); row[r * C:(r + 1) * C] = 1; A_eq.append(row); b_eq.append(rs[r])
    for c in range(C):
        row = np.zeros(R * C); row[c::C] = 1; A_eq.append(row); b_eq.append(cs[c])
    w = rng.integers(-3, 4, size=R * C)
    v, sol = lp_max(w, A_eq=np.array(A_eq), b_eq=b_eq, bounds=[(0, None)] * (R * C))
    nonint_margins += abs(v - round(v)) > 1e-7
report("N10a integer endpoints, margins only", nonint_margins == 0,
       f"300 random integer tables (2-3 x 2-4), random integer targets: {nonint_margins} fractional LP optima")

# search for a 2x2 unit with realised ordering where LP and integer endpoints differ
found = None
for Nn in range(5, 60):
    for Vv in range(1, Nn):
        for Yy in range(0, Nn + 1):
            lo_lp = max(Vv * Yy / Nn, max(0, Vv + Yy - Nn))
            lo_int = max(int(np.ceil(Vv * Yy / Nn - 1e-12)), max(0, Vv + Yy - Nn))
            hi = min(Vv, Yy)
            if lo_lp <= hi and abs(lo_lp - round(lo_lp)) > 1e-9 and lo_int <= hi:
                found = (Nn, Vv, Yy, lo_lp, lo_int); break
        if found: break
    if found: break
report("N10b ordering row gives fractional LP endpoint", found is not None,
       f"N={found[0]}, V={found[1]}, yes={found[2]}: LP lower veteran-yes {found[3]:.3f}, integer {found[4]} "
       "(single unit: rounding fixes it; with rows coupling units it need not)")

# a coupled case: three units, a published veteran-yes total for all three (equality) plus
# two pooled realised ordering rows with different unit weights (as the realised delta_B and
# delta_W rows of 2c M1 weight units differently):  sum_i w_i (y_i/V_i - (Y_i-y_i)/(N_i-V_i)) >= 0.
# Search for an LP-feasible, integer-infeasible system, or an integer endpoint above the rounded
# LP endpoint. (One coupling row plus one total cannot do it: box x total is TU and the row's
# max over it is attained at an integer vertex.)
def coupled(Ns, Vs, Ys, W, total):
    I = len(Ns)
    rows = np.array([[w[i] * (1 / Vs[i] + 1 / (Ns[i] - Vs[i])) for i in range(I)] for w in W])
    rhs = np.array([sum(w[i] * Ys[i] / (Ns[i] - Vs[i]) for i in range(I)) for w in W])
    b = [dd_count(Ns[i], Vs[i], Ys[i]) for i in range(I)]
    c = np.zeros(I); c[0] = 1           # target: unit 1's veteran-yes count
    r = linprog(c, A_ub=-rows, b_ub=-rhs, A_eq=np.ones((1, I)), b_eq=[total], bounds=b, method="highs")
    if r.status != 0: return None
    cons = [LinearConstraint(rows, rhs, np.inf), LinearConstraint(np.ones((1, I)), total, total)]
    m = milp(c, constraints=cons, integrality=np.ones(I),
             bounds=Bounds([bb[0] for bb in b], [bb[1] for bb in b]))
    return (r.fun, None) if m.status != 0 else (r.fun, m.fun)
gap_found = None
for _ in range(20000):
    I = 3
    Ns = rng.integers(6, 30, size=I); Vs = np.array([rng.integers(1, n) for n in Ns])
    Ys = np.array([rng.integers(0, n + 1) for n in Ns]); W = rng.dirichlet(np.ones(I), size=2)
    b = [dd_count(Ns[i], Vs[i], Ys[i]) for i in range(I)]
    total = rng.integers(sum(x[0] for x in b), sum(x[1] for x in b) + 1)
    g = coupled(Ns, Vs, Ys, W, total)
    if g is None: continue
    if g[1] is None:
        gap_found = ("LP feasible but no integer table", Ns, Vs, Ys, total, g); break
    if np.ceil(g[0] - 1e-9) < g[1] - 1e-9:
        gap_found = ("integer endpoint above the rounded LP endpoint", Ns, Vs, Ys, total, g); break
report("N10c coupled rows: LP and integer answers can differ beyond rounding", gap_found is not None,
       (f"{gap_found[0]}: N={[int(v) for v in gap_found[1]]}, V={[int(v) for v in gap_found[2]]}, yes={[int(v) for v in gap_found[3]]}, "
        f"published veteran-yes total {int(gap_found[4])}, LP lower end for unit 1 = {gap_found[5][0]:.3f}, "
        f"integer = {gap_found[5][1]}") if gap_found else "no case found")

# N11: ballot patterns. One precinct, 100 voters, 50 in each group A,B; two binary contests.
def pattern_bounds(counts, nA):
    # variables: group-A counts in patterns 00,01,10,11; target: (n_A10 - n_A01)/nA
    pats = ["00", "01", "10", "11"]
    c = np.zeros(4); c[pats.index("10")] = 1 / nA; c[pats.index("01")] = -1 / nA
    bnds = [(0, counts[p]) for p in pats]
    A_eq = np.ones((1, 4)); b_eq = [nA]
    hi = lp_max(c, A_eq=A_eq, b_eq=b_eq, bounds=bnds)[0]
    lo = -lp_max(-c, A_eq=A_eq, b_eq=b_eq, bounds=bnds)[0]
    # one-contest support of A in contest 1 = (n_A10 + n_A11)/nA
    c1 = np.zeros(4); c1[2] = c1[3] = 1 / nA
    s_hi = lp_max(c1, A_eq=A_eq, b_eq=b_eq, bounds=bnds)[0]; s_lo = -lp_max(-c1, A_eq=A_eq, b_eq=b_eq, bounds=bnds)[0]
    return (lo, hi), (s_lo, s_hi)
d, s = pattern_bounds({"00": 50, "01": 0, "10": 0, "11": 50}, 50)
d4, s4 = pattern_bounds({"00": 48, "01": 0, "10": 4, "11": 48}, 50)
report("N11 ballot patterns", d == (0, 0) and s == (0, 1) and abs(d4[1] - .08) < TOL,
       f"no split tickets: A's support in contest 1 in [{s[0]:.0f},{s[1]:.0f}] but its contest-1 minus contest-2 "
       f"difference = [{d[0]:.0f},{d[1]:.0f}]; with 4 split (10) ballots: difference in [{d4[0]:.2f},{d4[1]:.2f}], "
       f"support in [{s4[0]:.2f},{s4[1]:.2f}]")

# N12: finite-composition contrast (Kitagawa) identity and the review's numbers
maxdiff = 0.0
for _ in range(500):
    R = rng.integers(2, 5)
    x, z = rng.dirichlet(np.ones(R)), rng.dirichlet(np.ones(R))
    mx, mz = rng.uniform(size=R), rng.uniform(size=R)
    lhs = z @ mz - x @ mx
    rhs = (z - x) @ ((mx + mz) / 2) + ((x + z) / 2) @ (mz - mx)
    maxdiff = max(maxdiff, abs(lhs - rhs))
report("N12a symmetric decomposition identity", maxdiff < 1e-12, f"500 random cases, max |diff| = {maxdiff:.1e}")

# under LCE (m(x) = H x, one column) with z - x = t (e_r - e_s): mix term = t * O_rs(midpoint) and
# rate term = t * L_rs(midpoint) (PACKAGE_DESIGN P-S objects): the secant signs are the local
# signs of M-SIGN/M-LOC evaluated at the midpoint only.
maxdiff = 0.0
for _ in range(500):
    R = rng.integers(2, 5); H = rng.normal(size=(R, R))
    r_, s_ = rng.choice(R, 2, replace=False); x = rng.dirichlet(np.ones(R) * 3)
    t = min(x[s_], 1) * rng.uniform(.1, 1)
    z = x.copy(); z[r_] += t; z[s_] -= t
    mx, mz, xb = H @ x, H @ z, (x + z) / 2
    mix = (z - x) @ ((mx + mz) / 2); rate = xb @ (mz - mx)
    O = (H @ xb)[r_] - (H @ xb)[s_]; L = (xb @ H)[r_] - (xb @ H)[s_]
    maxdiff = max(maxdiff, abs(mix - t * O), abs(rate - t * L))
report("N12c LCE: secant terms = local P-S objects at the midpoint", maxdiff < 1e-12,
       f"500 random H, pairs and steps: max |diff| = {maxdiff:.1e}")

# bins L: x=.25 (group-1 share), yes .45 ; H: x=.75, yes .55. Unknowns: (a_L, b_L, a_H, b_H) =
# group-1 and group-0 rates in each bin. Constraints: accounting; rates in [0,1].
# Target: average gap g = ((a_L-b_L)+(a_H-b_H))/2.
def contrast_bounds(signs):
    xL, xH, tL, tH = .25, .75, .45, .55
    A_eq = np.array([[xL, 1 - xL, 0, 0], [0, 0, xH, 1 - xH]]); b_eq = [tL, tH]
    A_ub, b_ub = [], []
    t = xH - xL
    if signs:
        # mix term t*gbar >= 0  and rate term xbar.(m_H - m_L) >= 0
        A_ub.append(-t * np.array([.5, -.5, .5, -.5])); b_ub.append(0)
        xb = np.array([(xL + xH) / 2, 1 - (xL + xH) / 2])
        A_ub.append(-np.array([-xb[0], -xb[1], xb[0], xb[1]])); b_ub.append(0)
    c = np.array([.5, -.5, .5, -.5])
    kw = dict(A_eq=A_eq, b_eq=b_eq, bounds=[(0, 1)] * 4)
    if A_ub: kw.update(A_ub=np.array(A_ub), b_ub=b_ub)
    hi = lp_max(c, **kw)[0]; lo = -lp_max(-c, **kw)[0]
    return lo, hi
nos = contrast_bounds(False); withs = contrast_bounds(True)
report("N12b contrast worked example", abs(withs[0]) < TOL and abs(withs[1] - .2) < TOL,
       f"average group gap: accounting only [{nos[0]:.3f},{nos[1]:.3f}]; with both signs [{withs[0]:.3f},{withs[1]:.3f}]")

# N13: exception-mass worked example on the two-precinct data (ordering may fail in
#      precincts holding at most a share pi of voters)
def exc_lower(pi):
    best = None
    for k in range(3):
        for S in itertools.combinations(range(2), k):
            if sum(N[list(S)]) / N.sum() > pi + 1e-12: continue
            lo = 0
            for i in range(2):
                ddl = dd_count(N[i], V[i], Y[i])[0]
                lo += ddl if i in S else max(ddl, ordlo[i])
            best = lo if best is None else min(best, lo)
    return best / V.sum()
fr = {pi: float(exc_lower(pi)) for pi in [0, .5, 1]}
report("N13 exception-mass example", fr == {0: .45, .5: .2, 1: .125},
       f"lower end of veteran yes-rate by exception share pi: {fr} "
       "(pi = .5: either precinct may be the exception; the worst is P2)")

# N14: misclassification example (P2: 100 men, recorded 60 V, 50 yes; true V in [58, 62])
cands = []
for Vt in range(58, 63):
    lo_, hi_ = dd_count(100, Vt, 50)
    cands.append((lo_ / Vt, hi_ / Vt))
mlo, mhi = min(a for a, b in cands), max(b for a, b in cands)
report("N14 misclassification example", abs(mlo - 8 / 58) < TOL and abs(mhi - 50 / 58) < TOL,
       f"recorded labels: [{10/60:.3f},{50/60:.3f}]; up to 2 mislabelled: [{mlo:.3f},{mhi:.3f}]")

# N15: E4 transfer worked example (2x2, pooled LCE).  tau(x) = .3 + .3x + .1x^2 is identified;
# m0 = .3 + k x, m1 = (.6 - k) + (.1 + k) x, true k = .1.  Electorate A has precincts at x = .1, .3;
# electorate B at x = .6, .9 (equal sizes).  A's veteran rate is published (model-implied, +/- d).
def e4_example(d):
    xs_all = [.1, .3, .6, .9]; xa = np.array([.1, .3]); xb = np.array([.6, .9])
    m1 = lambda x, k: (.6 - k) + (.1 + k) * x
    m0 = lambda x, k: .3 + k * x
    ks = np.linspace(-2, 2, 400001)
    ok = np.ones_like(ks, bool)
    for x in xs_all:
        ok &= (m0(x, ks) >= 0) & (m0(x, ks) <= 1) & (m1(x, ks) >= 0) & (m1(x, ks) <= 1)
    rate = lambda xv, k: np.sum(xv * m1(xv, k)) / np.sum(xv)      # group-1 weighted
    if d is not None:
        trueA = rate(xa, .1)
        rA = np.array([rate(xa, k) for k in ks[ok]]) if False else (np.sum(xa * (.6 + .1 * xa)) + ks * np.sum(xa * (xa - 1))) / np.sum(xa)
        ok &= np.abs(rA - trueA) <= d + 1e-12
    rB = (np.sum(xb * (.6 + .1 * xb)) + ks[ok] * np.sum(xb * (xb - 1))) / np.sum(xb)
    return rB.min(), rB.max(), rate(xb, .1)
lo0, hi0, tB = e4_example(None); lo1, hi1, _ = e4_example(0.0); lo2, hi2, _ = e4_example(0.02)
report("N15 E4 transfer worked example", abs(lo1 - tB) < 1e-4 and abs(hi1 - tB) < 1e-4 and lo2 < tB < hi2,
       f"electorate B veteran rate (model-implied, truth {tB:.3f}): LCE alone [{lo0:.3f},{hi0:.3f}]; "
       f"with A's table, same pattern exactly: [{lo1:.3f},{hi1:.3f}]; pattern may differ (A's rate +/- .02): [{lo2:.3f},{hi2:.3f}]")

# N16: Avenue D (polynomial ladder) worked example. 2x2; truth m0 = .3 + .1x, m1 = .5 + .2x
# (so the straight line is true); precinct types at x = .1, .3, .5, .7, .9.  Contextual effects of
# degree <= J: m_r(x) = sum_j a_rj x^j.  The accounting identity tau(x) = x m1 + (1-x) m0 is imposed as
# a polynomial identity (tau identified, degree J+1), leaving J free directions. Outer version: rates
# valid at the support types and on a fine grid of the hull. Target: group-1 rate, group-1 weighted.
def poly_bounds(J, grid):
    xs = np.array([.1, .3, .5, .7, .9])
    tau = lambda x: x * (.5 + .2 * x) + (1 - x) * (.3 + .1 * x)   # .3 + .3x + .1x^2
    tc = np.array([.3, .3, .1] + [0.0] * (J - 1))                    # coefficients up to degree J+1
    nv = 2 * (J + 1)                                                  # a1_0..a1_J, a0_0..a0_J
    A_eq = np.zeros((J + 2, nv)); b_eq = tc[:J + 2]
    for j in range(J + 1):
        A_eq[j + 1, j] += 1                  # x * a1_j x^j -> degree j+1
        A_eq[j, J + 1 + j] += 1              # a0_j x^j -> degree j
        A_eq[j + 1, J + 1 + j] -= 1          # -x a0_j x^j -> degree j+1
    pts = np.unique(np.concatenate([xs, grid]))
    A_ub, b_ub = [], []
    for x in pts:
        p = x ** np.arange(J + 1)
        for off in (0, J + 1):
            row = np.zeros(nv); row[off:off + J + 1] = p
            A_ub.append(row); b_ub.append(1); A_ub.append(-row); b_ub.append(0)
    w = xs / xs.sum()
    c = np.zeros(nv); c[:J + 1] = (w[:, None] * xs[:, None] ** np.arange(J + 1)).sum(0)
    kw = dict(A_ub=np.array(A_ub), b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * nv)
    return -lp_max(-c, **kw)[0], lp_max(c, **kw)[0]
grid = np.linspace(.1, .9, 81)
l1, u1 = poly_bounds(1, grid); l2, u2 = poly_bounds(2, grid); l3, u3 = poly_bounds(3, grid)
truth = float(np.sum((np.array([.1,.3,.5,.7,.9]) / 2.5) * (.5 + .2 * np.array([.1,.3,.5,.7,.9]))))
dd = [dd_rate_2x2(x, .3 + .3 * x + .1 * x * x) for x in [.1, .3, .5, .7, .9]]
wv = np.array([.1, .3, .5, .7, .9]) / 2.5
ddl, ddu = float(wv @ [a for a, b in dd]), float(wv @ [b for a, b in dd])
report("N16 polynomial ladder widens", l1 >= l2 - 1e-9 and u1 <= u2 + 1e-9 and l2 >= l3 - 1e-9 and u2 <= u3 + 1e-9,
       f"group-1 rate (truth {truth:.3f}): straight lines [{l1:.3f},{u1:.3f}]; quadratic [{l2:.3f},{u2:.3f}]; "
       f"cubic [{l3:.3f},{u3:.3f}]; Duncan-Davis (no model) [{ddl:.3f},{ddu:.3f}]")

# N17: odds-ratio dial (NEW_AVENUES 12.3(c)) on precinct P2 (100 men, 60 veterans, 50 yes).
# Cell y = veteran yes; OR(y) = y (y - 10) / ((60 - y)(50 - y)), increasing in y on [10, 50].
def or_range(lo_or, hi_or):
    ys = np.linspace(10, 50, 400001)[1:-1]
    orr = ys * (ys - 10) / ((60 - ys) * (50 - ys))
    ok = (orr >= lo_or) & (orr <= hi_or)
    return ys[ok].min() / 60, ys[ok].max() / 60
o1 = or_range(1, np.inf); o2 = or_range(1, 4)
report("N17 odds-ratio dial example", abs(o1[0] - .5) < 1e-3 and o2[1] < 50 / 60,
       f"P2 veteran yes rate: DD [{10/60:.3f},{50/60:.3f}]; OR >= 1 (ordering) [{o1[0]:.3f},{o1[1]:.3f}]; "
       f"1 <= OR <= 4 [{o2[0]:.3f},{o2[1]:.3f}]")

# N18: M-CONTRAST "hull" implies the model between signs. delta_B(r,s) = E_p[X_r X_s (m_r - m_s)]
# has non-negative weights, so O_rs(x) = m_r(x) - m_s(x) >= 0 at every support point implies
# delta_B(r,s) >= 0 (redundant when both are imposed; checked on random LCE designs).
viol = 0
for _ in range(2000):
    R = rng.integers(2, 5); I = rng.integers(3, 12)
    X = rng.dirichlet(np.ones(R), size=I); p = rng.dirichlet(np.ones(I)); H = rng.normal(size=(R, R))
    M = X @ H.T                                   # m(x_i) = H x_i, rows = units
    r_, s_ = 0, 1
    O = M[:, r_] - M[:, s_]
    if np.all(O >= 0):
        dB = np.sum(p * X[:, r_] * X[:, s_] * O)
        viol += dB < -1e-12
report("N18 hull contrasts imply the model between sign", viol == 0, f"{viol} violations among random designs with O_rs >= 0 at all support points")

# N19: estimated (probabilistically matched) composition. If each of a unit's N men is a veteran
# with known match probability q_j, independently, the true share X = qbar + u with E[u | q] = 0
# (Berkson error) and Var(u | q) = sum_j q_j (1 - q_j) / N^2.  Then for one column (2 groups,
# X = veteran share) E[T | q] = qbar' Q qbar + Var(u) * (Q_vv - 2 Q_vn + Q_nn): a KNOWN correction
# to the regression that involves only the symmetric part Q, so K stays the only unknown.
# Assumes the match probabilities carry no information about rates beyond the true composition.
Hm = np.array([[.55, .40], [.20, .30]])         # m_v(x) = .55 x + .40 (1-x) etc. (rows v, n; cols v-share, n-share)
Qm = (Hm + Hm.T) / 2
maxdiff = 0.0
for _ in range(20):
    Nn = rng.integers(20, 200); q = rng.beta(.5, .5, size=Nn) * rng.uniform(.2, 1)
    qbar = q.mean(); var_u = np.sum(q * (1 - q)) / Nn ** 2
    sims = 200000
    V = (rng.uniform(size=(sims, Nn)) < q).sum(1) / Nn       # true veteran share draws
    Xs = np.stack([V, 1 - V], 1)
    ET_mc = np.mean(np.einsum("ij,jk,ik->i", Xs, Qm, Xs))   # E[X'QX | q]; rates' noise averages out
    xb = np.array([qbar, 1 - qbar])
    ET_formula = xb @ Qm @ xb + var_u * (Qm[0, 0] - 2 * Qm[0, 1] + Qm[1, 1])
    maxdiff = max(maxdiff, abs(ET_mc - ET_formula))
report("N19 Berkson correction for matched composition", maxdiff < 5e-4,
       f"20 random units (20-200 men, match probabilities in [0,1]): max |MC - formula| = {maxdiff:.1e} (MC error ~1e-4)")

# N20: count layer with an estimated group margin.  P2: 100 men, 50 yes.  Matching gives 50 men
# with probability 1 and 20 with probability .5 (expected veterans 60, sd 2.24).  An exact
# (Poisson-binomial) central 95% set for the true count is used as the margin set.
from math import comb
pmf = np.array([comb(20, k) * .5 ** 20 for k in range(21)])            # extra veterans among the 20
cdf = np.cumsum(pmf); lo_k = int(np.searchsorted(cdf, .025)); hi_k = int(np.searchsorted(cdf, .975))
Vset = range(50 + lo_k, 50 + hi_k + 1)
rates = [(dd_count(100, Vv, 50)[0] / Vv, dd_count(100, Vv, 50)[1] / Vv) for Vv in Vset]
rl, ru = min(a for a, b in rates), max(b for a, b in rates)
report("N20 estimated margin widens the count-layer range", rl < 10 / 60 and ru > 50 / 60,
       f"expected count 60 treated as exact: [{10/60:.3f},{50/60:.3f}]; true count in its 95% set "
       f"[{min(Vset)},{max(Vset)}]: [{rl:.3f},{ru:.3f}] (covers the actual rate with prob >= .95 under calibrated, independent matches)")

# N21: relative-magnitude (Rambachan-Roth-style) causal bounds when pre-period trends are themselves
# EI bounds. Two precincts (100 men each; 20 and 60 veterans), three elections t = -1, 0, 1 (post).
# Latent veteran-yes totals s_t = y_1t + y_2t lie in the sum of the precincts' DD boxes. The gap
# g_t = B_v,t - B_n,t = s_t (1/80 + 1/120) - Y_t / 120 is linear in s_t. Pre-period DiD D0 = g0 - g_-1,
# post DiD D1 = g1 - g0 (they SHARE s_0). ATT = D1 - b1 with |b1| <= Mbar |b0|, b0 = D0 (no treatment
# before t = 1). Joint: optimise over (s_-1, s_0, s_1) together. Endpoint-wise: take D0's and D1's
# intervals separately and combine extremes (can pair incompatible histories).
Yv = {-1: (30, 45), 0: (32, 48), 1: (30, 60)}
def sbox(t):
    (a, b) = Yv[t]
    l1, u1 = max(0, 20 + a - 100), min(20, a); l2, u2 = max(0, 60 + b - 100), min(60, b)
    return l1 + l2, u1 + u2
cf = 1 / 80 + 1 / 120
g = lambda t, st: st * cf - sum(Yv[t]) / 120
def att_bounds(Mbar, boxes, joint=True, grid=161):
    S = {t: np.linspace(*boxes[t], grid) for t in (-1, 0, 1)}
    if joint:
        sm1, s0, s1 = np.meshgrid(S[-1], S[0], S[1], indexing="ij")
        D0 = g(0, s0) - g(-1, sm1); D1 = g(1, s1) - g(0, s0)
        return float((D1 - Mbar * np.abs(D0)).min()), float((D1 + Mbar * np.abs(D0)).max()), \
               (float(D0.min()), float(D0.max())), (float(D1.min()), float(D1.max()))
    sm1, s0 = np.meshgrid(S[-1], S[0], indexing="ij"); D0 = g(0, s0) - g(-1, sm1)
    s0b, s1 = np.meshgrid(S[0], S[1], indexing="ij"); D1 = g(1, s1) - g(0, s0b)
    a0 = float(np.abs(D0).max())
    return float(D1.min() - Mbar * a0), float(D1.max() + Mbar * a0), None, None
B = {t: sbox(t) for t in (-1, 0, 1)}
jl, ju, D0r, D1r = att_bounds(1.0, B, True); el, eu, _, _ = att_bounds(1.0, B, False)
# trend dial with M = .10 for comparison, and a "narrow pre-trend" variant: a model rung pins the
# pre-period totals s_-1, s_0 to +/- 2 votes around a central value (illustrative).
td = (D1r[0] - .10, D1r[1] + .10)
Bn = dict(B); Bn[-1] = (B[-1][0] + (B[-1][1] - B[-1][0]) / 2 - 2, B[-1][0] + (B[-1][1] - B[-1][0]) / 2 + 2)
Bn[0] = (B[0][0] + (B[0][1] - B[0][0]) / 2 - 2, B[0][0] + (B[0][1] - B[0][0]) / 2 + 2)
nl, nu, nD0, nD1 = att_bounds(1.0, Bn, True)
report("N21 set-valued pre-trends in relative-magnitude bounds", jl >= el - 1e-9 and ju <= eu + 1e-9,
       f"arithmetic only: pre DiD D0 in [{D0r[0]:.2f},{D0r[1]:.2f}], post DiD D1 in [{D1r[0]:.2f},{D1r[1]:.2f}]; "
       f"ATT with Mbar=1: joint [{jl:.2f},{ju:.2f}], endpoint-wise [{el:.2f},{eu:.2f}]; trend dial M=.10: [{td[0]:.2f},{td[1]:.2f}]. "
       f"Narrow pre-trends (pre totals pinned to +/-2 votes): D0 in [{nD0[0]:.2f},{nD0[1]:.2f}], ATT joint [{nl:.2f},{nu:.2f}]")

# N22: covariate-benchmarked kappa (NEW_AVENUES 14.8) with a KNOWN true contextual violation.
# Units: covariate Z ~ U(-1,1); composition X = c0 + g Z + U(-w, w) (g, w set how much Z explains X).
# Rates: m1 = .5 + b Z + D1, m0 = .3 + .5 b Z + D0, plus unit noise (sd .08 per group), with the
# invisible contextual violation of 14.3: D1 = k (1 - X), D0 = -k X (x'D = 0: E[Y|X,Z] is untouched).
# Benchmark (drop Z, the only covariate): Delta^(-Z)_1 = b (E[Z|X] - E Z); kappa_Z = E[(.)^2]/Var(Y).
# Truth: CAR(Z) violation for group 1 = -k (X - E[X|Z]); kappa_true = E[(.)^2]/Var(Y).
# m* = kappa_true / kappa_Z: the multiple of the benchmark the truth needs.
# Overlap factor (14.5, 2x2): mean over Z of (2 - a - b)/(b - a) for the within-cell range [a, b] of X.
def kappa_case(c0, g, w, b, k, n=400000):
    Z = rng.uniform(-1, 1, n); X = c0 + g * Z + rng.uniform(-w, w, n)
    e1, e0 = rng.normal(0, .08, n), rng.normal(0, .08, n)
    m1 = .5 + b * Z + k * (1 - X) + e1; m0 = .3 + .5 * b * Z - k * X + e0
    Y = X * m1 + (1 - X) * m0; VY = Y.var()
    EZgX = np.polyval(np.polyfit(X, Z, 3), X)
    kZ = np.mean((b * (EZgX - Z.mean())) ** 2) / VY
    EXgZ = c0 + g * Z
    kt = np.mean((k * (X - EXgZ)) ** 2) / VY
    A = np.column_stack([np.ones(n), X, Z, X * Z, X ** 2])
    r2Y = 1 - np.var(Y - A @ np.linalg.lstsq(A, Y, rcond=None)[0]) / VY
    r2X = 1 - np.var(X - EXgZ) / X.var()
    a_, b_ = EXgZ - w, EXgZ + w
    ovl = np.mean((2 - a_ - b_) / (b_ - a_))
    # bias of the CAR(Z) answer for group 1's rate: CAR recovers eta_1 = .5 + b Z (the violation is
    # invisible), the truth is B_1 = E[X m1]/E[X]
    B1_true = np.sum(X * m1) / np.sum(X); B1_car = np.sum(X * (.5 + b * Z)) / np.sum(X)
    return dict(R2_Y=r2Y, R2_X=r2X, kZ=kZ, kt=kt, m=kt / kZ, ovl=ovl, bias=B1_car - B1_true)
cases = {"A covariate shares composition, small contextual effect": (.4, .15, .25, .4, .05),
         "B covariate explains outcome, not composition, large contextual effect": (.4, .0, .25, .4, .30),
         "C covariate absorbs most composition variation, large contextual effect": (.5, .35, .10, .4, .30),
         "D modest covariate tied to composition, small contextual effect": (.4, .15, .25, .1, .03)}
res = {nm: kappa_case(*v) for nm, v in cases.items()}
rA, rB, rC = res[list(cases)[0]], res[list(cases)[1]], res[list(cases)[2]]
ok = rA["m"] < 1 and rB["m"] > 10 and rC["kt"] < rB["kt"] / 5 and rC["ovl"] > rA["ovl"] and abs(rC["bias"]) > 4 * abs(rA["bias"])
report("N22 covariate-benchmarked kappa: when the benchmark is informative", ok,
       "; ".join(f"{nm}: R2(Y~X,Z) {r['R2_Y']:.2f}, R2(X~Z) {r['R2_X']:.2f}, kappa_Z {r['kZ']:.4f}, "
                 f"kappa_true {r['kt']:.4f}, m* {r['m']:.2f}, overlap factor {r['ovl']:.1f}, bias of CAR(Z) for B_1 {r['bias']:+.3f}" for nm, r in res.items()))

print(f"\n{sum(results)}/{len(results)} PASS")
