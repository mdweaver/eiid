"""
Checks for the reviewer's proposed avenues (NEW_AVENUES.md section 13).
 (i)   Multiple validation tables: rank of the stacked map K -> B*(K) over J known tables
       with different compositions (shared K, table-specific Q).
 (ii)  Crossing partitions: the same individuals aggregated two ways (e.g. precincts and wards whose
       boundaries cross).  Bounds on an aggregate rate from each partition's DD, their intersection,
       and the joint atom-level LP.
 (iii) Covariate-stratified margins inside units: group x age and outcome x age known per unit.
       DD for the group x outcome cell with and without the age stratification.
 (iv)  Optional (needs the owner's WI file as argv[1]): covariance of unit-level group rates identified
       from Var(T | x) = x' Sigma x under homoskedastic deviations (the CAR2-type idea), 2x2 yes-rate.
"""
import sys
import numpy as np
from scipy.optimize import linprog
rng = np.random.default_rng(12)
out = ["# Checks for the review items", ""]

# ---------- (i)
def kmap(R, C, mu):
    pairs = [(r, s) for r in range(R) for s in range(r + 1, R)]
    M = np.zeros((R * C, (C - 1) * len(pairs))); mur = mu.sum(1)
    for ci in range(C - 1):
        for pi, (r, s) in enumerate(pairs):
            col = ci * len(pairs) + pi
            for cc, sg in ((ci, 1.0), (C - 1, -1.0)):
                M[r * C + cc, col] += sg * mu[r, s] / mur[r]
                M[s * C + cc, col] -= sg * mu[r, s] / mur[s]
    return M
out.append("## (i) multiple validation tables, shared K")
out.append("| R | C | dim K | rank with J = 1, 2, 3, 4 tables | tables needed | overidentifying restrictions at that J |")
out.append("|---|---|---|---|---|---|")
for R, C in [(3, 2), (3, 3), (4, 2), (4, 3), (5, 3)]:
    ranks = []; Ms = []
    for J in range(1, 5):
        a = rng.dirichlet(np.ones(R) * 2) * 3                       # each table: its own composition distribution
        X = rng.dirichlet(a, 400); N = rng.integers(100, 1000, 400)
        Ms.append(kmap(R, C, np.einsum("i,ir,is->rs", N, X, X) / 400))
        ranks.append(np.linalg.matrix_rank(np.vstack(Ms), tol=1e-9))
    d = (C - 1) * R * (R - 1) // 2
    J = next((j + 1 for j, rk in enumerate(ranks) if rk == d), None)
    over = (J * (R - 1) * (C - 1) - d) if J else None
    out.append(f"| {R} | {C} | {d} | {', '.join(map(str, ranks))} | {J} | {over} |")
out.append("")

# ---------- (ii)
def fh(nr, nc, tot):        # bounds on cell (row A, col yes) given margins in one unit
    return max(0.0, nr + nc - tot), min(nr, nc)
out.append("## (ii) crossing partitions of the same individuals (2 groups x yes/no; 60 atoms; 12 units each way)")
out.append("| design | DD, partition 1 | DD, partition 2 | intersection | joint atom LP | truth |")
out.append("|---|---|---|---|---|---|")
for lab, het in (("homogeneous atoms within units", 0.05), ("heterogeneous atoms", 0.25)):
    na = 60
    gA = rng.integers(20, 200, na); gB = rng.integers(20, 200, na)
    base = rng.uniform(0.2, 0.8, na)
    pA = np.clip(base + 0.2 + het * rng.standard_normal(na), 0.01, 0.99); pB = np.clip(base - 0.2 + het * rng.standard_normal(na), 0.01, 0.99)
    yA = rng.binomial(gA, pA); yB = rng.binomial(gB, pB)
    P1 = np.repeat(np.arange(12), 5); P2 = rng.permutation(P1)          # two crossing partitions of atoms
    totA = gA.sum(); truth = yA.sum() / totA
    def dd(P):
        lo = hi = 0.0
        for u in range(12):
            m = P == u; nr, nc, tot = gA[m].sum(), (yA[m] + yB[m]).sum(), (gA[m] + gB[m]).sum()
            l, h = fh(nr, nc, tot); lo += l; hi += h
        return lo / totA, hi / totA
    d1, d2 = dd(P1), dd(P2)
    # joint LP: variables yA_a, yB_a (yes counts by group in each atom), 0<=yA_a<=gA_a etc
    Aeq, beq = [], []
    for P in (P1, P2):
        for u in range(12):
            m = P == u; row = np.zeros(2 * na); row[np.where(m)[0]] = 1; row[na + np.where(m)[0]] = 1
            Aeq.append(row); beq.append((yA[m] + yB[m]).sum())
    bnds = [(0, g) for g in gA] + [(0, g) for g in gB]
    c = np.concatenate([np.ones(na), np.zeros(na)]) / totA
    lo = linprog(c, A_eq=Aeq, b_eq=beq, bounds=bnds, method="highs").fun
    hi = -linprog(-c, A_eq=Aeq, b_eq=beq, bounds=bnds, method="highs").fun
    out.append(f"| {lab} | [{d1[0]:.3f}, {d1[1]:.3f}] | [{d2[0]:.3f}, {d2[1]:.3f}] | [{max(d1[0], d2[0]):.3f}, {min(d1[1], d2[1]):.3f}] | [{lo:.3f}, {hi:.3f}] | {truth:.3f} |")
out.append("Atoms (e.g. census blocks) have KNOWN group counts; outcomes are observed only as unit totals, in both partitions. DD per partition uses each unit's group and outcome totals.")
out.append("")

# ---------- (iii)
out.append("## (iii) age-stratified margins inside units (2 groups x yes/no x 4 age bands), 200 units")
wid_plain = wid_strat = 0.0; tot = 0.0
for u in range(200):
    nAk = rng.integers(5, 80, 4); nBk = rng.integers(5, 80, 4)
    age_eff = np.array([-0.25, -0.05, 0.1, 0.3]); base = rng.uniform(0.3, 0.7)
    yAk = rng.binomial(nAk, np.clip(base + 0.15 + age_eff, 0, 1)); yBk = rng.binomial(nBk, np.clip(base - 0.15 + age_eff, 0, 1))
    l, h = fh(nAk.sum(), (yAk + yBk).sum(), (nAk + nBk).sum()); wid_plain += h - l
    for k in range(4):
        l, h = fh(nAk[k], yAk[k] + yBk[k], nAk[k] + nBk[k]); wid_strat += h - l
    tot += nAk.sum()
out.append(f"DD width for group A's aggregate yes rate: one-way margins {wid_plain / tot:.3f}; with group x age and outcome x age margins {wid_strat / tot:.3f} (ratio {wid_strat / wid_plain:.2f}).")
out.append("")

# ---------- (iv)
if len(sys.argv) > 1:
    import os, pandas as pd
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from apply_B_suffrage import prep_wi
    wi = pd.read_csv(sys.argv[1])
    out.append("## (iv) Wisconsin: covariance of unit-level yes rates from Var(T|x) = x' Sigma x (units >= 100 men)")
    out.append("| year | sd veterans | sd non-veterans | correlation | PSD | cluster-bootstrap 90% interval for the correlation |")
    out.append("|---|---|---|---|---|---|")
    for y in ("1857", "1865"):
        d = prep_wi(wi, y); d = d[d.N >= 100].reset_index(drop=True)
        cl = wi.loc[wi.County_Cluster.isin(d.county), "County_Cluster"]
        def est(dd):
            x = (dd.V / dd.N).values; T = (dd.yes / dd.N).values; n = dd.N.values
            X = np.column_stack([np.ones_like(x), x, x * x]); r = T - X @ np.linalg.lstsq(X, T, rcond=None)[0]
            r2 = r ** 2 - T * (1 - T) / n                              # remove (an upper bound on) coin-flip noise
            Z = np.column_stack([x * x, 2 * x * (1 - x), (1 - x) ** 2])
            s = np.linalg.lstsq(Z, r2, rcond=None)[0]
            S = np.array([[s[0], s[1]], [s[1], s[2]]])
            return S
        S = est(d)
        cor = S[0, 1] / np.sqrt(max(S[0, 0], 1e-12) * max(S[1, 1], 1e-12))
        cors = []
        counties = d.county.unique()
        for _ in range(500):
            pick = rng.choice(counties, len(counties))
            db = pd.concat([d[d.county == c] for c in pick])
            Sb = est(db)
            if Sb[0, 0] > 0 and Sb[1, 1] > 0: cors.append(Sb[0, 1] / np.sqrt(Sb[0, 0] * Sb[1, 1]))
        lo, hi = np.quantile(cors, [0.05, 0.95])
        psd = bool(np.all(np.linalg.eigvalsh(S) >= -1e-12))
        sdv = np.sqrt(S[0, 0]) if S[0, 0] > 0 else float("nan"); sdn = np.sqrt(S[1, 1]) if S[1, 1] > 0 else float("nan")
        out.append(f"| {y} | {sdv:.3f} | {sdn:.3f} | {cor:.2f} | {psd} | [{lo:.2f}, {hi:.2f}] |")
txt = "\n".join(out); print(txt); open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
