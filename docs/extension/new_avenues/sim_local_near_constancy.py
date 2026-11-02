"""
Candidate assumption 4: local near-constancy with a data-implied floor (2x2).
Units i have composition x_i (share of group W) and outcome share Y_i = x_i bW_i + (1-x_i) bB_i.
Neighbours: the k nearest units in composition (chosen without looking at outcomes).
  sup version : |bW_i - bW_j| <= eps and |bB_i - bB_j| <= eps for every neighbour pair
  mean version: average over neighbour pairs of (|dW| + |dB|)/2 <= eps   (robust to a few odd pairs)
eps_min = smallest eps with a feasible LP (data-implied floor).  Target: aggregate rate of group W.
Rates here are the units' realised shares (large units, so no binomial noise) unless stated.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix
rng = np.random.default_rng(4)

def make(p, s, n_mean=None):
    x = np.sort(rng.uniform(0.05, 0.95, p))
    bW = np.clip(0.20 + 0.30 * x + s * rng.standard_normal(p), 0, 1)
    bB = np.clip(0.70 - 0.20 * x + s * rng.standard_normal(p), 0, 1)
    n = np.full(p, 1000)
    if n_mean is not None:                      # binomial realisation of the shares
        n = rng.integers(n_mean // 2, 3 * n_mean // 2, p)
        nW = np.round(n * x).astype(int); x = nW / n
        bW = rng.binomial(nW, bW) / np.maximum(nW, 1); bB = rng.binomial(n - nW, bB) / np.maximum(n - nW, 1)
    y = x * bW + (1 - x) * bB
    truth = np.sum(n * x * bW) / np.sum(n * x)
    make.last = (bW, bB)
    return x, y, n, truth

def pairs_knn(x, k):
    P = set()
    for i in range(len(x)):
        for j in np.argsort(np.abs(x - x[i]))[1:k + 1]:
            P.add((min(i, j), max(i, j)))
    return sorted(P)

def solve(x, y, n, P, eps=None, mode="sup", target=None):
    """Variables: bW (p), bB (p), [eps], and for mode='mean' slack t_W, t_B per pair."""
    p, m = len(x), len(P)
    nv = 2 * p + 1 + (2 * m if mode == "mean" else 0)
    ie = 2 * p
    Aeq = lil_matrix((p, nv)); beq = y.copy()
    for i in range(p):
        Aeq[i, i] = x[i]; Aeq[i, p + i] = 1 - x[i]
    rows = []; b = []
    for q, (i, j) in enumerate(P):
        for off in (0, p):
            if mode == "sup":
                rows.append({off + i: 1, off + j: -1, ie: -1}); b.append(0)
                rows.append({off + i: -1, off + j: 1, ie: -1}); b.append(0)
            else:
                t = 2 * p + 1 + 2 * q + (0 if off == 0 else 1)
                rows.append({off + i: 1, off + j: -1, t: -1}); b.append(0)
                rows.append({off + i: -1, off + j: 1, t: -1}); b.append(0)
    if mode == "mean":
        r = {2 * p + 1 + t: 1.0 / (2 * m) for t in range(2 * m)}; r[ie] = -1; rows.append(r); b.append(0)
    A = lil_matrix((len(rows), nv))
    for k, r in enumerate(rows):
        for c, v in r.items(): A[k, c] = v
    bounds = [(0, 1)] * (2 * p) + [((0, None) if eps is None else (eps, eps))] + [(0, None)] * (nv - 2 * p - 1)
    c = np.zeros(nv)
    if target is None:
        c[ie] = 1
    else:
        w = n * x / np.sum(n * x); c[:p] = w * target
    res = linprog(c, A_ub=csr_matrix(A), b_ub=np.array(b), A_eq=csr_matrix(Aeq), b_eq=beq, bounds=bounds, method="highs")
    if res.status != 0: return None
    return res.fun if target is None else target * res.fun

out = ["# Local near-constancy with a data-implied floor (2x2)", ""]
for lab, s, nm in (("low heterogeneity (sd 0.03), large units", 0.03, None),
                   ("high heterogeneity (sd 0.08), large units", 0.08, None),
                   ("low heterogeneity (sd 0.03), units of ~200 (binomial noise)", 0.03, 200)):
    x, y, n, truth = make(200, s, nm)
    P = pairs_knn(x, 5)
    w = n * x / np.sum(n * x)
    dd = (np.sum(w * np.maximum(0, (y - (1 - x)) / x)), np.sum(w * np.minimum(1, y / x)))
    out.append(f"## {lab}: 200 units, 5 nearest-composition neighbours; truth {truth:.3f}; DD [{dd[0]:.3f}, {dd[1]:.3f}]")
    out.append("| version | eps_min | eps | bound on group W's aggregate rate | width / DD | covers |")
    out.append("|---|---|---|---|---|---|")
    for mode in ("sup", "mean"):
        emin = solve(x, y, n, P, None, mode)
        for mult in (1.0, 1.5, 2.0, 3.0):
            e = emin * mult
            lo = solve(x, y, n, P, e, mode, target=1.0); hi = solve(x, y, n, P, e, mode, target=-1.0)
            if lo is None or hi is None:
                out.append(f"| {mode} | {emin:.3f} | {mult}x | infeasible | | |"); continue
            out.append(f"| {mode} | {emin:.3f} | {mult}x = {e:.3f} | [{lo:.3f}, {hi:.3f}] | {(hi-lo)/(dd[1]-dd[0]):.2f} | {lo-1e-7 <= truth <= hi+1e-7} |")
    out.append("")
txt = "\n".join(out); print(txt); open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")

# ---- All-pairs version: |b_i - b_j| <= e + L |x_i - x_j| for ALL pairs (no chain drift).
def solve_all(x, y, n, e, L, target):
    p = len(x); iu, ju = np.triu_indices(p, 1); m = len(iu)
    A = lil_matrix((4 * m, 2 * p)); b = np.empty(4 * m); r = 0
    for off in (0, p):
        for s in (1, -1):
            A[np.arange(r, r + m), off + iu] = s; A[np.arange(r, r + m), off + ju] = -s
            b[r:r + m] = e + L * np.abs(x[iu] - x[ju]); r += m
    Aeq = lil_matrix((p, 2 * p)); Aeq[np.arange(p), np.arange(p)] = x; Aeq[np.arange(p), p + np.arange(p)] = 1 - x
    w = n * x / np.sum(n * x); c = np.concatenate([w * target, np.zeros(p)])
    res = linprog(c, A_ub=csr_matrix(A), b_ub=b, A_eq=csr_matrix(Aeq), b_eq=y, bounds=[(0, 1)] * (2 * p), method="highs")
    return None if res.status != 0 else target * res.fun

def emin_all(x, y, n, L, hi=1.0):
    lo = 0.0
    for _ in range(25):
        mid = (lo + hi) / 2
        if solve_all(x, y, n, mid, L, 1.0) is None: lo = mid
        else: hi = mid
    return hi

out2 = ["", "## All-pairs version: |b_i - b_j| <= e + L|x_i - x_j| for every pair (100 units)",
        "Truth: group slopes +0.30 and -0.20 in x, unit heterogeneity sd s.  e_min(L) = smallest feasible e.",
        "| design | L | e_min(L) | e used | bound on W's aggregate rate | width / DD | covers |", "|---|---|---|---|---|---|---|"]
for lab, s, nm in (("sd 0.03, large units", 0.03, None), ("sd 0.08, large units", 0.08, None),
                   ("sd 0.03, units ~200", 0.03, 200)):
    x, y, n, truth = make(100, s, nm)
    w = n * x / np.sum(n * x)
    dd = (np.sum(w * np.maximum(0, (y - (1 - x)) / x)), np.sum(w * np.minimum(1, y / x)))
    for L in (0.0, 0.3, 0.6):
        em = emin_all(x, y, n, L)
        bW, bB = make.last; iu, ju = np.triu_indices(len(x), 1)
        etrue = max(np.max(np.abs(bW[iu] - bW[ju]) - L * np.abs(x[iu] - x[ju])), np.max(np.abs(bB[iu] - bB[ju]) - L * np.abs(x[iu] - x[ju])))
        out2.append(f"| {lab}: e the truth needs at L={L} = {etrue:.3f} (rows with e below this assume something false) | | | | | | |")
        for mult in (1.0, 1.5, 2.0):
            e = em * mult
            lo, hi = solve_all(x, y, n, e, L, 1.0), solve_all(x, y, n, e, L, -1.0)
            if lo is None or hi is None:
                out2.append(f"| {lab} | {L} | {em:.3f} | {mult}x | infeasible | | |"); continue
            out2.append(f"| {lab} (truth {truth:.3f}) | {L} | {em:.3f} | {mult}x = {e:.3f} | [{lo:.3f}, {hi:.3f}] | {(hi-lo)/(dd[1]-dd[0]):.2f} | {lo-1e-7 <= truth <= hi+1e-7} |")
txt2 = "\n".join(out2); print(txt2)
open(__file__.replace(".py", "_output.md"), "a").write(txt2 + "\n")
