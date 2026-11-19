"""
Exploratory diagnostics on the ei.Datasets New Zealand tables (3x3 only; 450 tables,
about 31,000 polling places).  Input: CSVs written by export_nz_units.R.
Usage: python3 nz_diagnostics.py <dir with nz_units.csv and nz_truth.csv>

D1 (Avenue A feasibility).  Under the binomial-propensity model,
      Var(T_c | x, n) = A_c(x) + B_c(x)/n,   0 <= B_c(x) <= sum_r x_r m_rc (1 - m_rc).
    We regress leverage-corrected squared residuals of the within-table quadratic fit
    on table fixed effects, a global quadratic in x (for A), and b_i/n_i with
    b_i = sum_r x_ir B_rc (1 - B_rc) built from the TRUE aggregate rates.  Binomial
    voting with little poll-level heterogeneity predicts a slope a bit below 1; a slope
    well above 1 signals within-poll dependence (over-dispersion), which would make
    the unit-size extrapolation invalid without a design-effect allowance.
D2 (Avenue E feasibility).  For each table, the part of the true rate that M0 must
    attribute to K:  v_rc = B_rc(truth) - sum_s Qhat_c,rs mu_rs / mu_r.  We report its
    size and its correlation across consecutive elections in the same electorate name
    (exploratory: v is a pseudo-true, model-dependent quantity and boundaries change).
"""
import sys
import numpy as np, pandas as pd
d = sys.argv[1] if len(sys.argv) > 1 else "."
U = pd.read_csv(f"{d}/nz_units.csv"); T = pd.read_csv(f"{d}/nz_truth.csv")
rows = ["LAB", "NAT", "OTH"]
U["n"] = U[[f"x_{r}" for r in rows]].sum(axis=1)
U["nt"] = U[[f"t_{r}" for r in rows]].sum(axis=1)
U = U[(U.n > 0) & (U.nt > 0)].copy()
for r in rows:
    U[f"X_{r}"] = U[f"x_{r}"] / U.n; U[f"T_{r}"] = U[f"t_{r}"] / U.nt
out = ["# NZ exploratory diagnostics", "",
       f"tables: {U.table.nunique()}, polling places: {len(U)}, median size {U.n.median():.0f}, "
       f"10%/90% sizes {U.n.quantile(.1):.0f}/{U.n.quantile(.9):.0f}, share with n<100: {(U.n<100).mean():.2f}",
       f"share of units whose party-vote and candidate-vote totals differ by >2%: {(abs(U.n-U.nt)/U.n>0.02).mean():.2f}", ""]

def quad_design(X):
    a, b, c = X[:, 0], X[:, 1], X[:, 2]
    return np.column_stack([a * a, b * b, c * c, a * b, a * c, b * c])

truthB = {}
for k, g in T.groupby("table"):
    M = g.pivot(index="row", columns="col", values="count").loc[rows, rows].values.astype(float)
    truthB[k] = M / M.sum(1, keepdims=True)

recs, vrecs = [], []
for k, g in U.groupby("table"):
    if len(g) < 15: continue
    X = g[[f"X_{r}" for r in rows]].values; Tm = g[[f"T_{r}" for r in rows]].values
    n = g.n.values.astype(float); Z = quad_design(X)
    H = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T; lev = np.clip(np.diag(H), 0, 0.9)
    coef = np.linalg.lstsq(Z, Tm, rcond=None)[0]           # 6 x 3
    res = Tm - Z @ coef
    B = truthB[k]
    for c in range(3):
        b = X @ (B[:, c] * (1 - B[:, c]))
        for i in range(len(g)):
            recs.append((k, c, res[i, c] ** 2 / (1 - lev[i]), b[i] / n[i], *X[i]))
    # D2: Q from coefficients. x'Qx with Q sym: diag = coef[0..2], off = coef[3..5]/2
    mu = np.einsum("i,ir,is->rs", n, X, X); mur = mu.sum(1)
    for c in range(3):
        Q = np.diag(coef[:3, c]); Q[0, 1] = Q[1, 0] = coef[3, c] / 2; Q[0, 2] = Q[2, 0] = coef[4, c] / 2; Q[1, 2] = Q[2, 1] = coef[5, c] / 2
        for r in range(3):
            v = B[r, c] - (Q[r] @ mu[r]) / mur[r]
            vrecs.append((k, g.election.iloc[0], g.district.iloc[0], rows[r], rows[c], v))

R = pd.DataFrame(recs, columns=["table", "c", "r2", "b_over_n", "x1", "x2", "x3"])
out.append("## D1: is the 1/n (binomial) term visible, and is it binomial-sized?")
out.append("| column | slope on b/n | cluster SE | units | reading |")
out.append("|---|---|---|---|---|")
for c in range(3):
    s = R[R.c == c].copy()
    # within-table demeaning (table FE) of y and regressors
    regs = ["b_over_n", "x1", "x2", "x1x1", "x2x2", "x1x2"]
    s["x1x1"], s["x2x2"], s["x1x2"] = s.x1 ** 2, s.x2 ** 2, s.x1 * s.x2
    cols = ["r2"] + regs
    dm = s[cols] - s.groupby("table")[cols].transform("mean")
    Xr = dm[regs].values; y = dm.r2.values
    XtXi = np.linalg.pinv(Xr.T @ Xr); beta = XtXi @ Xr.T @ y; e = y - Xr @ beta
    meat = np.zeros((len(regs), len(regs)))
    for _, idx in s.groupby("table").indices.items():
        sc = Xr[idx].T @ e[idx]; meat += np.outer(sc, sc)
    se = np.sqrt(np.diag(XtXi @ meat @ XtXi))
    reading = "binomial-sized" if beta[0] - 2 * se[0] <= 1 <= beta[0] + 2 * se[0] or beta[0] < 1 else "over-dispersed (>1)"
    out.append(f"| {rows[c]} | {beta[0]:.3f} | {se[0]:.3f} | {len(s)} | {reading} |")
out.append("")

V = pd.DataFrame(vrecs, columns=["table", "election", "district", "row", "col", "v"])
out.append("## D2: size and persistence of the part of the truth that M0 attributes to K")
out.append("| cell (row->col) | median abs v | corr with same electorate, previous election | pairs |")
out.append("|---|---|---|---|")
elections = sorted(V.election.unique())
for (r, c), g in V.groupby(["row", "col"]):
    g = g.set_index(["district", "election"]).v
    pairs = []
    for (dist, el), val in g.items():
        j = elections.index(el)
        if j > 0 and (dist, elections[j - 1]) in g.index:
            pairs.append((g[(dist, elections[j - 1])], val))
    P = np.array(pairs)
    cor = np.corrcoef(P[:, 0], P[:, 1])[0, 1] if len(P) > 3 else np.nan
    out.append(f"| {r}->{c} | {g.abs().median():.3f} | {cor:.2f} | {len(P)} |")
txt = "\n".join(out); print(txt)
open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
