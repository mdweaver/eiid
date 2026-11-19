"""
Low-rank (single-index) contextual effects.  Under LCE the stacked mean tables vec(M_H(x_i)) are
linear in x_i, so the units x (R*C) matrix of means has rank <= R; rank 1 is Goodman.
Intermediate rung: rank 2, i.e. H_j = m_j 1' + d_j c'  (rates depend on composition only through
one index c'x).  Numerical Jacobian rank of the map (m, d, c) -> Q = sym(H) at a random point
gives the number of directions the regression leaves unidentified (generic local count).
"""
import numpy as np
rng = np.random.default_rng(1)
def counts(R, C):
    one = np.ones(R)
    B = np.linalg.svd(np.eye(R) - np.outer(one, one) / R)[0][:, :R - 1]
    def Qmap(th):
        i = 0; ms, ds = [], []
        for j in range(C - 1):
            ms.append(th[i:i + R]); i += R; ds.append(th[i:i + R]); i += R
        cdir = B @ np.concatenate([[1.0], th[i + 1:i + R - 1]]) if R > 2 else B[:, 0]
        out = []
        for m, d in zip(ms, ds):
            H = np.outer(m, one) + np.outer(d, cdir); Q = (H + H.T) / 2
            out += list(Q[np.triu_indices(R)])
        return np.array(out)
    npar = 2 * R * (C - 1) + (R - 1)
    th = rng.normal(size=npar); eps = 1e-6
    J = np.array([(Qmap(th + eps * e) - Qmap(th - eps * e)) / (2 * eps) for e in np.eye(npar)]).T
    rk = np.linalg.matrix_rank(J, tol=1e-6)
    eff = npar - 1                      # one coordinate is a placeholder (scale of c fixed)
    return eff, rk, J.shape[0]
lines = ["| R | C | single-index params | identified equations (Q) | rank | unidentified (single index) | unidentified (LCE) | testable restrictions |", "|---|---|---|---|---|---|---|---|"]
for R, C in [(2, 2), (2, 3), (3, 2), (3, 3), (4, 2), (4, 3), (5, 3)]:
    p, rk, neq = counts(R, C)
    lines.append(f"| {R} | {C} | {p} | {neq} | {rk} | {p - rk} | {(C - 1) * R * (R - 1) // 2} | {neq - rk} |")
txt = "\n".join(lines); print(txt); open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
