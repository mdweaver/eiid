"""Reference numbers for the M-SIGN module tests (package/eiPartialID/tests/testthat).

Writes three CSV files next to the R tests:
  msign_example_3x2.csv     noise-free 3x2 design: X1..X3, N, T1 = tau(X), T2 = 1 - T1
  msign_reference_3x2.csv   affine-in-k forms (c, a1, a2, a3) of delta_B / delta_W (W1, W2,
                            pairwise), ER / NM / gamma, D, and LP intervals for B~_r under
                            LCE alone and under CR (branches +, -), from verify_R3.py's
                            assoc_3x2 / assoc_alt
  msign_gap_3x2.csv         a 3x2 design (found with the verify_R3 S5 search) where CR with
                            a non-W1 pairing has two disjoint branch images; its reference rows
                            are appended to msign_reference_3x2.csv with case = "gap"
The k-parametrisation K[r,s] = k_(r,s), pairs (0,1),(0,2),(1,2), is the package's z for
column 1 (K_1[r,s] = z, rows/pairs in the same order).
numpy/scipy only.  Run: PYTHONDONTWRITEBYTECODE=1 python3 make_msign_reference.py
"""
import os
import sys

import numpy as np
from scipy.optimize import linprog

sys.dont_write_bytecode = True
import verify_package_design as V  # noqa: E402
import verify_R3 as R3  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "..", "package", "eiPartialID", "tests", "testthat")
FREE = [(None, None)] * 3


def lp_interval(a, c, A, b):
    out = []
    for s in (1, -1):
        res = linprog(s * a, A_ub=A, b_ub=b, bounds=FREE, method="highs")
        if res.status == 2:
            return None
        out.append(s * res.fun + c)
    return tuple(out)


def tau(P):
    return np.einsum("ij,jk,ik->i", P["Xs"], P["Q"], P["Xs"])


def descriptive(P, r):
    Xs, p = P["Xs"], P["p"]
    mu = p @ Xs
    t = tau(P)
    var = p @ (Xs[:, r] - mu[r]) ** 2
    cov = p @ ((Xs[:, r] - mu[r]) * t)
    ER = cov / var
    NM = cov / (mu[r] * (1 - mu[r]))
    gam = var / (mu[r] * (1 - mu[r]))
    return ER, NM, gam


def D_form(P, r):
    Xs, p = P["Xs"], P["p"]
    mu = p @ Xs
    a = np.zeros(3); c = 0.0
    for s in range(3):
        ta, tc = V.target_3x2(P, s)
        wgt = 1.0 if s == r else -mu[s] / (1 - mu[r])
        a = a + wgt * ta[0]; c = c + wgt * tc[0]
    return a, c


def forms(P):
    rows = []
    for r in range(3):
        (aB, cB), (aW, cW) = R3.assoc_3x2(P, r)
        rows.append(("dB_W1", r, -1, cB[0], *aB[0]))
        rows.append(("dW_W1", r, -1, cW[0], *aW[0]))
        _, (a2, c2) = R3.assoc_alt(P, r, (r + 1) % 3)
        rows.append(("dW_W2", r, -1, c2[0], *a2[0]))
        ER, NM, gam = descriptive(P, r)
        rows.append(("ER", r, -1, ER, 0, 0, 0))
        rows.append(("NM", r, -1, NM, 0, 0, 0))
        rows.append(("gamma", r, -1, gam, 0, 0, 0))
        aD, cD = D_form(P, r)
        rows.append(("D", r, -1, cD, *aD))
        ta, tc = V.target_3x2(P, r)
        rows.append(("Btilde", r, -1, tc[0], *ta[0]))
    for r, s in V.PAIRS:
        (aBp, cBp), _ = R3.assoc_alt(P, r, s)
        rows.append(("dB_pair", r, s, cBp[0], *aBp[0]))
    return rows


def pairing(P, version, r, other):
    (aB1, cB1), (aW1, cW1) = R3.assoc_3x2(P, r)
    (aBp, cBp), (aW2, cW2) = R3.assoc_alt(P, r, other)
    if version == "W1":
        return (aB1, cB1), (aW1, cW1)
    if version == "W2":
        return (aB1, cB1), (aW2, cW2)
    return (aBp, cBp), (aW2, cW2)


def intervals(P, which, version, r, other):
    """LP interval of B~_t under LCE alone and under the two CR branches (rows-sign +1/-1)."""
    A, b = V.rows_3x2(P, which)
    (aB, cB), (aW, cW) = pairing(P, version, r, other)
    out = []
    for t in range(3):
        a, c = V.target_3x2(P, t)
        I0 = lp_interval(a[0], c[0], A, b)
        Ip = lp_interval(a[0], c[0], np.vstack([A, -aB, -aW]), np.concatenate([b, cB, cW]))
        Im = lp_interval(a[0], c[0], np.vstack([A, aB, aW]), np.concatenate([b, -cB, -cW]))
        res = linprog(np.zeros(3), A_ub=A, b_ub=b, A_eq=np.vstack([aB, aW]),
                      b_eq=-np.concatenate([cB, cW]), bounds=FREE, method="highs")
        nan2 = (np.nan, np.nan)
        out.append((f"LCE_{version}", t, -1, *(I0 or nan2), 0, 0))
        out.append((f"CRplus_{version}", t, -1, *(Ip or nan2), 0, 0))
        out.append((f"CRminus_{version}", t, -1, *(Im or nan2), 0, 0))
        out.append((f"cert_{version}", t, -1, float(res.status == 0), 0, 0, 0))
    return out


def write_data(P, path):
    Xs, t = P["Xs"], tau(P)
    N = np.round(P["p"] * P["Ntot"]).astype(int)
    with open(path, "w") as f:
        f.write("X1,X2,X3,N,T1,T2\n")
        for i in range(len(t)):
            f.write(",".join(repr(float(v)) for v in Xs[i]) + f",{N[i]},{float(t[i])!r},{float(1 - t[i])!r}\n")


def with_integer_N(P):
    # the package takes integer-like N; re-derive p from rounded N so both sides agree
    N = np.round(P["p"] * 1e6)
    P = dict(P); P["p"] = N / N.sum(); P["Ntot"] = N.sum()
    return P


if __name__ == "__main__":
    RNG = np.random.default_rng(20261105)
    V.RNG = RNG
    # --- fixed example ------------------------------------------------------
    P = with_integer_N(V.make_3x2(n=60))
    write_data(P, os.path.join(OUT, "msign_example_3x2.csv"))
    rows = [("main",) + r for r in forms(P)]
    for ver in ("W1", "W2", "pair"):
        rows += [("main",) + r for r in intervals(P, ["DD"], ver, 0, 1)]
    # --- gap search (verify_R3 S5 approach: 2-d / near-1-d support, DD or DD + ordering) ---
    found = None
    for it in range(6000):
        support = ("2d", "near-1d")[it % 2]
        which = (["DD"], ["DD", "O"])[(it // 2) % 2]
        Pg = with_integer_N(V.make_3x2(n=60, support=support))
        A, b = V.rows_3x2(Pg, which)
        for ver in ("W2", "pair"):
            (aB, cB), (aW, cW) = pairing(Pg, ver, 0, 1)
            for t in range(3):
                a, c = V.target_3x2(Pg, t)
                Ip = lp_interval(a[0], c[0], np.vstack([A, -aB, -aW]), np.concatenate([b, cB, cW]))
                Im = lp_interval(a[0], c[0], np.vstack([A, aB, aW]), np.concatenate([b, -cB, -cW]))
                if Ip is None or Im is None:
                    continue
                lo, hi = sorted([Ip, Im])
                gap = hi[0] - lo[1]
                hull = max(Ip[1], Im[1]) - min(Ip[0], Im[0])
                if gap > 0.01 * hull:
                    found = (Pg, ver, t, gap, hull)
                    break
            if found:
                break
        if found:
            break
    if found is None:
        raise SystemExit("no gap found")
    Pg, ver, t, gap, hull = found
    print(f"gap case: iteration {it}, support {support}, rows {which}, pairing {ver}, "
          f"target B~_{t}, gap {gap:.4f}, hull {hull:.4f}")
    # ordering rows used (package: restrictions$order, column 1), as reference rows
    for r, s in V.PAIRS:
        g = Pg["Xs"] @ (Pg["H"][r] - Pg["H"][s])
        if "O" in which and (np.all(g >= 0) or np.all(g <= 0)):
            hi, lo = (r, s) if np.all(g >= 0) else (s, r)
            rows.append(("gap", "order", hi, lo, 0, 0, 0, 0))
    write_data(Pg, os.path.join(OUT, "msign_gap_3x2.csv"))
    rows += [("gap",) + r for r in forms(Pg)]
    for v2 in ("W1", "W2", "pair"):
        rows += [("gap",) + r for r in intervals(Pg, which, v2, 0, 1)]
    with open(os.path.join(OUT, "msign_reference_3x2.csv"), "w") as f:
        f.write("case,name,r,s,c,a1,a2,a3\n")
        for row in rows:
            f.write(f"{row[0]},{row[1]},{row[2] + 1},{row[3] + 1 if row[3] >= 0 else 0},"
                    + ",".join(repr(float(v)) for v in row[4:]) + "\n")
    print(f"wrote {len(rows)} reference rows")
