"""Checks for proof-register item R2 (PHASE5a_sign_module.md §2): when does an
added restriction tighten a target interval?

  T1  face test vs brute force (3x2, LCE base, random targets; random rows and
      degenerate rows passing through an optimal vertex).  Also counts how
      often the solver's dual marginal is > 0 although nothing tightens.
  T2  'implied' (max_F a'k <= b)  <=>  tightens no linear target; witness
      target phi = a'k tightens whenever the row is not implied.
  T3  2x2: sign rows on delta*_W, delta*_B either are implied or tighten D;
      closed forms [L_J, min(U_J, D_ER)] etc. match the LP.
numpy/scipy only.  Seed 20261105.
"""
import sys

import numpy as np
from scipy.optimize import linprog

sys.dont_write_bytecode = True
import verify_package_design as V  # noqa: E402

RNG = np.random.default_rng(20261105)
V.RNG = RNG
FREE = [(None, None)] * 3


def lpmax(c, A, b):
    r = linprog(-c, A_ub=A, b_ub=b, bounds=FREE, method="highs")
    return (None, None) if r.status != 0 else (-r.fun, r)


def face_min(a, gam, U, A, b, eps=1e-9):
    """min a'k over the optimal face {k in F : gam'k >= U - eps}."""
    r = linprog(a, A_ub=np.vstack([A, -gam]), b_ub=np.append(b, -(U - eps)),
                bounds=FREE, method="highs")
    return r.fun


def T1(reps=400):
    agree = n = active_redundant = dual_pos_redundant = 0
    for i in range(reps):
        P = V.make_3x2(support=("2d", "near-1d")[i % 2])
        A, b = V.rows_3x2(P, ["DD"])
        gam = RNG.normal(size=3)
        U, res = lpmax(gam, A, b)
        if U is None:
            continue
        kstar = res.x
        for kind in ("random", "through-vertex"):
            a = RNG.normal(size=3)
            if kind == "random":
                lo = linprog(a, A_ub=A, b_ub=b, bounds=FREE, method="highs").fun
                hi = -linprog(-a, A_ub=A, b_ub=b, bounds=FREE, method="highs").fun
                bb = RNG.uniform(lo, hi)
            else:
                bb = a @ kstar                       # row active at the solver's optimum
            A2, b2 = np.vstack([A, a]), np.append(b, bb)
            U2, res2 = lpmax(gam, A2, b2)
            brute = (U - U2) > 1e-7 * max(1, abs(U))
            test = face_min(a, gam, U, A, b) > bb + 1e-7 * max(1, abs(bb))
            agree += brute == test
            n += 1
            lam = -res2.ineqlin.marginals[-1]        # HiGHS: marginals <= 0 for max via min(-c)
            if not brute and abs(a @ res2.x - bb) < 1e-8:
                active_redundant += 1
                dual_pos_redundant += lam > 1e-9
    return n, agree, active_redundant, dual_pos_redundant


def T2(reps=400):
    ok = n = 0
    for _ in range(reps):
        P = V.make_3x2()
        A, b = V.rows_3x2(P, ["DD"])
        a = RNG.normal(size=3)
        hi, _ = lpmax(a, A, b)
        lo = -lpmax(-a, A, b)[0]
        bb = RNG.uniform(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
        if bb < lo:                                   # would be empty; skip
            continue
        implied = hi <= bb
        U2, _ = lpmax(a, np.vstack([A, a]), np.append(b, bb))
        witness_tightens = (hi - U2) > 1e-9
        ok += implied != witness_tightens
        n += 1
    return n, ok


def T3(reps=2000):
    agree = n = middle = 0
    worst = 0.0
    for _ in range(reps):
        P = V.make_2x2()
        J = V.jkst_rows(P)
        IJ = V.interval(J)
        LJ, UJ = V.map_interval(P, IJ, "D")
        f = V.funcs_2x2(P, 0.0)
        DER, DNM = f["D_ER"], f["D_NM"]
        aW, bW = V.affine(lambda w: V.funcs_2x2(P, w)["dW"])
        aB, bB = V.affine(lambda w: V.funcs_2x2(P, w)["dB"])
        for (a, c), closed, implied in (
            ((aW, bW), (LJ, min(UJ, DER)), DER >= UJ),        # delta_W >= 0
            ((-aW, -bW), (max(LJ, DER), UJ), DER <= LJ),      # delta_W <= 0
            ((aB, bB), (max(LJ, DNM), UJ), DNM <= LJ),        # delta_B >= 0
            ((-aB, -bB), (LJ, min(UJ, DNM)), DNM >= UJ),      # delta_B <= 0
        ):
            I = V.interval(J + [(a, c)])
            if I is None:
                assert closed[0] > closed[1] - 1e-9
                continue
            lp = V.map_interval(P, I, "D")
            worst = max(worst, abs(lp[0] - closed[0]), abs(lp[1] - closed[1]))
            tightens = (UJ - LJ) - (lp[1] - lp[0]) > 1e-9
            # in 2x2: not implied <=> tightens (no 'cuts set but not target' case)
            sl, _ = V.affine(lambda w: V.funcs_2x2(P, w)["D"])
            if abs(sl) > 1e-12:
                agree += tightens == (not implied)
                n += 1
            middle += (not implied) and not tightens
    return n, agree, middle, worst


if __name__ == "__main__":
    n, agree, ar, dpr = T1()
    print(f"T1 face test = brute force in {agree}/{n} (3x2, random + through-vertex rows); "
          f"active-but-redundant rows {ar}, of which solver dual > 0 in {dpr}  "
          f"{'PASS' if agree == n else 'FAIL'}")
    n, ok = T2()
    print(f"T2 implied <=> witness target does not tighten: {ok}/{n}  {'PASS' if ok == n else 'FAIL'}")
    n, agree, mid, worst = T3()
    print(f"T3 2x2 sign rows: closed forms match LP (max err {worst:.1e}); not-implied <=> tightens D "
          f"in {agree}/{n}; 'cuts set but not D' cases {mid}  {'PASS' if agree == n and worst < 1e-8 else 'FAIL'}")


def T4(reps=300):
    """3x2, true-sign rows on delta^{0*}_W and delta^{0*}_B: per target B_0, B_1, B_2,
    classify as implied / cuts set only / tightens target."""
    import verify_R3 as R3
    counts = {"implied": 0, "cuts set, not this target": 0, "tightens target": 0}
    for i in range(reps):
        P = V.make_3x2(support=("2d", "near-1d")[i % 2])
        A, b = V.rows_3x2(P, ["DD"])
        (aB, cB), (aW, cW) = R3.assoc_3x2(P, 0)
        for ag, cg in ((aB[0], cB[0]), (aW[0], cW[0])):
            s = np.sign(ag @ P["k_true"] + cg)
            a, bb = -s * ag, s * cg                     # row a'k <= bb  <=>  s*g >= 0
            hi, _ = lpmax(a, A, b)
            if hi is None:
                continue
            for r in range(3):
                if hi <= bb + 1e-9:
                    counts["implied"] += 1
                    continue
                g, c = V.target_3x2(P, r)
                U, _ = lpmax(g[0], A, b); L = -lpmax(-g[0], A, b)[0]
                A2, b2 = np.vstack([A, a]), np.append(b, bb)
                U2, _ = lpmax(g[0], A2, b2); L2 = -lpmax(-g[0], A2, b2)[0]
                key = "tightens target" if (U - U2) + (L2 - L) > 1e-7 else "cuts set, not this target"
                counts[key] += 1
    return counts


if __name__ == "__main__":
    c = T4()
    tot = sum(c.values())
    print("T4 3x2 true delta^0 sign rows x targets B_0..B_2: " +
          "; ".join(f"{k} {v}/{tot}" for k, v in c.items()))
