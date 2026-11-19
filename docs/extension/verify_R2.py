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
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # sibling imports also under python3 -I
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
    skipped = {}
    for i in range(reps):
        P = V.make_3x2(support=("2d", "near-1d")[i % 2])
        A, b = V.rows_3x2(P, ["DD"])
        (aB, cB), (aW, cW) = R3.assoc_3x2(P, 0)
        for ag, cg in ((aB[0], cB[0]), (aW[0], cW[0])):
            s = np.sign(ag @ P["k_true"] + cg)
            a, bb = -s * ag, s * cg                     # row a'k <= bb  <=>  s*g >= 0
            hi, res_hi = lpmax(a, A, b)
            if hi is None:                              # record why, rather than skip silently
                st = linprog(-a, A_ub=A, b_ub=b, bounds=FREE, method="highs").status
                skipped[f"base solve status {st}"] = skipped.get(f"base solve status {st}", 0) + 3
                continue
            for r in range(3):
                if hi <= bb + 1e-9:
                    counts["implied"] += 1
                    continue
                g, c = V.target_3x2(P, r)
                A2, b2 = np.vstack([A, a]), np.append(b, bb)
                sols = [lpmax(g[0], A, b)[0], lpmax(-g[0], A, b)[0], lpmax(g[0], A2, b2)[0], lpmax(-g[0], A2, b2)[0]]
                if any(v is None for v in sols):
                    skipped["target solve failed"] = skipped.get("target solve failed", 0) + 1
                    continue
                U, L, U2, L2 = sols[0], -sols[1], sols[2], -sols[3]
                key = "tightens target" if (U - U2) + (L2 - L) > 1e-7 else "cuts set, not this target"
                counts[key] += 1
    return counts, skipped


def T1b(reps=150):
    """Non-singleton optimal faces (added after review).  Target normal = a facet
    normal A[j], so argmax_F A[j]'k is that facet (when the row is irredundant).
    Rows a'k <= bb with bb below / inside / above the range of a'k over the face."""
    out = {"below": [0, 0], "inside": [0, 0], "above": [0, 0]}
    faces = 0
    for i in range(reps):
        P = V.make_3x2(support=("2d", "near-1d")[i % 2])
        A, b = V.rows_3x2(P, ["DD"])
        for j in range(len(b)):
            gam = A[j]
            U, _ = lpmax(gam, A, b)
            if U is None or U < b[j] - 1e-9:          # row j redundant: no facet
                continue
            a = RNG.normal(size=3)
            Af, bf = np.vstack([A, -gam]), np.append(b, -(U - 1e-10))
            fmin = linprog(a, A_ub=Af, b_ub=bf, bounds=FREE, method="highs").fun
            fmax = -linprog(-a, A_ub=Af, b_ub=bf, bounds=FREE, method="highs").fun
            if fmax - fmin < 1e-6:                    # face is (numerically) a vertex
                continue
            faces += 1
            gmin = linprog(a, A_ub=A, b_ub=b, bounds=FREE, method="highs").fun
            cases = [("inside", fmin + RNG.uniform(0.1, 0.9) * (fmax - fmin)),
                     ("above", fmax + 0.1 * (fmax - fmin))]
            if fmin - gmin > 1e-6:                    # room to cut the whole face
                cases.append(("below", gmin + RNG.uniform(0.1, 0.9) * (fmin - gmin)))
            for kind, bb in cases:
                U2, _ = lpmax(gam, np.vstack([A, a]), np.append(b, bb))
                if U2 is None:
                    continue
                brute = (U - U2) > 1e-7 * max(1, abs(U))
                test = face_min(a, gam, U, A, b) > bb + 1e-7 * max(1, abs(bb))
                out[kind][0] += brute == test
                out[kind][1] += 1
    return faces, out


if __name__ == "__main__":
    faces, out = T1b()
    ok = all(a == n for a, n in out.values())
    print(f"T1b non-singleton optimal faces (facet-normal targets): {faces} faces; face test = brute force: " +
          "; ".join(f"row {k} face {a}/{n}" for k, (a, n) in out.items()) + f"  {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    c, sk = T4()
    tot = sum(c.values())
    print("T4 3x2 true delta^0 sign rows x targets B_0..B_2: " +
          "; ".join(f"{k} {v}/{tot}" for k, v in c.items()) +
          f"; not classified {sum(sk.values())}" + (f" ({'; '.join(f'{k}: {v}' for k, v in sk.items())})" if sk else "") +
          f"  (total {tot + sum(sk.values())} = 300 designs x 2 rows x 3 targets)")


# ---------------------------------------------------------------------------
# Added after the Batch B review (Nov 2026).
# ---------------------------------------------------------------------------
SQ = [(0, 1), (0, 1)]                                   # the unit square [0,1]^2


def _sq_max(c, A=None, b=None):
    kw = dict(bounds=SQ, method="highs")
    if A is not None:
        kw.update(A_ub=np.atleast_2d(A), b_ub=np.atleast_1d(b))
    r = linprog(-np.asarray(c, float), **kw)
    return None if r.status != 0 else -r.fun


def _sq_feasible(A, b):
    return linprog(np.zeros(2), A_ub=np.atleast_2d(A), b_ub=np.atleast_1d(b), bounds=SQ, method="highs").status == 0


def T1c():
    """R2 step 4 with a target intercept. F = [0,1]^2, phi = 2 + x, U = 3. The face LP must be
    {k in F, h'k >= U - phi0}; the old text wrote h'k >= U, here x >= 3, which is infeasible and so
    reports 'moves' for every added row. Rows: x + y <= 1.2 (keeps part of the facet x = 1) and
    x + y <= 0.8 (removes all of it)."""
    phi0, h = 2.0, np.array([1.0, 0.0])
    U = phi0 + _sq_max(h)
    rows = []
    for rhs in (1.2, 0.8):
        a = np.array([1.0, 1.0])
        brute_moves = (U - (phi0 + _sq_max(h, a, rhs))) > 1e-9
        old_moves = not _sq_feasible(np.vstack([-h, a]), [-U, rhs])            # h'k >= U
        new_moves = not _sq_feasible(np.vstack([-h, a]), [-(U - phi0), rhs])   # h'k >= U - phi0
        rows.append((rhs, brute_moves, old_moves, new_moves))
    old_ok = all(bm == om for _, bm, om, _ in rows)
    new_ok = all(bm == nm for _, bm, _, nm in rows)
    return rows, old_ok, new_ok


def T4b():
    """R2(e) wording. F = [0,1]^2, target x, rows N1: x + y <= 0.8 and N2: x <= 0.5. N1 added last
    has zero gain for x, yet it is not implied by F cap N2 (it removes (0.5, 1)) and it lowers the
    upper endpoint of x + y from 1.5 to 0.8. The old sentence 'redundant given the others' read as
    'implied' is false; 'zero incremental gain for this target' is true."""
    x, xy = np.array([1.0, 0.0]), np.array([1.0, 1.0])
    N1, N2 = (xy, 0.8), (x, 0.5)
    gain_x = _sq_max(x, [N2[0]], [N2[1]]) - _sq_max(x, [N2[0], N1[0]], [N2[1], N1[1]])
    gain_xy = _sq_max(xy, [N2[0]], [N2[1]]) - _sq_max(xy, [N2[0], N1[0]], [N2[1], N1[1]])
    implied = _sq_max(N1[0], [N2[0]], [N2[1]]) <= N1[1] + 1e-12      # R2(c) implication test
    old_claim = (abs(gain_x) < 1e-12) <= implied                    # 'zero last-added gain => redundant (implied)'
    new_claim = abs(gain_x) < 1e-12 and not implied and gain_xy > 0.5
    return gain_x, gain_xy, implied, old_claim, new_claim


if __name__ == "__main__":
    rows, old_ok, new_ok = T1c()
    print("T1c [F] target intercept (phi = 2 + x on [0,1]^2): " +
          "; ".join(f"row x+y<={r}: brute {'moves' if bm else 'keeps'}, old face LP {'moves' if om else 'keeps'}, "
                    f"corrected face LP {'moves' if nm else 'keeps'}" for r, bm, om, nm in rows) +
          f"; old step 4 agrees with brute force: {old_ok}  {'PASS' if new_ok and not old_ok else 'FAIL'}")
    gx, gxy, imp, oc, nc = T4b()
    print(f"T4b [F] target-specific redundancy: last-added gain for x {gx:.3f}, for x+y {gxy:.3f}; implied by the others: {imp}; "
          f"old 'redundant given the others' holds: {oc}  {'PASS' if nc and not oc else 'FAIL'}")
