"""Checks for proof-register item R3 (PHASE5a_sign_module.md §1).

R3: LCE intersected with closed convex restrictions on k is sharp for the
model-implied target; CR is a union whose hull need not be sharp.

  S1  constructive sharpness: for random t in [L, U] (3x2, LCE + delta rows),
      find k with phi(k) = t, build beta_i = M_H(X_i) (2a Prop. 2 step 4) and
      check beta_i in [0,1], x' beta = tau(x), restrictions hold, phi = t.
  S2  2x2: under CR on (delta*_B, delta*_W) exactly one branch is non-empty
      (or both = the point D*_ER = D*_NM), so CR is a single LP.
  S3  3x2 (W1, row r): frequency with which the two CR branches are both
      non-empty and their target intervals are disjoint (a hole).
numpy/scipy only.  Seed 20261105.
"""
import numpy as np
from scipy.optimize import linprog

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # sibling imports also under python3 -I
import verify_package_design as V

RNG = np.random.default_rng(20261105)
V.RNG = RNG
TOL = 1e-9


def assoc_3x2(P, r):
    """delta^{r*}_B and delta^{r*}_W (W1) as affine functions of k: (A, c)."""
    Xs, p = P["Xs"], P["p"]
    mu = p @ Xs

    def m_of(k):
        return Xs @ (P["Q"] + V.Kmat(k)).T          # m[i, s] = (H X_i)_s

    def dB(k):
        m = m_of(k)
        return np.array([sum(p @ (Xs[:, r] * Xs[:, s] * (m[:, r] - m[:, s]))
                             for s in range(3) if s != r)])

    def dW(k):
        m = m_of(k)
        xr = Xs[:, r]
        wr = p * xr / mu[r]
        wn = p * (1 - xr) / (1 - mu[r])
        cov_r = wr @ (m[:, r] * (xr - wr @ xr))
        num = sum(Xs[:, s] * m[:, s] for s in range(3) if s != r)   # (1-mu_r) w^{-r} m^{-r} = p * num
        cov_n = (p * num) @ (xr - wn @ xr) / (1 - mu[r])
        return np.array([mu[r] * cov_r + (1 - mu[r]) * cov_n])

    return V.linear_in_k(dB), V.linear_in_k(dW)


def lp_interval(a, c, A, b):
    out = []
    for s in (1, -1):
        res = linprog(s * a, A_ub=A, b_ub=b, bounds=[(None, None)] * 3, method="highs")
        if res.status == 2:
            return None
        if res.status == 3:
            return (-np.inf, np.inf)
        out.append(s * res.fun + c)
    return tuple(out)


def S1(reps=200):
    worst = 0.0
    checked = 0
    for _ in range(reps):
        P = V.make_3x2()
        A, b = V.rows_3x2(P, ["DD"])
        (aB, cB), (aW, cW) = assoc_3x2(P, 0)
        kt = P["k_true"]
        sB = np.sign(aB[0] @ kt + cB[0]); sW = np.sign(aW[0] @ kt + cW[0])
        A = np.vstack([A, -sB * aB, -sW * aW]); b = np.concatenate([b, sB * cB, sW * cW])
        a, c = V.target_3x2(P, 0)
        I = lp_interval(a[0], c[0], A, b)
        for t in RNG.uniform(I[0], I[1], 5):
            res = linprog(np.zeros(3), A_ub=A, b_ub=b, A_eq=a, b_eq=[t - c[0]],
                          bounds=[(None, None)] * 3, method="highs")
            k = res.x
            H = P["Q"] + V.Kmat(k)
            beta = P["Xs"] @ H.T
            tau = np.einsum("ij,jk,ik->i", P["Xs"], P["Q"], P["Xs"])
            e = max(0.0, -beta.min(), beta.max() - 1,
                    np.abs(np.einsum("ij,ij->i", P["Xs"], beta) - tau).max(),
                    -(sB * (aB[0] @ k + cB[0])), -(sW * (aW[0] @ k + cW[0])),
                    abs(a[0] @ k + c[0] - t))
            worst = max(worst, e); checked += 1
    return checked, worst


def S2(reps=2000):
    """Classify each branch as empty, degenerate (a point) or positive-width, and count the pairs.
    (Relabelled after the Batch B review: the old code counted 'non-empty' only above width 1e-9 and called
    'no positive-width branch' 'both degenerate', which also allowed two empty branches.)"""
    pairs = {}
    for _ in range(reps):
        P = V.make_2x2()
        J = V.jkst_rows(P)
        f = lambda w1: V.funcs_2x2(P, w1)
        aB, bB = V.affine(lambda w: f(w)["dB"])
        aW, bW = V.affine(lambda w: f(w)["dW"])
        Ip = V.interval(J + [(aB, bB), (aW, bW)])
        Im = V.interval(J + [(-aB, -bB), (-aW, -bW)])
        lab = sorted("empty" if I is None else ("point" if I[1] - I[0] <= 1e-9 else "interval")
                     for I in (Ip, Im))
        key = " + ".join(lab)
        pairs[key] = pairs.get(key, 0) + 1
    return pairs


def S3(reps=1000):
    both = holes = 0
    gaps = []
    for support in ("2d", "near-1d"):
        for _ in range(reps // 2):
            P = V.make_3x2(support=support)
            A, b = V.rows_3x2(P, ["DD"])
            (aB, cB), (aW, cW) = assoc_3x2(P, 0)
            a, c = V.target_3x2(P, 0)
            Ip = lp_interval(a[0], c[0], np.vstack([A, -aB, -aW]), np.concatenate([b, cB, cW]))
            Im = lp_interval(a[0], c[0], np.vstack([A, aB, aW]), np.concatenate([b, -cB, -cW]))
            if Ip is None or Im is None:
                continue
            both += 1
            lo, hi = sorted([Ip, Im])
            gap = hi[0] - lo[1]
            if gap > 1e-7:
                holes += 1
                gaps.append(gap / (max(Ip[1], Im[1]) - min(Ip[0], Im[0])))
    return both, holes, (np.median(gaps) if gaps else 0.0), (max(gaps) if gaps else 0.0)


if __name__ == "__main__":
    n, w = S1()
    print(f"S1 constructive sharpness (3x2, LCE + true delta^0_B, delta^0_W signs): "
          f"{n} target values, worst violation {w:.1e}  {'PASS' if w < 1e-7 else 'FAIL'}")
    pairs = S2()
    bad = pairs.get("interval + interval", 0) + pairs.get("empty + empty", 0)
    print("S2 2x2 CR branches (empty / point / positive-width interval): " +
          "; ".join(f"{k} {v}" for k, v in sorted(pairs.items())) +
          f"; frequencies are descriptive (reused design generator)  {'PASS' if bad == 0 else 'FAIL'}")
    both, holes, med, mx = S3()
    print(f"S3 3x2 CR (W1, row 0, target B_0): both branches non-empty {both}; "
          f"disjoint images (hole) {holes}; gap/hull width median {med:.3f}, max {mx:.3f}")


def assoc_alt(P, r, s):
    """Pairwise delta*_B(r,s) and W2 component Delta*_r, affine in k."""
    Xs, p = P["Xs"], P["p"]
    mu = p @ Xs
    m_of = lambda k: Xs @ (P["Q"] + V.Kmat(k)).T

    def dBp(k):
        m = m_of(k)
        return np.array([p @ (Xs[:, r] * Xs[:, s] * (m[:, r] - m[:, s]))])

    def D2(k):
        m = m_of(k)
        tot = 0.0
        for g in range(3):
            w = p * Xs[:, g] / mu[g]
            tot += mu[g] * (w @ (m[:, g] * (Xs[:, r] - w @ Xs[:, r])))
        return np.array([tot])

    return V.linear_in_k(dBp), V.linear_in_k(D2)


def S4(reps=1000):
    """CR with other pairings: (delta^0_B, Delta_0) and (delta_B(0,1), Delta_0)."""
    out = {}
    for name in ("oneVrestB+W2", "pairB+W2"):
        both = holes = 0
        gaps = []
        for i in range(reps):
            P = V.make_3x2(support="2d" if i % 2 == 0 else "near-1d")
            A, b = V.rows_3x2(P, ["DD"])
            (aB1, cB1), _ = assoc_3x2(P, 0)
            (aBp, cBp), (aW, cW) = assoc_alt(P, 0, 1)
            aB, cB = (aB1, cB1) if name == "oneVrestB+W2" else (aBp, cBp)
            a, c = V.target_3x2(P, 0)
            Ip = lp_interval(a[0], c[0], np.vstack([A, -aB, -aW]), np.concatenate([b, cB, cW]))
            Im = lp_interval(a[0], c[0], np.vstack([A, aB, aW]), np.concatenate([b, -cB, -cW]))
            if Ip is None or Im is None:
                continue
            both += 1
            lo, hi = sorted([Ip, Im])
            gap = hi[0] - lo[1]
            if gap > 1e-7:
                holes += 1
                gaps.append(gap / (max(Ip[1], Im[1]) - min(Ip[0], Im[0])))
        out[name] = (both, holes, np.median(gaps) if gaps else 0.0, max(gaps) if gaps else 0.0)
    return out


if __name__ == "__main__":
    for name, (both, holes, med, mx) in S4().items():
        print(f"S4 3x2 CR pairing {name} (target B_0): both branches non-empty {both}/1000; "
              f"hole {holes}; gap/hull width median {med:.3f}, max {mx:.3f}")


def S5(reps=600):
    """Adversarial hole search + connectivity certificate.

    Certificate: if Theta ∩ {g_B = 0, g_W = 0} is non-empty, the two branch
    images share phi(k0) and their union is an interval (R3(c)).
    """
    tried = both = cert = holes = 0
    gaps = []
    for i in range(reps):
        P = V.make_3x2(support=("2d", "near-1d")[i % 2])
        which = (["DD"], ["DD", "O"], ["DD", "O", "L"])[i % 3]
        A, b = V.rows_3x2(P, which)
        (aB1, cB1), _ = assoc_3x2(P, 0)
        (aBp, cBp), (aW, cW) = assoc_alt(P, 0, 1)
        for aB, cB in ((aB1, cB1), (aBp, cBp)):
            Aeq = np.vstack([aB, aW]); beq = -np.concatenate([cB, cW])
            res = linprog(np.zeros(3), A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq,
                          bounds=[(None, None)] * 3, method="highs")
            for r in range(3):
                tried += 1
                a, c = V.target_3x2(P, r)
                Ip = lp_interval(a[0], c[0], np.vstack([A, -aB, -aW]), np.concatenate([b, cB, cW]))
                Im = lp_interval(a[0], c[0], np.vstack([A, aB, aW]), np.concatenate([b, -cB, -cW]))
                if Ip is None or Im is None:
                    continue
                both += 1
                cert += res.status == 0
                lo, hi = sorted([Ip, Im])
                gap = hi[0] - lo[1]
                if gap > 1e-7:
                    holes += 1
                    gaps.append((gap / (max(Ip[1], Im[1]) - min(Ip[0], Im[0])), res.status))
    return tried, both, cert, holes, gaps


if __name__ == "__main__":
    tried, both, cert, holes, gaps = S5()
    print(f"S5 adversarial (3 bases x 2 pairings x 3 targets): {tried} problems; both branches "
          f"non-empty {both}; certificate (Theta meets g_B=g_W=0) holds {cert}; holes {holes}")
    if gaps:
        g = np.array([x for x, _ in gaps])
        print(f"   hole gap/hull median {np.median(g):.3f}, max {g.max():.3f}; "
              f"certificate held in a hole case: {any(s == 0 for _, s in gaps)}")


# ---------------------------------------------------------------------------
# Added after the Batch B review (Nov 2026).
# 2x2 LCE on one column: m1(x) = H11 x + k (1-x), m0(x) = -k x + H22 (1-x), with Q12 = 0 here; in general
# m1 = q11 x + (q01 + k)(1-x), m0 = (q01 - k) x + q00 (1-x), tau = q11 x^2 + 2 q01 x(1-x) + q00 (1-x)^2.
# ---------------------------------------------------------------------------
def _lce_2x2(x, p, q00, q01, q11, verts):
    """Return (k-interval of Theta_O on the vertices, D~(k) as (const, slope), ER, gamma)."""
    rows = []
    for v in verts:                                   # 0 <= m_g(v) <= 1: rows s + a k >= 0
        rows += [(q11 * v + q01 * (1 - v), 1 - v), (1 - q11 * v - q01 * (1 - v), -(1 - v)),
                 (q01 * v + q00 * (1 - v), -v), (1 - q01 * v - q00 * (1 - v), v)]
    lo, hi = -np.inf, np.inf
    for s_, a in rows:
        if abs(a) < 1e-15:
            if s_ < -1e-12:
                return None, None, None, None
            continue
        if a > 0:
            lo = max(lo, -s_ / a)
        else:
            hi = min(hi, -s_ / a)
    kint = (lo, hi) if lo <= hi + 1e-12 else None
    w1, w0 = p * x / (p @ x), p * (1 - x) / (p @ (1 - x))
    m1c, m0c = q11 * x + q01 * (1 - x), q01 * x + q00 * (1 - x)
    D_c = w1 @ m1c - w0 @ m0c
    D_s = w1 @ (1 - x) + w0 @ x                       # dD/dk: m1 gains k(1-x), m0 loses k x
    tau = q11 * x ** 2 + 2 * q01 * x * (1 - x) + q00 * (1 - x) ** 2
    mu = p @ x
    varx = p @ (x - mu) ** 2
    ER = p @ ((x - mu) * tau) / varx
    gamma = varx / (mu * (1 - mu))
    return kint, (D_c, D_s), ER, gamma


def S6():
    """Full law versus regression (R3(a) versus (b)). X = 1/4, 1/2, 3/4 equiprobable, equal sizes, T a fair
    binary total at every composition. Full law: T in {0, 1} forces beta_1 = beta_0 = T in every unit, so both
    row means are .5 everywhere and k = 0. The regression alone (tau = .5) admits m1 = .5 + k(1-x),
    m0 = .5 - k x. With C = {k in [.2, .4]}: full-law identified set empty; regression gap set [1/6, 1/3]."""
    x = np.array([.25, .5, .75]); p = np.ones(3) / 3
    kint, (Dc, Ds), ER, gamma = _lce_2x2(x, p, .5, .5, .5, (x.min(), x.max()))
    lo, hi = max(kint[0], .2), min(kint[1], .4)
    reg = (Dc + Ds * lo, Dc + Ds * hi) if lo <= hi else None
    full_k = {0.0}                                    # the only k compatible with the full law
    full = [k for k in full_k if .2 <= k <= .4]
    return reg, full


def S7():
    """LCE refutes CR although LCE alone fits (5a §3.5 item 2, old 'unrefutable'). X = 0, 1/4, 1 equiprobable,
    A = [0, 1], m0 = 1 - x, m1 = (49/55) x, deterministic T. The Bernstein cross coefficient q01 is 0, so the
    vertex rows force k = 0: LCE fits with a single point. There ER = 0 and D~ = -261/1925, while CR needs
    D~ in hull{gamma ER, ER} = {0}: both CR branches are empty under LCE. In the model-free continuous
    table problem (2c M3) CR is feasible: beta_i1 = beta_i2 = T_i gives delta_B = 0."""
    x = np.array([0., .25, 1.]); p = np.ones(3) / 3
    kint, (Dc, Ds), ER, gamma = _lce_2x2(x, p, 1.0, 0.0, 49 / 55, (0.0, 1.0))
    lce_fits = kint is not None
    Dk = [Dc + Ds * k for k in (kint if lce_fits else ())]
    lo_cr, hi_cr = min(gamma * ER, ER), max(gamma * ER, ER)
    cr_under_lce = lce_fits and max(Dk) >= lo_cr - 1e-12 and min(Dk) <= hi_cr + 1e-12
    # model-free witness: beta_i1 = beta_i0 = T_i, so delta_B = 0 (D = NM), which lies in the CR hull
    tau = (49 / 55) * x ** 2 + (1 - x) ** 2
    mu = p @ x
    dB_witness = p @ (x * (1 - x) * (tau - tau))
    return lce_fits, kint, Dc, ER, cr_under_lce, dB_witness


if __name__ == "__main__":
    reg, full = S6()
    print(f"S6 [V] full law vs regression (C = [.2, .4] on k): regression gap set [{reg[0]:.6f}, {reg[1]:.6f}] "
          f"(1/6, 1/3 = {1 / 6:.6f}, {1 / 3:.6f}); full-law set {'empty' if not full else full}  "
          f"{'PASS' if abs(reg[0] - 1 / 6) < 1e-12 and abs(reg[1] - 1 / 3) < 1e-12 and not full else 'FAIL'}")
    lce_fits, kint, Dc, ER, cr_lce, dBw = S7()
    print(f"S7 [F] LCE excludes CR: LCE fits ({lce_fits}, k in [{kint[0]:.3g}, {kint[1]:.3g}]); ER = {ER:.2e}; "
          f"D~ = {Dc:.6f} (-261/1925 = {-261 / 1925:.6f}); CR feasible under LCE: {cr_lce}; "
          f"model-free witness delta_B = {dBw:.1f} (CR feasible without LCE). Old 'unrefutable' as applied "
          f"in 5a §3.5 fails  {'PASS' if lce_fits and not cr_lce and abs(Dc + 261 / 1925) < 1e-12 and abs(ER) < 1e-12 else 'FAIL'}")
