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
    one_branch = both_point = bad = 0
    for _ in range(reps):
        P = V.make_2x2()
        J = V.jkst_rows(P)
        f = lambda w1: V.funcs_2x2(P, w1)
        aB, bB = V.affine(lambda w: f(w)["dB"])
        aW, bW = V.affine(lambda w: f(w)["dW"])
        Ip = V.interval(J + [(aB, bB), (aW, bW)])
        Im = V.interval(J + [(-aB, -bB), (-aW, -bW)])
        ne = [I is not None and I[1] - I[0] > 1e-9 for I in (Ip, Im)]
        if sum(ne) == 1:
            one_branch += 1
        elif sum(ne) == 0:
            both_point += 1
        else:
            bad += 1
    return one_branch, both_point, bad


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
    o, z, bad = S2()
    print(f"S2 2x2 CR branches: exactly one non-empty {o}, both degenerate {z}, "
          f"both non-empty {bad}  {'PASS' if bad == 0 else 'FAIL'}")
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
