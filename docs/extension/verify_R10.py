"""Monte Carlo for proof-register item R10 (PHASE5a_sign_module.md §3):
confidence intervals for the model-implied gap D~ in 2x2 under LCE + CR.

DGP: unit rates beta_ig ~ Beta with conditional mean m_g(X) linear in X (LCE
holds), T_i = (1-X_i) beta_i0 + X_i beta_i1.  Regression of T on the
homogeneous quadratic in (1-X, X); HC3 covariance; conservative CI
(P3 Thm 5(a)) at 95% with the k-rows appended.

Methods compared (all target D~, alpha_E = 0):
  JKST       base rows only
  CR-slab    R10(b): parallel-normal disjunction, each piece inflated, then max
  CR-union   R10(c): union of the two branch CIs, each at the full 95% level
  CR-naive   branch picked by sign(ER^), CR rows plug-in (not inflated)
Designs put the truth on the CR boundary (delta*_W = 0) with ER* far from or
near 0, plus one interior design.  numpy/scipy only.  Seed 20261105.
"""
import numpy as np

RNG = np.random.default_rng(20261105)
ALPHA = 0.05
NDRAW = 4000


def design(n):
    X = np.clip(RNG.beta(2, 2, n) * 0.8 + 0.1, 0.02, 0.98)
    N = RNG.integers(300, 3000, n).astype(float)
    return X, N / N.sum()


def linear_parts(X, p):
    """Every functional below is linear in theta = (q00, q01, q11, k).
    Return coefficient vectors (length 4) for: rows at hull ends, D, dW, dB, ER, and gamma."""
    mu = p @ X
    xL, xU = X.min(), X.max()

    def m(theta, x):
        q00, q01, q11, k = theta
        return q00 * (1 - x) + (q01 + k) * x, (q01 - k) * (1 - x) + q11 * x   # m0, m1

    def funcs(theta):
        m0, m1 = m(theta, X)
        tau = (1 - X) * m0 + X * m1
        q1, q0 = p * X / mu, p * (1 - X) / (1 - mu)
        cov = lambda q, a, b: q @ (a * b) - (q @ a) * (q @ b)
        VarXN = p @ (X - mu) ** 2
        out = dict(
            D=p @ (X * m1) / mu - p @ ((1 - X) * m0) / (1 - mu),
            dW=mu * cov(q1, X, m1) + (1 - mu) * cov(q0, X, m0),
            dB=p @ (X * (1 - X) * (m1 - m0)),
            ER=p @ ((X - mu) * tau) / VarXN,
        )
        rows = []
        for x in (xL, xU):
            a0, a1 = m(theta, x)
            rows += [a0, 1 - a0, a1, 1 - a1]                  # each >= 0
        out["rows"] = np.array(rows)
        return out

    E = np.eye(4)
    z = funcs(np.zeros(4))
    lin = {key: np.array([funcs(E[j])[key] - z[key] for j in range(4)]).T for key in z}
    const = z
    gamma = (p @ (X - mu) ** 2) / (mu * (1 - mu))
    return lin, const, gamma, funcs


def interval_k(rows):
    """rows: list of (s, a) meaning s + a*k >= 0; return feasible [lo, hi] or None."""
    lo, hi = -1e6, 1e6
    for s, a in rows:
        if abs(a) < 1e-12:
            if s < -1e-12:
                return None
            continue
        r = -s / a
        if a > 0:
            lo = max(lo, r)
        else:
            hi = min(hi, r)
    return (lo, hi) if lo <= hi else None


def make_truth(X, p, kind):
    """Linear m_g(x) = a_g + b_g x with delta*_W = 0 (boundary) or interior; CR holds."""
    mu = p @ X
    q1, q0 = p * X / mu, p * (1 - X) / (1 - mu)
    var = lambda q: q @ X ** 2 - (q @ X) ** 2
    if kind == "interior":
        a0, b0, b1, Dt = 0.30, 0.10, 0.15, 0.20
    else:
        a0, b0 = 0.40, 0.20
        b1 = -b0 * (1 - mu) * var(q0) / (mu * var(q1))      # delta*_W = 0
        Dt = 0.15 if kind == "boundary, ER far from 0" else 0.005
    a1 = Dt + a0 + b0 * (q0 @ X) - b1 * (q1 @ X)
    # theta = (q00, q01, q11, k): m0 = q00 (1-x) + (q01+k) x  =>  q00 = a0, q01 + k = a0 + b0
    # m1 = (q01-k)(1-x) + q11 x  =>  q01 - k = a1, q11 = a1 + b1
    q00, q11 = a0, a1 + b1
    q01 = (a0 + b0 + a1) / 2
    k = (a0 + b0 - a1) / 2
    return np.array([q00, q01, q11, k]), (a0, b0, a1, b1)


def one_rep(X, p, theta0, ab, conc, parts):
    lin, const, gamma, funcs = parts
    a0, b0, a1, b1 = ab
    n = len(X)
    m0, m1 = a0 + b0 * X, a1 + b1 * X
    beta0 = RNG.beta(m0 * conc, (1 - m0) * conc)
    beta1 = RNG.beta(m1 * conc, (1 - m1) * conc)
    T = (1 - X) * beta0 + X * beta1
    F = np.column_stack([(1 - X) ** 2, 2 * X * (1 - X), X ** 2])
    XtX_inv = np.linalg.inv(F.T @ F)
    bhat = XtX_inv @ F.T @ T
    e = T - F @ bhat
    h = np.einsum("ij,jk,ik->i", F, XtX_inv, F)
    S = XtX_inv @ (F.T * (e / (1 - h)) ** 2) @ F @ XtX_inv          # HC3

    def piece(key, sign=1.0, idx=None):
        """Return (s_hat, sigma, a_k, grad) for sign * functional, k-coefficient separated."""
        L = lin[key] if idx is None else lin[key][idx]
        c = const[key] if idx is None else const[key][idx]
        g = sign * L[:3]
        return sign * c + g @ bhat, np.sqrt(g @ S @ g), sign * L[3], g

    base = [piece("rows", 1.0, j) for j in range(8)]
    Dp = piece("D")
    ER = piece("ER")
    # CR slab pieces (R10(b)): lower: g_D k + max(s1, s2) >= 0 ; upper: -g_D k + max(s3, s4) >= 0
    gD = Dp[2]
    def lin_beta(gvec, c):
        return c + gvec @ bhat, np.sqrt(gvec @ S @ gvec), gvec
    cD, gDb = const["D"], lin["D"][:3]
    cE, gE = const["ER"], lin["ER"][:3]
    slab = [lin_beta(gDb - gE, cD - cE), lin_beta(gDb - gamma * gE, cD - gamma * cE),
            lin_beta(gE - gDb, cE - cD), lin_beta(gamma * gE - gDb, gamma * cE - cD)]
    plus = [piece("dW", 1.0), piece("dB", 1.0)]
    minus = [piece("dW", -1.0), piece("dB", -1.0)]

    def crit(grads, two_sided_grad):
        """(1-ALPHA) quantile of max(one-sided Z_j, |Z_c|)."""
        G = np.array(list(grads) + [two_sided_grad])
        sd = np.sqrt(np.einsum("ij,jk,ik->i", G, S, G))
        keep = sd > 1e-12
        Z = RNG.multivariate_normal(np.zeros(3), S, NDRAW) @ G[keep].T / sd[keep]
        Zc = np.abs(Z[:, -1]) if keep[-1] else np.zeros(NDRAW)
        Zs = Z[:, :-1] if keep[-1] else Z
        return np.quantile(np.maximum(Zs.max(1) if Zs.shape[1] else -np.inf, Zc), 1 - ALPHA)

    def ci(rows_k, extra_grads):
        grads = [r[3] for r in base] + extra_grads
        kc = crit(grads, Dp[3])
        R = [(s + kc * sd, a) for s, sd, a, _ in base] + rows_k(kc)
        I = interval_k(R)
        if I is None:
            return None
        cD_hat, sdD = Dp[0], Dp[1]
        vals = [cD_hat + gD * I[0], cD_hat + gD * I[1]]
        return (min(vals) - kc * sdD, max(vals) + kc * sdD)

    out = {}
    out["JKST"] = ci(lambda kc: [], [])
    out["CR-slab"] = ci(lambda kc: [(max(slab[0][0] + kc * slab[0][1], slab[1][0] + kc * slab[1][1]), gD),
                                   (max(slab[2][0] + kc * slab[2][1], slab[3][0] + kc * slab[3][1]), -gD)],
                        [s[2] for s in slab])
    Ip = ci(lambda kc: [(s + kc * sd, a) for s, sd, a, _ in plus], [r[3] for r in plus])
    Im = ci(lambda kc: [(s + kc * sd, a) for s, sd, a, _ in minus], [r[3] for r in minus])
    out["CR-union"] = [I for I in (Ip, Im) if I is not None] or None
    br = plus if ER[0] >= 0 else minus
    out["CR-naive"] = ci(lambda kc: [(s, a) for s, sd, a, _ in br], [])
    return out


def run(kind, n=150, conc=15, reps=600):
    X, p = design(n)
    parts = linear_parts(X, p)
    theta0, ab = make_truth(X, p, kind)
    f0 = parts[3](theta0)
    assert f0["rows"].min() >= 0, "truth violates [0,1]"
    assert f0["dW"] * f0["dB"] >= -1e-12, "truth violates CR"
    D0 = f0["D"]
    stats = {m: [0, 0, []] for m in ("JKST", "CR-slab", "CR-union", "CR-naive")}
    for _ in range(reps):
        out = one_rep(X, p, theta0, ab, conc, parts)
        for mth, I in out.items():
            if I is None:
                stats[mth][1] += 1
                continue
            Is = I if isinstance(I, list) else [I]
            stats[mth][0] += any(lo - 1e-12 <= D0 <= hi + 1e-12 for lo, hi in Is)
            stats[mth][2].append(sum(hi - lo for lo, hi in Is))
    print(f"{kind}: true D~ = {D0:.3f}, ER* = {f0['ER']:.3f}, delta*_W = {f0['dW']:.1e}, "
          f"delta*_B = {f0['dB']:.3f}; n = {n}, {reps} reps")
    for mth, (cov, emp, w) in stats.items():
        print(f"   {mth:9s} coverage {cov / reps:.3f}  empty {emp / reps:.3f}  "
              f"median length {np.median(w) if w else float('nan'):.3f}")


if __name__ == "__main__":
    for kind in ("boundary, ER far from 0", "boundary, ER near 0", "interior"):
        run(kind)
