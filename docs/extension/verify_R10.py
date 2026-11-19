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

Corrected after the Batch B review (Nov 2026):
  * lengths of a union are now reported as the occupied (merged) length, with the
    hull span and the old unmerged sum beside it (the old script summed the two branch
    lengths even when they overlapped);
  * V9  analytic bound for R10(e) (the referee's repeated three-point design);
  * V10 the containment claimed in 5a §3.5 item 4 fails (the referee's three intervals);
  * V11 an overlapping pair: unmerged sum > Lebesgue length;
  * V12 deliberately false CR (the EI referee's matched-regression experiment).
What these runs do NOT test: the CLR procedure, coverage of the actual-rate target
(residual allowance), the crossing emptiness rule, or a growing-dimension regime.
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


def merge(Is):
    """Merge a list of closed intervals; return the list of disjoint pieces."""
    out = []
    for lo, hi in sorted(Is):
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out


def run(kind, n=150, conc=15, reps=600):
    X, p = design(n)
    parts = linear_parts(X, p)
    theta0, ab = make_truth(X, p, kind)
    f0 = parts[3](theta0)
    assert f0["rows"].min() >= 0, "truth violates [0,1]"
    assert f0["dW"] * f0["dB"] >= -1e-12, "truth violates CR"
    D0 = f0["D"]
    stats = {m: [0, 0, [], [], [], 0] for m in ("JKST", "CR-slab", "CR-union", "CR-naive")}
    for _ in range(reps):
        out = one_rep(X, p, theta0, ab, conc, parts)
        for mth, I in out.items():
            if I is None:
                stats[mth][1] += 1
                continue
            Is = I if isinstance(I, list) else [I]
            stats[mth][0] += any(lo - 1e-12 <= D0 <= hi + 1e-12 for lo, hi in Is)
            pieces = merge(Is)
            stats[mth][2].append(sum(hi - lo for lo, hi in pieces))           # occupied (Lebesgue) length
            stats[mth][3].append(max(hi for _, hi in Is) - min(lo for lo, _ in Is))   # hull span
            stats[mth][4].append(sum(hi - lo for lo, hi in Is))               # old statistic: unmerged sum
            stats[mth][5] += len(Is) == 2 and len(pieces) == 1 and \
                min(Is[0][1], Is[1][1]) - max(Is[0][0], Is[1][0]) > 0     # runs with positive overlap
    print(f"{kind}: true D~ = {D0:.3f}, ER* = {f0['ER']:.3f}, delta*_W = {f0['dW']:.1e}, "
          f"delta*_B = {f0['dB']:.3f}; n = {n}, {reps} reps")
    for mth, (cov, emp, w, span, raw, nov) in stats.items():
        line = (f"   {mth:9s} coverage {cov / reps:.3f}  empty {emp / reps:.3f}  "
                f"median length {np.median(w) if w else float('nan'):.3f}")
        if mth == "CR-union":
            line += (f"  (hull span {np.median(span):.3f}; unmerged sum, old statistic, {np.median(raw):.3f}; "
                     f"positive overlap in {nov}/{reps} runs)")
        print(line)
    return stats


if __name__ == "__main__":
    for kind in ("boundary, ER far from 0", "boundary, ER near 0", "interior"):
        run(kind)


# ---------------------------------------------------------------------------
# Added after the Batch B review (Nov 2026).
# ---------------------------------------------------------------------------
def V9():
    """R10(e), analytic (the partial-identification referee's argument). Repeated design X = 1/4, 1/2, 3/4,
    equal sizes, basis ((1-x)^2, 2x(1-x), x^2). For any coefficient vector q, ER(q) = q11 - q00 and
    c_D(q) = (7/12)(q11 - q00), so se(c_D)/se(ER^) = 7/12 exactly. With the union bound over 8 base rows and
    2 target tails, k_crit <= z_.995, so limsup P(cover) <= Phi(7/12 z_.995). With 4000 critical-value draws,
    add P{Bin(4000, .04) >= 200}."""
    from scipy.stats import norm, binom
    x = np.array([.25, .5, .75]); p = np.ones(3) / 3
    F = np.array([(1 - x) ** 2, 2 * x * (1 - x), x * x]).T
    cov = lambda w, a, b: (w / w.sum()) @ ((a - (w / w.sum()) @ a) * (b - (w / w.sum()) @ b))
    er = np.array([cov(p, x, f) / cov(p, x, x) for f in F.T])
    # c_D(q): D at k = 0 for m1 = q01 (1-x) + q11 x, m0 = q00 (1-x) + q01 x, per basis vector
    cD = np.array([-((1 - x) ** 2).sum() / (1 - x).sum(), 0.0, (x * x).sum() / x.sum()])
    ratio = (cD @ er) / (er @ er)
    pop = norm.cdf(ratio * norm.ppf(.995))
    draws = norm.cdf(ratio * norm.ppf(.996)) + binom.sf(199, 4000, .04)
    # far-edge design of run(): X = .1 + .8 U, U ~ Beta(2,2): ratio .5 + 2 Var X in the design limit
    varX = .64 * (1 / 20)
    far = norm.cdf((.5 + 2 * varX) * norm.ppf(.995))
    return er, cD, ratio, pop, draws, .5 + 2 * varX, far


def V10():
    """5a §3.5 item 4. Raw base interval [0, .1], bracket [.15, .2], target allowance .1. The naive rule that
    intersects the raw base with the bracket and then expands is empty; bracket cap (base expanded) is not.
    So 'the naive interval contains bracket cap JKST CI' fails."""
    base, bracket, allow = (0.0, 0.1), (0.15, 0.2), 0.1
    raw = (max(base[0], bracket[0]), min(base[1], bracket[1]))
    naive = None if raw[0] > raw[1] else (raw[0] - allow, raw[1] + allow)
    jk = (base[0] - allow, base[1] + allow)
    br_jk = (max(jk[0], bracket[0]), min(jk[1], bracket[1]))
    br_jk = None if br_jk[0] > br_jk[1] else br_jk
    contains = br_jk is None or (naive is not None and naive[0] <= br_jk[0] and br_jk[1] <= naive[1])
    return naive, br_jk, contains


def V11():
    """Two overlapping branch intervals: the old statistic (sum of lengths) exceeds the Lebesgue length."""
    Is = [(0.0, 0.3), (0.2, 0.5)]
    return sum(h - l for l, h in Is), sum(h - l for l, h in merge(Is))


def V12(reps=400):
    """The EI referee's matched-regression experiment: m0 = .40 - s x, m1 = .50 + s - s x. The aggregate
    regression is .40 + .10 x for every s; the between sign is positive throughout; for s > 0 both contextual
    slopes are negative, so CR is false. The design is fixed across s at each n. Same generator as the
    referee (seed 1739), so the counts reproduce the report."""
    global RNG
    RNG = np.random.default_rng(1739)
    rows = []
    for n in (150, 600):
        X, p = design(n)
        parts = linear_parts(X, p)
        for s_ in (0, .08, .25):
            ab = (.40, -s_, .50 + s_, -s_)
            theta = np.array([.40, .45, .50, -.05 - s_])
            truth = parts[3](theta)["D"]
            res = {name: [0, 0, []] for name in ("JKST", "CR-slab")}
            for _ in range(reps):
                out = one_rep(X, p, theta, ab, 15, parts)
                for name in res:
                    ci_ = out[name]
                    if ci_ is None:
                        res[name][1] += 1
                    else:
                        res[name][0] += ci_[0] <= truth <= ci_[1]
                        res[name][2].append(ci_[1] - ci_[0])
            rows.append((n, s_, truth, res))
    return rows


if __name__ == "__main__":
    er, cD, ratio, pop, draws, far_ratio, far = V9()
    print(f"V9 [V] R10(e) analytic (three-point design): ER(q) = {er}, c_D(q) = {np.round(cD, 6)}; "
          f"se ratio = {ratio:.6f} (7/12 = {7 / 12:.6f}); limsup coverage <= {pop:.6f}; with 4000 draws <= {draws:.6f}; "
          f"far-edge design ratio {far_ratio:.3f}, bound {far:.3f}  "
          f"{'PASS' if abs(ratio - 7 / 12) < 1e-12 and pop < .95 and draws < .95 else 'FAIL'}")
    naive, br_jk, contains = V10()
    print(f"V10 [F] 5a §3.5 item 4: naive (intersect, then expand) = {naive}; bracket cap JKST CI = {br_jk}; "
          f"old containment claim holds: {contains}  {'PASS' if not contains else 'FAIL'}")
    raw, merged = V11()
    print(f"V11 [V] union length: unmerged sum {raw:.2f} > merged (Lebesgue) {merged:.2f}  "
          f"{'PASS' if raw > merged + 1e-12 else 'FAIL'}")
    for n, s_, truth, res in V12():
        print(f"V12 [I] false CR, n = {n}, s = {s_:.2f}, true gap {truth:.3f}: " +
              "; ".join(f"{name} covers {c}/400, empty {e}, median width {np.median(w):.3f}"
                        for name, (c, e, w) in res.items()))


def V13(reps=600):
    """Round-2 EI referee: which plug-in procedure the 8-10% miss rate describes. Same 600 replications of the
    150-precinct positive-edge design as run() (seed 20261105; the two beta draws of each replication are
    recorded, so the random stream is unchanged). Literal intersection of the fitted bracket hull{gamma ER, ER}
    with the completed JKST interval, versus CR-naive (bracket imposed, then the target's own inflation)."""
    global RNG
    gen = np.random.default_rng(20261105)
    RNG = gen
    x, p = design(150)
    theta, ab = make_truth(x, p, "boundary, ER far from 0")
    parts = linear_parts(x, p)
    truth = parts[3](theta)["D"]
    gamma = parts[2]
    F = np.column_stack(((1 - x) ** 2, 2 * x * (1 - x), x ** 2))
    fit = np.linalg.solve(F.T @ F, F.T)
    er_weights = (p * (x - p @ x)) @ F @ fit / (p @ ((x - p @ x) ** 2))

    class CaptureRates:
        def __init__(self, generator):
            self.generator, self.rates = generator, []
        def beta(self, *args):
            out = self.generator.beta(*args); self.rates.append(out); return out
        def __getattr__(self, name):
            return getattr(self.generator, name)

    RNG = CaptureRates(gen)
    literal = expanded = 0
    for _ in range(reps):
        RNG.rates.clear()
        out = one_rep(x, p, theta, ab, 15, parts)
        b0, b1 = RNG.rates
        er = er_weights @ ((1 - x) * b0 + x * b1)
        lo, hi = sorted([gamma * er, er])
        base = out["JKST"]
        literal += base is not None and max(lo, base[0]) <= truth <= min(hi, base[1])
        ci_ = out["CR-naive"]
        expanded += ci_ is not None and ci_[0] <= truth <= ci_[1]
    RNG = gen
    return literal, expanded


if __name__ == "__main__":
    lit, exp_ = V13()
    print(f"V13 [V] literal intersection (fitted bracket cap completed JKST CI) covers {lit}/600; "
          f"CR-naive (bracket, then target inflation) covers {exp_}/600; the 8-10% miss rate refers to CR-naive only  "
          f"{'PASS' if lit == 304 and exp_ == 550 else 'FAIL'}")
