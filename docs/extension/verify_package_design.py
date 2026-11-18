"""Numerical checks for PACKAGE_DESIGN.md.

All quantities are population (mean-level, model-implied) objects under linear
contextual effects (LCE); no sampling noise. Each check prints one block that is
cited in PACKAGE_DESIGN.md as [V1]..[V6].

  V1  Monotone-EI identities under LCE and the JKST -> (delta_W, delta_B) segment
  V2  2x2: global sign / CR restrictions layered on JKST (tighten / redundant / empty)
  V3  2x2: nesting  sign-rows  subset  weighted-local  subset  delta_W>=0 (when x* in hull)
  V4  3x2: sandwich lemma (ordering row normal = - local row normal in k-space)
  V5  3x2: rank of the restriction normals and tightening frequencies
  V6  2x3: joint column restrictions vs the best single-column restriction
  V7  Rx2: rank of pairwise restriction functionals (one pair: R-1; spanning tree: all)

Run: python3 verify_package_design.py > verify_package_design_output.txt
"""
import itertools

import numpy as np
from scipy.optimize import linprog

RNG = np.random.default_rng(20261105)
TOL = 1e-9


# ---------------------------------------------------------------- 2x2 helpers
def make_2x2(n=200):
    """Random 2x2 LCE population. Group 1 share X; group means linear in X."""
    a, b = RNG.uniform(0.5, 4, 2)
    X = np.clip(RNG.beta(a, b, n), 0.02, 0.98)
    N = RNG.integers(200, 2000, n).astype(float)
    p = N / N.sum()
    xL, xU = X.min(), X.max()
    # endpoints of true lines inside [0,1] on [0,1] -> valid everywhere
    m0e, m1e = RNG.uniform(0, 1, 2), RNG.uniform(0, 1, 2)
    w0, w1 = m0e[0], m0e[1] - m0e[0]
    u, v = m1e[0], m1e[1] - m1e[0]
    c1, d1 = u - w0 + w1, v - w1          # tau(x) = w0 + c1 x + d1 x^2
    return dict(X=X, p=p, xL=xL, xU=xU, w0=w0, c1=c1, d1=d1, w1_true=w1)


def funcs_2x2(P, w1):
    X, p, w0, c1, d1 = P["X"], P["p"], P["w0"], P["c1"], P["d1"]
    m0 = w0 + w1 * X
    m1 = w0 + c1 - w1 + (w1 + d1) * X
    tau = X * m1 + (1 - X) * m0
    mu = p @ X
    B1 = p @ (X * m1) / mu
    B0 = p @ ((1 - X) * m0) / (1 - mu)
    D = B1 - B0
    VarXN = p @ (X - mu) ** 2
    D_ER = p @ ((X - mu) * tau) / VarXN
    D_NM = p @ (X * tau) / mu - p @ ((1 - X) * tau) / (1 - mu)
    dB = p @ (X * (1 - X) * (m1 - m0))
    q1, q0 = p * X / mu, p * (1 - X) / (1 - mu)
    cov = lambda q, a, b: q @ (a * b) - (q @ a) * (q @ b)
    dW = mu * cov(q1, X, m1) + (1 - mu) * cov(q0, X, m0)
    return dict(D=D, B1=B1, D_ER=D_ER, D_NM=D_NM, dB=dB, dW=dW,
                VarX=mu * (1 - mu), VarXN=VarXN, mu=mu, q0=q0, q1=q1)


def affine(f):
    """Return (slope, intercept) of an affine scalar function of w1."""
    f0, f1 = f(0.0), f(1.0)
    return f1 - f0, f0


def interval(rows, lo=-50.0, hi=50.0):
    """Feasible w1 set for rows a*w1 + b >= 0. Returns (lo, hi) or None."""
    for a, b in rows:
        if abs(a) < 1e-14:
            if b < -TOL:
                return None
            continue
        r = -b / a
        if a > 0:
            lo = max(lo, r)
        else:
            hi = min(hi, r)
    return (lo, hi) if lo <= hi + TOL else None


def jkst_rows(P):
    """Mean-level DD: m0, m1 in [0,1] at the hull endpoints (linear => enough)."""
    rows = []
    w0, c1, d1 = P["w0"], P["c1"], P["d1"]
    for x in (P["xL"], P["xU"]):
        a0, b0 = x, w0                          # m0 = w0 + w1 x
        a1, b1 = (x - 1), w0 + c1 + d1 * x      # m1 = w0+c1+d1 x + w1 (x-1)
        rows += [(a0, b0), (-a0, 1 - b0), (a1, b1), (-a1, 1 - b1)]
    return rows


def map_interval(P, I, key):
    s, c = affine(lambda w: funcs_2x2(P, w)[key])
    vals = [s * I[0] + c, s * I[1] + c]
    return min(vals), max(vals)


def width(I):
    return 0.0 if I is None else I[1] - I[0]


# ---------------------------------------------------------------- V1
def V1(reps=500):
    err_W = err_B = err_seg = err_back = 0.0
    for _ in range(reps):
        P = make_2x2()
        w = P["w1_true"]
        F = funcs_2x2(P, w)
        # EGGHM Prop. 4 and Prop. 6 (arXiv v1 numbering, see PHASE1b A1)
        err_W = max(err_W, abs((F["D_ER"] - F["D"]) * F["VarXN"] - F["dW"]))
        err_B = max(err_B, abs((F["D_NM"] - F["D"]) * F["VarX"] + F["dB"]))
        # segment: dW/VarXN + dB/VarX = D_ER - D_NM for every w1
        for ww in RNG.uniform(-2, 2, 5):
            G = funcs_2x2(P, ww)
            err_seg = max(err_seg, abs(G["dW"] / G["VarXN"] + G["dB"] / G["VarX"]
                                       - (G["D_ER"] - G["D_NM"])))
        # re-imposing the JKST-implied marginal intervals for (dW, dB) gives back J
        J = interval(jkst_rows(P))
        rows = jkst_rows(P)
        for key in ("dW", "dB"):
            s, c = affine(lambda ww: funcs_2x2(P, ww)[key])
            lo, hi = map_interval(P, J, key)
            rows += [(s, c - lo), (-s, hi - c)]
        J2 = interval(rows)
        DJ, DJ2 = map_interval(P, J, "D"), map_interval(P, J2, "D")
        err_back = max(err_back, abs(DJ[0] - DJ2[0]) + abs(DJ[1] - DJ2[1]))
    print("V1  Monotone-EI identities under LCE (500 random 2x2 populations)")
    print(f"    max |(D_ER-D)Var(X_N) - delta_W|            = {err_W:.2e}")
    print(f"    max |(D_NM-D)Var(X)   + delta_B|            = {err_B:.2e}")
    print(f"    max |dW/Var(X_N)+dB/Var(X)-(D_ER-D_NM)|     = {err_seg:.2e}  (segment)")
    print(f"    max |D-interval change after re-imposing JKST-implied (dW,dB) marginals| = {err_back:.2e}")
    print()


# ---------------------------------------------------------------- V2
def V2(reps=2000):
    stats = {k: dict(tight=0, redund=0, empty=0, cover=0, ratio=[])
             for k in ("dW true sign", "dB true sign", "CR (true)", "dW wrong sign", "dB wrong sign")}
    nJ = 0
    for _ in range(reps):
        P = make_2x2()
        J = interval(jkst_rows(P))
        if J is None:
            continue
        nJ += 1
        F = funcs_2x2(P, P["w1_true"])
        DJ = map_interval(P, J, "D")
        sW = affine(lambda w: funcs_2x2(P, w)["dW"])
        sB = affine(lambda w: funcs_2x2(P, w)["dB"])
        sgnW, sgnB = np.sign(F["dW"]), np.sign(F["dB"])

        def run(name, sets):
            # sets: list of row-lists; union of the resulting intervals
            Is = [interval(jkst_rows(P) + r) for r in sets]
            Ds = [map_interval(P, I, "D") for I in Is if I is not None]
            st = stats[name]
            if not Ds:
                st["empty"] += 1
                return
            D = (min(d[0] for d in Ds), max(d[1] for d in Ds))
            wJ, w = DJ[1] - DJ[0], D[1] - D[0]
            if w < wJ - 1e-7:
                st["tight"] += 1
            else:
                st["redund"] += 1
            st["ratio"].append(w / wJ if wJ > 0 else 1.0)
            st["cover"] += D[0] - 1e-7 <= F["D"] <= D[1] + 1e-7

        sg = lambda s, a: [(s * a[0], s * a[1])]
        run("dW true sign", [sg(sgnW, sW)])
        run("dB true sign", [sg(sgnB, sB)])
        if sgnW == sgnB:                     # CR is imposed only where it is true
            run("CR (true)", [sg(1, sW) + sg(1, sB), sg(-1, sW) + sg(-1, sB)])
        run("dW wrong sign", [sg(-sgnW, sW)])
        run("dB wrong sign", [sg(-sgnB, sB)])
    print(f"V2  2x2: Monotone-EI restrictions layered on JKST ({nJ} populations with non-empty JKST set)")
    print("    restriction        tightens  redundant  empty  truth covered  median width ratio (non-empty)")
    for k, st in stats.items():
        r = np.median(st["ratio"]) if st["ratio"] else float("nan")
        ne = st["tight"] + st["redund"]
        print(f"    {k:18s} {st['tight']:8d} {st['redund']:10d} {st['empty']:6d} {st['cover']:8d}/{ne:<5d} {r:10.3f}")
    print("    (CR rows: only populations where CR is true; imposed as the union of both same-sign LPs,")
    print("     i.e. the common sign is NOT assumed known.)")
    print()


# ---------------------------------------------------------------- V3
def V3(reps=2000):
    sign_in_local = local_in_global = global_in_local = xstar_in = 0
    local_tightens_global = 0
    n = 0
    for _ in range(reps):
        P = make_2x2()
        F = funcs_2x2(P, 0.0)
        X, mu, q0, q1, d1 = P["X"], F["mu"], F["q0"], F["q1"], P["d1"]
        var = lambda q: q @ X ** 2 - (q @ X) ** 2
        a, b = mu * var(q1), (1 - mu) * var(q0)
        xstar = a / (a + b)
        # lower bounds on w1 implied by each restriction (all are w1 >= threshold)
        t_sign = max(0.0, -d1)                       # w1>=0 and w1+d1>=0
        t_local = max(-d1 * P["xL"], -d1 * P["xU"])  # w1 + d1 x >= 0 on hull
        t_glob = -d1 * xstar                         # delta_W >= 0 under LCE (M5(c))
        # check the closed form for delta_W against direct computation
        sW = affine(lambda w: funcs_2x2(P, w)["dW"])
        assert abs(sW[0] - (a + b)) < 1e-9 and abs(sW[1] - a * d1) < 1e-9
        n += 1
        sign_in_local += t_sign >= t_local - 1e-12
        inside = P["xL"] - 1e-12 <= xstar <= P["xU"] + 1e-12
        xstar_in += inside
        if inside:
            local_in_global += t_local >= t_glob - 1e-12
            local_tightens_global += t_local > t_glob + 1e-9
        else:
            global_in_local += t_glob >= t_local - 1e-12
    print(f"V3  2x2 nesting of contextual-direction rows ({n} populations)")
    print(f"    closed form delta_W = (a+b) w1 + a d1 verified in all populations (assert)")
    print(f"    sign rows  =>  weighted-local rows      : {sign_in_local}/{n}")
    print(f"    x* = a/(a+b) inside the composition hull: {xstar_in}/{n}")
    print(f"    ... and then weighted-local => delta_W>=0 : {local_in_global}/{xstar_in}"
          f"  (strictly stronger in {local_tightens_global})")
    print(f"    x* outside hull: delta_W>=0 => local in {global_in_local}/{n - xstar_in} (non-nested otherwise)")
    # adversarial search: 2-4 support points, extreme weights (conjecture C1: x* in hull)
    worst = -np.inf
    for _ in range(50000):
        m = RNG.integers(2, 5)
        X = RNG.uniform(0, 1, m)
        p = RNG.dirichlet(np.ones(m) * RNG.choice([0.1, 1.0]))
        mu = p @ X
        q1, q0 = p * X / mu, p * (1 - X) / (1 - mu)
        var = lambda q: q @ X ** 2 - (q @ X) ** 2
        a, b = mu * var(q1), (1 - mu) * var(q0)
        if a + b < 1e-10:
            continue
        xs = a / (a + b)
        worst = max(worst, X.min() - xs, xs - X.max())
    print(f"    adversarial search (50,000 discrete designs): max hull violation of x* = {worst:.1e}")
    print()


# ---------------------------------------------------------------- 3x2 helpers
PAIRS = [(0, 1), (0, 2), (1, 2)]


def Kmat(k):
    K = np.zeros((3, 3))
    for (r, s), v in zip(PAIRS, k):
        K[r, s], K[s, r] = v, -v
    return K


def linear_in_k(f):
    """f: k -> vector, affine. Return (A, c) with f(k) = A k + c."""
    c = f(np.zeros(3))
    A = np.column_stack([f(np.eye(3)[j]) - c for j in range(3)])
    return A, c


def make_3x2(n=150, support="2d"):
    H = RNG.uniform(0, 1, (3, 3))          # m_r(x) = (H x)_r in [0,1] on simplex
    if support == "2d":
        Xs = RNG.dirichlet([2, 2, 2], n)
    else:                                  # near-1d: points close to a segment
        A, B = RNG.dirichlet([1, 1, 1], 2)
        t = RNG.uniform(0, 1, n)[:, None]
        Xs = t * A + (1 - t) * B + RNG.normal(0, 0.005, (n, 3))
        Xs = np.clip(Xs, 1e-3, None)
        Xs /= Xs.sum(1, keepdims=True)
    N = RNG.integers(200, 2000, n).astype(float)
    p = N / N.sum()
    Q = (H + H.T) / 2
    Kt = (H - H.T) / 2
    k_true = np.array([Kt[r, s] for r, s in PAIRS])
    return dict(Xs=Xs, p=p, Q=Q, k_true=k_true, H=H)


def rows_3x2(P, which, truth_only=True):
    """Return (A_ub, b_ub) for A k <= b."""
    Xs, Q = P["Xs"], P["Q"]
    Hof = lambda k: Q + Kmat(k)
    A, b = [], []

    def add(f, sign=1.0):
        Ak, c = linear_in_k(f)               # f(k) = Ak k + c, require sign*f >= 0
        A.append(-sign * Ak)
        b.append(sign * c)

    if "DD" in which:
        for x in Xs:
            add(lambda k, x=x: Hof(k) @ x)
            add(lambda k, x=x: 1 - Hof(k) @ x)
    Ht = P["H"]
    for r, s in PAIRS:
        if "O" in which:                    # ordering m_r >= m_s (or reverse), all points
            g = Xs @ (Ht[r] - Ht[s])
            if truth_only and not (np.all(g >= 0) or np.all(g <= 0)):
                continue
            sg = 1.0 if np.all(g >= 0) else -1.0
            for x in Xs:
                add(lambda k, x=x: np.array([(Hof(k) @ x)[r] - (Hof(k) @ x)[s]]), sg)
        if "L" in which:                    # weighted-local: (x'H)_r - (x'H)_s sign, all points
            g = Xs @ (Ht[:, r] - Ht[:, s])
            if truth_only and not (np.all(g >= 0) or np.all(g <= 0)):
                continue
            sg = 1.0 if np.all(g >= 0) else -1.0
            for x in Xs:
                add(lambda k, x=x: np.array([(x @ Hof(k))[r] - (x @ Hof(k))[s]]), sg)
    if not A:
        return np.zeros((0, 3)), np.zeros(0)
    return np.vstack(A), np.concatenate(b)


def target_3x2(P, r):
    Xs, p, Q = P["Xs"], P["p"], P["Q"]
    return linear_in_k(lambda k: np.array([p @ (Xs[:, r] * (Xs @ (Q + Kmat(k))[r])) / (p @ Xs[:, r])]))


def bound_3x2(P, which, r):
    A, b = rows_3x2(P, which)
    a, c = target_3x2(P, r)
    out = []
    for s in (1, -1):
        res = linprog(s * a[0], A_ub=A, b_ub=b, bounds=[(None, None)] * 3, method="highs")
        if res.status == 2:
            return None
        if res.status == 3:
            return (-np.inf, np.inf)
        out.append(s * res.fun + c[0])
    return (out[0], out[1])


def V4(reps=200):
    err = 0.0
    for _ in range(reps):
        P = make_3x2(30)
        for x in P["Xs"]:
            for r, s in PAIRS:
                fo = lambda k: np.array([((Kmat(k)) @ x)[r] - ((Kmat(k)) @ x)[s]])
                fl = lambda k: np.array([(x @ Kmat(k))[r] - (x @ Kmat(k))[s]])
                Ao, _ = linear_in_k(fo)
                Al, _ = linear_in_k(fl)
                err = max(err, np.abs(Ao + Al).max())
                # identified half-slope q_rs(x) = (Qx)_r - (Qx)_s ; O + L parts sum to 2 q
                H = P["H"]
                lhs = (H @ x)[r] - (H @ x)[s] + (x @ H)[r] - (x @ H)[s]
                err = max(err, abs(lhs - 2 * ((P["Q"] @ x)[r] - (P["Q"] @ x)[s])))
    print("V4  3x2 sandwich lemma: in k-space the ordering normal for pair (r,s) at x equals minus the")
    print("    weighted-local normal, and (between gap)+(weighted contextual term) = identified directional slope")
    print(f"    max abs error over 200 populations x 30 points x 3 pairs = {err:.2e}")
    print()


def V5(reps=200):
    for support in ("2d", "near-1d"):
        ranks = {"O": [], "L": [], "O+L": []}
        res = {w: dict(tight=0, n=0, cover=0, ratio=[]) for w in ("DD+O", "DD+L", "DD+O+L")}
        unb = 0
        for _ in range(reps):
            P = make_3x2(60, support)
            for w in ("O", "L", "O+L"):
                A, _ = rows_3x2(P, w)
                ranks[w].append(np.linalg.matrix_rank(A, tol=1e-8) if len(A) else 0)
            base = [bound_3x2(P, "DD", r) for r in range(3)]
            if any(b_ is None or not np.isfinite(b_[0]) for b_ in base):
                unb += 1
                continue
            for w in res:
                for r in range(3):
                    B = bound_3x2(P, w, r)
                    a, c = target_3x2(P, r)
                    truth = (a @ P["k_true"] + c)[0]
                    st = res[w]
                    st["n"] += 1
                    if B is None:
                        continue
                    wb, wn = base[r][1] - base[r][0], B[1] - B[0]
                    st["tight"] += wn < wb - 1e-7
                    st["ratio"].append(wn / wb if wb > 0 else 1.0)
                    st["cover"] += B[0] - 1e-7 <= truth <= B[1] + 1e-7
        print(f"V5  3x2, support = {support} ({reps} populations; true restrictions only; dim K = 3)")
        for w, rk in ranks.items():
            cnt = np.bincount(rk, minlength=4)
            print(f"    rank of imposed {w:4s} normals: " + ", ".join(f"rank {i}: {cnt[i]}" for i in range(4)))
        print(f"    DD-only set unbounded/empty in {unb} populations (excluded)")
        print("    rungs     group-rate targets tightened   truth covered   median width ratio vs DD")
        for w, st in res.items():
            r = np.median(st["ratio"]) if st["ratio"] else float("nan")
            print(f"    {w:8s}  {st['tight']:6d}/{st['n']:<6d}                {st['cover']:6d}/{st['n']:<6d}  {r:8.3f}")
        print()


# ---------------------------------------------------------------- V6: 2x3
def V6(reps=1000):
    better_joint = n = 0
    ratios = []
    cover = 0
    for _ in range(reps):
        X = np.clip(RNG.beta(*RNG.uniform(0.5, 4, 2), 200), 0.02, 0.98)
        N = RNG.integers(200, 2000, 200).astype(float)
        p = N / N.sum()
        xL, xU = X.min(), X.max()
        m0_0, m0_1 = RNG.dirichlet([1, 1, 1], 2)      # group-0 probs at x=0, x=1
        m1_0, m1_1 = RNG.dirichlet([1, 1, 1], 2)
        a_, t_true = m0_0, m0_1 - m0_0
        u, v = m1_0, m1_1 - m1_0
        b_, d_ = u - a_ + t_true, v - t_true
        mu = p @ X
        q1, q0 = p * X / mu, p * (1 - X) / (1 - mu)
        var = lambda q: q @ X ** 2 - (q @ X) ** 2
        aw, bw = mu * var(q1), (1 - mu) * var(q0)
        # free k = (t0, t1); t2 = -t0 - t1
        T = np.array([[1, 0], [0, 1], [-1, -1]], float)
        A, b = [], []
        for x in (xL, xU):
            for c in range(3):
                # m0c = a_c + t_c x ; m1c = a_c + b_c + d_c x + t_c (x - 1)
                for coef, const in ((x * T[c], a_[c]), ((x - 1) * T[c], a_[c] + b_[c] + d_[c] * x)):
                    A += [-coef, coef]
                    b += [const, 1 - const]
        A, b = np.array(A), np.array(b)
        # target: group-1 rate of choice 0
        tgt_c = p @ (X * (a_[0] + b_[0] + d_[0] * X)) / mu
        tgt_a = (p @ (X * (X - 1)) / mu) * T[0]
        truth = tgt_a @ np.array([t_true[0], t_true[1]]) + tgt_c

        def rows_dW(c):
            # delta_W,c = aw (t_c + d_c) + bw t_c ; impose its true sign
            val = aw * (t_true[c] + d_[c]) + bw * t_true[c]
            s = 1.0 if val >= 0 else -1.0
            return [-s * (aw + bw) * T[c]], [s * aw * d_[c]]

        def bound(extra):
            AA, bb = A, b
            for c in extra:
                ra, rb = rows_dW(c)
                AA, bb = np.vstack([AA, ra]), np.concatenate([bb, rb])
            out = []
            for s in (1, -1):
                r = linprog(s * tgt_a, A_ub=AA, b_ub=bb, bounds=[(None, None)] * 2, method="highs")
                if r.status != 0:
                    return None
                out.append(s * r.fun + tgt_c)
            return out
        base = bound([])
        if base is None:
            continue
        singles = [bound([c]) for c in range(3)]
        joint = bound([0, 1, 2])
        n += 1
        ws = min(s_[1] - s_[0] for s_ in singles)
        wj = joint[1] - joint[0]
        better_joint += wj < ws - 1e-7
        ratios.append(wj / (base[1] - base[0]))
        cover += joint[0] - 1e-7 <= truth <= joint[1] + 1e-7
    print(f"V6  2x3: true delta_W signs on all three columns jointly vs best single column ({n} populations)")
    print(f"    joint strictly narrower than the best single-column restriction: {better_joint}/{n}")
    print(f"    truth covered by joint bound: {cover}/{n};  median joint width / JKST width = {np.median(ratios):.3f}")
    print()


# ---------------------------------------------------------------- V7: rank lemma, R groups
def V7():
    print("V7  Rank of pairwise ordering / weighted-local functionals g_rs(x;K) = (Kx)_r - (Kx)_s")
    print("    (binary outcome, full-dimensional composition support, 200 Dirichlet points)")
    for R in (3, 4, 5):
        pairs = list(itertools.combinations(range(R), 2))
        d = len(pairs)
        Xs = RNG.dirichlet(np.ones(R), 200)

        def Km(k):
            K = np.zeros((R, R))
            for (r, s), v in zip(pairs, k):
                K[r, s], K[s, r] = v, -v
            return K

        def rank(prs):
            rows = [[(Km(e) @ x)[r] - (Km(e) @ x)[s] for e in np.eye(d)] for x in Xs for r, s in prs]
            return np.linalg.matrix_rank(np.array(rows), 1e-9)
        print(f"    R={R}: dim K={d}; one pair: {rank([pairs[0]])}; star tree: "
              f"{rank([(0, s) for s in range(1, R)])}; chain tree: {rank([(i, i + 1) for i in range(R - 1)])}; "
              f"all pairs: {rank(pairs)}")
    print()


if __name__ == "__main__":
    V1()
    V2()
    V3()
    V4()
    V5()
    V6()
    V7()
