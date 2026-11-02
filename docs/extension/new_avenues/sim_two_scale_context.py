"""
Verification simulation for NEW_AVENUES.md, Avenue B (context measured at a second
scale).  2x2: rows P-voters (share x) / others, column: votes P.  Units (polls) sit in
neighbourhoods with composition z (here: the mean P share of the neighbourhood's polls).

Model:  E[beta_1 | x, z] = a1 + b1 x + g1 z   (loyalty)
        E[beta_0 | x, z] = a0 + b0 x + g0 z   (crossover)
Regression: E[T | x, z] = a0 + (a1-a0+b0) x + (b1-b0) x^2 + g0 z + (g1-g0) x z.
So g1, g0 (row-specific NEIGHBOURHOOD effects) are identified; only b0 (own-poll
"who is behind the pattern") is not.  Restrictions on b0 compared:
  O      : rates in [0,1] at every observed (x, z)        (outer set, conservative)
  RM(M)  : |b_r| <= M |g_r| for r = 0,1  (own-poll effect no larger than M times the
           identified neighbourhood effect; LP rows)
  PROP   : (b1, b0) = lambda (g1, g0)    (context acts the same way at both scales;
           point-identifies b0 when g1 != g0)
Truth: (b1, b0) = lam_true * (g1, g0).
"""
import numpy as np
rng = np.random.default_rng(7)
a1, a0, g1, g0 = 0.60, 0.03, 0.15, 0.25

def simulate(J, per, lam, sd_within=0.12, phi=60):
    b1, b0 = lam * g1, lam * g0
    z = np.clip(rng.beta(2, 2, J), 0.05, 0.95)
    zz = np.repeat(z, per)
    x = np.clip(zz + sd_within * rng.standard_normal(J * per), 0.02, 0.98)
    n = np.maximum(30, np.round(np.exp(np.log(300) + 0.6 * rng.standard_normal(J * per)))).astype(int)
    n1 = np.round(n * x).astype(int); x = n1 / n
    zbar = np.array([x[k * per:(k + 1) * per].mean() for k in range(J)])
    zz = np.repeat(zbar, per)                       # observed neighbourhood composition
    m1 = a1 + b1 * x + g1 * zz; m0 = a0 + b0 * x + g0 * zz
    p1 = rng.beta(m1 * phi, (1 - m1) * phi); p0 = rng.beta(m0 * phi, (1 - m0) * phi)
    t = (rng.binomial(n1, p1) + rng.binomial(n - n1, p0)) / n
    return x, zz, n, t

def fit(x, z, t):
    X = np.column_stack([np.ones_like(x), x, x * x, z, x * z])
    c = np.linalg.lstsq(X, t, rcond=None)[0]
    return dict(a0=c[0], c1=c[1], d=c[2], g0=c[3], g1=c[3] + c[4])

def b0_interval(f, x, z, M=None, prop=False):
    a0_, c1, d, G0, G1 = f["a0"], f["c1"], f["d"], f["g0"], f["g1"]
    lo, hi = -np.inf, np.inf
    # m1 = c1 + a0 - b0 + (d + b0) x + G1 z = base1 + b0 (x - 1);  m0 = a0 + b0 x + G0 z
    base1 = c1 + a0_ + d * x + G1 * z; base0 = a0_ + G0 * z
    for base, s in ((base1, x - 1), (base0, x)):
        # 0 <= base + s b0 <= 1
        pos, neg = s > 1e-12, s < -1e-12
        if pos.any():
            lo = max(lo, np.max((0 - base[pos]) / s[pos])); hi = min(hi, np.min((1 - base[pos]) / s[pos]))
        if neg.any():
            lo = max(lo, np.max((1 - base[neg]) / s[neg])); hi = min(hi, np.min((0 - base[neg]) / s[neg]))
    if M is not None:
        lo = max(lo, -M * abs(G0), -M * abs(G1) - d); hi = min(hi, M * abs(G0), M * abs(G1) - d)
    if prop:
        lam = d / (G1 - G0); b = lam * G0
        lo, hi = max(lo, b), min(hi, b)
    return lo, hi

def bstar_of_b0(f, x, z, n, b0):
    m1 = f["c1"] + f["a0"] + f["d"] * x + f["g1"] * z + b0 * (x - 1)
    return np.sum(n * x * m1) / np.sum(n * x)

def truth_bstar(x, z, n, lam):
    m1 = a1 + lam * g1 * x + g1 * z
    return np.sum(n * x * m1) / np.sum(n * x)

rules = [("O", {}), ("RM(M=2)", dict(M=2)), ("RM(M=1)", dict(M=1)), ("RM(M=0.5)", dict(M=0.5)),
         ("RM(M=0) = no own effect", dict(M=0)), ("PROP", dict(prop=True))]

lines = ["# Avenue B verification: context measured at a second scale (2x2)", "",
         f"a1={a1}, a0={a0}, g1={g1}, g0={g0}; own-poll effects (b1,b0) = lambda (g1,g0).", ""]
for lam in (0.4, 1.5):
    # identification: huge sample, true coefficients
    x, z, n, _ = simulate(4000, 10, lam)
    f_true = dict(a0=a0, c1=a1 - a0 + lam * g0, d=lam * (g1 - g0), g0=g0, g1=g1)
    tb = truth_bstar(x, z, n, lam)
    # DD for loyalty
    t_mean = a0 + f_true["c1"] * x + f_true["d"] * x * x + g0 * z + (g1 - g0) * x * z
    ddl = np.sum(n * np.maximum(0, t_mean - (1 - x))) / np.sum(n * x)
    ddu = np.sum(n * np.minimum(x, t_mean)) / np.sum(n * x)
    lines.append(f"## lambda = {lam} (truth B* = {tb:.3f}; DD at the mean shares [{ddl:.3f}, {ddu:.3f}])")
    lines.append("Identification (population coefficients):")
    lines.append("| restriction | bound on loyalty B* | width | contains truth |")
    lines.append("|---|---|---|---|")
    for lab, kw in rules:
        lo, hi = b0_interval(f_true, x, z, **kw)
        if lo > hi + 1e-12:
            lines.append(f"| {lab} | empty | — | False (refuted) |"); continue
        B = sorted([bstar_of_b0(f_true, x, z, n, lo), bstar_of_b0(f_true, x, z, n, hi)])
        lines.append(f"| {lab} | [{B[0]:.3f}, {B[1]:.3f}] | {B[1]-B[0]:.3f} | {B[0]-1e-9 <= tb <= B[1]+1e-9} |")
    lines.append("")
    lines.append("Estimation (93 neighbourhoods x 10 polls, 300 reps; plug-in):")
    lines.append("| restriction | median width | covers truth | empty |")
    lines.append("|---|---|---|---|")
    acc = {lab: [] for lab, _ in rules}
    for _ in range(300):
        x, z, n, t = simulate(93, 10, lam)
        f = fit(x, z, t); tb_s = truth_bstar(x, z, n, lam)
        for lab, kw in rules:
            lo, hi = b0_interval(f, x, z, **kw)
            if lo > hi + 1e-12:
                acc[lab].append((np.nan, np.nan, tb_s)); continue
            B = sorted([bstar_of_b0(f, x, z, n, lo), bstar_of_b0(f, x, z, n, hi)])
            acc[lab].append((B[0], B[1], tb_s))
    for lab, _ in rules:
        A = np.array(acc[lab]); emp = np.isnan(A[:, 0])
        w = np.median(A[~emp, 1] - A[~emp, 0]) if (~emp).any() else np.nan
        cov = np.mean((A[:, 0] - 1e-9 <= A[:, 2]) & (A[:, 2] <= A[:, 1] + 1e-9))
        lines.append(f"| {lab} | {w:.3f} | {cov:.2f} | {emp.mean():.2f} |")
    errs = []
    for _ in range(300):
        x, z, n, t = simulate(93, 10, lam); f = fit(x, z, t)
        b = f["d"] / (f["g1"] - f["g0"]) * f["g0"]
        errs.append(bstar_of_b0(f, x, z, n, b) - truth_bstar(x, z, n, lam))
    errs = np.array(errs)
    lines.append(f"PROP point estimate (ignoring the O rows): median abs error {np.median(np.abs(errs)):.3f}, "
                 f"5-95% error range [{np.quantile(errs,0.05):.3f}, {np.quantile(errs,0.95):.3f}]")
    lines.append("")
txt = "\n".join(lines); print(txt)
open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
