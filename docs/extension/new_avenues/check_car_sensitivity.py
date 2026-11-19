"""Checks for Section 14 of NEW_AVENUES.md: CAR, CAR-delta and omitted-confounder
(Cinelli-Hazlett / Chernozhukov et al.) sensitivity in ecological inference.

Part 1  Invisible CAR violation: E[B|X] = eta + k*(1-x1, -x1) leaves E[Y|X] exactly
        linear (passes any test of the partially linear form) but biases the CAR
        estimate by k*(1 - E[x1^2]/E[x1]).
Part 2  Omitted-confounder bound for a violation driven by a confounder A = x1 + tau*nu
        that tracks composition. Computes the Chernozhukov et al. (2024) quantities
        (C_gamma, C_alpha, rho, S) by Monte Carlo with quantile binning of A, and compares
        the bound obtained with a "benchmark-like" C_alpha to the true bias, to the
        Duncan-Davis (M0) aggregate bound and to the delta / kappa needed.
Part 3  CAR-delta (sup-norm) in the 2x2 case with composition confined to [a, b]:
        LP half-width of the vertex rate eta_1 versus the closed form
        delta * (2 - a - b) / (b - a), and the half-width for the aggregate target.

Run:  python3 check_car_sensitivity.py  > check_car_sensitivity_output.md
"""
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(20261102)
ETA = np.array([0.6, 0.3])
K = 0.2

print("# CAR sensitivity checks (Section 14)\n")

# ---------------------------------------------------------------- Part 1
print("## Part 1: invisible CAR violation\n")
m = 2_000_000
x1 = rng.beta(2, 2, m)
X = np.column_stack([x1, 1 - x1])
EB1 = ETA[0] + K * (1 - x1)
EB2 = ETA[1] - K * x1
EY = EB1 * x1 + EB2 * (1 - x1)
lin = ETA[0] * x1 + ETA[1] * (1 - x1)
print(f"max |E[Y|X] - eta'X| = {np.max(np.abs(EY - lin)):.2e}  (exactly linear: invisible)")
beta1 = np.sum(x1 * EB1) / np.sum(x1)
print(f"true beta_1 = {beta1:.4f};  CAR/Goodman limit = {ETA[0]:.4f};  bias = {beta1-ETA[0]:+.4f}")
print(f"formula k*(1 - E[x1^2]/E[x1]) = {K*(1-np.mean(x1**2)/np.mean(x1)):+.4f}\n")

# ---------------------------------------------------------------- Part 2
print("## Part 2: omitted-confounder bound when the confounder tracks composition\n")
print("DGP: x1 ~ Beta(2,2); A = x1 + tau*nu, nu ~ N(0,1); B = eta + k*(1-A, -A); Y = B'X.")
print("CAR holds given A; the short (no-A) CEF is exactly eta'X, so the short model is correct")
print("in form and the CAR estimate converges to eta_1. Target beta_1 (equal unit sizes).\n")


def dd_aggregate(x1, y):
    lo = np.maximum(0, (y - (1 - x1)) / x1)
    hi = np.minimum(1, y / x1)
    w = x1 / x1.sum()
    return np.sum(w * lo), np.sum(w * hi)


print("| tau | R2(x1~A) | C_alpha | rho | S | true bias | bound (rho=1, true C) "
      "| bound with C_alpha=0.65 | DD beta_1 interval | delta needed | kappa needed |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
m = 2_000_000
nbins = 400
for tau in [0.02, 0.05, 0.1, 0.2, 0.3]:
    x1 = rng.beta(2, 2, m)
    X = np.column_stack([x1, 1 - x1])
    A = x1 + tau * rng.standard_normal(m)
    B = np.column_stack([ETA[0] + K * (1 - A), ETA[1] - K * A])
    Y = np.sum(B * X, axis=1)
    u = x1 / x1.mean()
    # short model: no covariates. Riesz representer alpha0 = zeta'X, zeta = E[XX']^{-1} e1 E[u]
    Exx = X.T @ X / m
    zeta0 = np.linalg.solve(Exx, np.array([1.0, 0.0]))
    a0 = X @ zeta0
    g0 = X @ ETA                       # short CEF (exact)
    # long model given A: gamma_A = E[Y|X,A] = Y (deterministic); alpha_A = zeta(A)'X
    edges = np.quantile(A, np.linspace(0, 1, nbins + 1))
    idx = np.clip(np.searchsorted(edges, A, side="right") - 1, 0, nbins - 1)
    aA = np.empty(m)
    for b in range(nbins):
        s = idx == b
        Xb = X[s]
        Eb = Xb.T @ Xb / s.sum()
        zb = np.linalg.solve(Eb, np.array([1.0, 0.0]) * u[s].mean())
        aA[s] = Xb @ zb
    gA = Y
    Cg2 = np.mean((gA - g0) ** 2) / np.mean((Y - g0) ** 2)
    Ca2 = (np.mean(aA ** 2) - np.mean(a0 ** 2)) / np.mean(a0 ** 2)
    rho = np.corrcoef(gA - g0, aA - a0)[0, 1]
    S = np.sqrt(np.mean((Y - g0) ** 2) * np.mean(a0 ** 2))
    true_beta = np.sum(x1 * B[:, 0]) / np.sum(x1)
    bias = true_beta - ETA[0]
    bound_true = S * np.sqrt(Cg2) * np.sqrt(max(Ca2, 0))
    bound_bench = S * 1.0 * 0.65
    r2 = np.corrcoef(x1, A)[0, 1] ** 2
    lo, hi = dd_aggregate(x1, np.clip(Y, 0, 1))
    delta_need = K  # sup over x of |E[B|X]-eta| = k*max(1-x1, x1) -> k
    kappa_need = (K ** 2 * np.var(x1)) / np.var(Y)
    print(f"| {tau} | {r2:.3f} | {np.sqrt(max(Ca2,0)):.2f} | {rho:+.2f} | {S:.4f} | {bias:+.4f} "
          f"| {bound_true:.4f} | {bound_bench:.4f} | [{lo:.3f}, {hi:.3f}] | {delta_need:.2f} "
          f"| {kappa_need:.2f} |")
print()
print("Reading: the true bias (0.4k = 0.08) does not depend on tau. As tau -> 0 the confounder")
print("becomes composition itself (a contextual effect): S -> 0 but C_alpha -> infinity. A")
print("C_alpha calibrated on an ordinary covariate (0.65 here, i.e. R2 about 0.3) gives a bound")
print("far below the true bias when tau is small. The Duncan-Davis interval is always valid.\n")

# ---------------------------------------------------------------- Part 3
print("## Part 3: CAR-delta (sup-norm) width and overlap, 2x2, no covariates\n")
print("Composition confined to [a, b]; CAR holds in truth (g(x) = eta'x); |Delta| <= delta.")
print("Half-widths from an LP on a 41-point grid with equal weights.\n")
print("| a | b | delta | eta_1 half-width (LP) | closed form delta(2-a-b)/(b-a) | beta_1 half-width (LP) |")
print("|---|---|---|---|---|---|")


def car_delta_lp(a, b, delta, target):
    xs = np.linspace(a, b, 41)
    n = len(xs)
    g = ETA[0] * xs + ETA[1] * (1 - xs)
    # variables: eta1, eta2, D1[0..n-1], D2[0..n-1]
    nv = 2 + 2 * n
    Aeq = np.zeros((n, nv))
    for i, x in enumerate(xs):
        Aeq[i, 0] = x
        Aeq[i, 1] = 1 - x
        Aeq[i, 2 + i] = x
        Aeq[i, 2 + n + i] = 1 - x
    bounds = [(None, None), (None, None)] + [(-delta, delta)] * (2 * n)
    c = np.zeros(nv)
    if target == "eta1":
        c[0] = 1
    else:  # beta_1 = sum_i w_i x_i (eta1 + D1_i) / sum_i w_i x_i
        w = xs / xs.sum()
        c[0] = 1
        c[2:2 + n] = w
    lo = linprog(c, A_eq=Aeq, b_eq=g, bounds=bounds, method="highs").fun
    hi = -linprog(-c, A_eq=Aeq, b_eq=g, bounds=bounds, method="highs").fun
    centre = ETA[0]
    return (hi - lo) / 2, centre


for a, b in [(0.0, 1.0), (0.1, 0.9), (0.2, 0.8), (0.3, 0.7), (0.4, 0.6)]:
    for delta in [0.02, 0.05]:
        hw_eta, _ = car_delta_lp(a, b, delta, "eta1")
        hw_beta, _ = car_delta_lp(a, b, delta, "beta1")
        cf = delta * (2 - a - b) / (b - a)
        print(f"| {a} | {b} | {delta} | {hw_eta:.4f} | {cf:.4f} | {hw_beta:.4f} |")
print()
print("Reading: under CAR-delta, poor overlap multiplies delta by (2-a-b)/(b-a) for the")
print("vertex rate. This is the partial-identification counterpart of a large Riesz-representer")
print("norm (McCartan & Kuriwaki's positivity diagnostic). The aggregate target beta_1 is")
print("affected about as much: it too needs the two groups' rates separated at each composition.")
