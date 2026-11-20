"""How much does the full conditional law (Theta_X, via averaged unit DD limits) tighten the
regression-only LCE set (Theta_O)?  Population calculation for 2x2 pooled LCE (M0), with the
law of T given X known exactly at replicated composition types (an oracle: no sampling error).

Model: m0(x) = c0 + k x, m1(x) = (c0 + c1 - k) + (c2 + k) x, tau(x) = x m1 + (1-x) m0;
(c0, c1, c2) identified, k free.  Theta_O: 0 <= m_r(x) <= 1 at every support type.
Theta_X adds E[L_i | x] <= m1(x) <= E[U_i | x] (2x2 DD limits for group 1, unit by unit).
Unit rates: beta_ir ~ Beta with mean m_r(x) and concentration phi (heterogeneity across
precincts with the same composition); T_i = x beta_i1 + (1-x) beta_i0 (large precincts).
Target: B1 = population-average model-implied group-1 rate (group-1 weighted).
Reports the width ratio Theta_X / Theta_O for the target.  Run: python3 -I sharp_gain_study.py
"""
import numpy as np
rng = np.random.default_rng(7)

def bounds_k(xs, wts, c, extra=None):
    c0, c1, c2 = c
    lo, hi = -np.inf, np.inf
    def add(a, b, lo, hi, sense):  # constraint a + b*k (sense: >=0 or <= 0)
        # a + b k >= 0
        if abs(b) < 1e-14:
            return lo, hi, (a >= -1e-12)
        r = -a / b
        if b > 0: lo = max(lo, r)
        else: hi = min(hi, r)
        return lo, hi, True
    ok = True
    for j, x in enumerate(xs):
        # m0 = c0 + k x in [0,1]
        for (a, b) in [(c0, x), (1 - c0, -x)]:
            lo, hi, o = add(a, b, lo, hi, 0); ok &= o
        # m1 = (c0+c1+c2 x) + k (x - 1) in [0,1]
        base = c0 + c1 + c2 * x
        for (a, b) in [(base, x - 1), (1 - base, 1 - x)]:
            lo, hi, o = add(a, b, lo, hi, 0); ok &= o
        if extra is not None:
            L, U = extra[j]
            for (a, b) in [(base - L, x - 1), (U - base, 1 - x)]:
                lo, hi, o = add(a, b, lo, hi, 0); ok &= o
    return lo, hi, ok

def study(m0f, m1f, xs, phi, ndraw=20000):
    # identified coefficients from the truth
    X = np.array(xs)
    tau = X * m1f(X) + (1 - X) * m0f(X)
    A = np.vstack([np.ones_like(X), X, X ** 2]).T
    c = np.linalg.lstsq(A, tau, rcond=None)[0]
    extra = []
    for x in xs:
        def draw(m):
            a, b = m * phi, (1 - m) * phi
            return rng.beta(a, b, ndraw)
        b1, b0 = draw(m1f(x)), draw(m0f(x))
        T = x * b1 + (1 - x) * b0
        L = np.maximum(0, (T - (1 - x)) / x); U = np.minimum(1, T / x)
        extra.append((L.mean(), U.mean()))
    w = X / X.sum()  # group-1 weights (equal-size types)
    # target B1(k) = sum w_j m1(x_j) = const + k * sum w_j (x_j - 1)
    slope = np.sum(w * (X - 1))
    loO, hiO, _ = bounds_k(xs, w, c)
    loX, hiX, okX = bounds_k(xs, w, c, extra)
    return abs(slope) * (hiO - loO), abs(slope) * (hiX - loX), okX

designs = {
    "wide hull [.05,.95]": np.linspace(.05, .95, 10),
    "middle hull [.2,.6]": np.linspace(.2, .6, 9),
    "skewed [.02,.4]": np.linspace(.02, .4, 9),
}
truths = {
    "moderate gap": (lambda x: .30 + .20 * x, lambda x: .70 - .10 * x),
    "large gap": (lambda x: .10 + .10 * x, lambda x: .85 + .05 * x),
    "small gap": (lambda x: .45 + .10 * x, lambda x: .55 + .05 * x),
}
print("width of the population-average group-1 rate: Theta_O (regression only) vs Theta_X (full law)")
print(f"{'design':22s} {'truth':14s} {'phi':>5s} {'Theta_O':>8s} {'Theta_X':>8s} {'ratio':>6s}")
ratios = []
for dn, xs in designs.items():
    for tn, (m0f, m1f) in truths.items():
        for phi in [4, 10, 30, 100]:
            wO, wX, ok = study(m0f, m1f, xs, phi)
            ratios.append(wX / wO)
            print(f"{dn:22s} {tn:14s} {phi:5d} {wO:8.3f} {wX:8.3f} {wX/wO:6.2f}")
r = np.array(ratios)
print(f"\nratio Theta_X/Theta_O over {len(r)} cells: median {np.median(r):.2f}, min {r.min():.2f}, "
      f"share of cells with ratio < 0.9: {np.mean(r < .9):.2f}")
print("phi is the Beta concentration of precinct-level group rates around their mean: phi = 4 means\n"
      "very heterogeneous precincts (sd up to .22), phi = 100 nearly homogeneous (sd up to .05).")

# --- High-kink-mass designs (added after consultation with the PACKAGE_DESIGN author):
# bimodal precinct rates (Beta concentration phi < 2 gives U-shaped laws) and small precincts
# (binomial voting with n voters per precinct), so returns at a fixed composition straddle the
# Duncan-Davis kinks.  Same oracle construction; DD limits use the realised share T = votes/n.
def study_small(m0f, m1f, xs, phi, n, ndraw=20000):
    X = np.array(xs)
    tau = X * m1f(X) + (1 - X) * m0f(X)
    c = np.linalg.lstsq(np.vstack([np.ones_like(X), X, X ** 2]).T, tau, rcond=None)[0]
    extra = []
    for x in xs:
        p1 = rng.beta(m1f(x) * phi, (1 - m1f(x)) * phi, ndraw)
        p0 = rng.beta(m0f(x) * phi, (1 - m0f(x)) * phi, ndraw)
        n1 = int(round(x * n)); n0 = n - n1
        T = (rng.binomial(n1, p1) + rng.binomial(n0, p0)) / n
        xe = n1 / n
        L = np.maximum(0, (T - (1 - xe)) / xe); U = np.minimum(1, T / xe)
        extra.append((L.mean(), U.mean()))
    w = X / X.sum(); slope = np.sum(w * (X - 1))
    loO, hiO, _ = bounds_k(xs, w, c); loX, hiX, _ = bounds_k(xs, w, c, extra)
    return abs(slope) * (hiO - loO), abs(slope) * (hiX - loX)

print("\nhigh-kink-mass designs: bimodal precinct rates (phi < 2) and/or small precincts (n voters)")
print(f"{'design':22s} {'truth':14s} {'phi':>5s} {'n':>5s} {'Theta_O':>8s} {'Theta_X':>8s} {'ratio':>6s}")
r2 = []
xs10 = np.linspace(.1, .9, 9)   # compositions that are exact multiples of 1/n for n = 10, 20, 50
for tn, (m0f, m1f) in truths.items():
    for phi in [0.5, 1.0, 4]:
        for n in [10, 20, 50, 1000]:
            wO, wX = study_small(m0f, m1f, xs10, phi, n)
            r2.append(wX / wO)
            print(f"{'grid [.1,.9]':22s} {tn:14s} {phi:5.1f} {n:5d} {wO:8.3f} {wX:8.3f} {wX/wO:6.2f}")
r2 = np.array(r2)
print(f"\nhigh-kink-mass cells: median ratio {np.median(r2):.2f}, min {r2.min():.2f}, share < 0.9: {np.mean(r2 < .9):.2f}")
