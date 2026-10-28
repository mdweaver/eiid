# Illustration for the Phase 1 clarification on district (riding) dummies.
# 2x2 for transparency: X = share of "group b" (e.g. previous party voters),
# T = share doing the outcome. 93 ridings x 10 polls. Riding baselines
# (intercepts) are correlated with the riding's mean X; the contextual slopes
# w1, b1 are common to all ridings (model M1 holds; pooled model M0 does not).
# Run: python3 -I illustrate_district_dummies.py
import numpy as np

D, PD, N = 93, 10, 400
W1, B1 = -0.10, 0.20          # common contextual slopes (truth)


def simulate(rng, corr=0.4):
    m = rng.uniform(0.1, 0.9, D)                       # riding mean of X
    w0 = 0.55 + corr * (m - 0.5) + 0.05 * rng.normal(size=D)
    b0 = 0.35 + corr * (m - 0.5) + 0.05 * rng.normal(size=D)
    d = np.repeat(np.arange(D), PD)
    x = np.clip(m[d] + rng.uniform(-0.15, 0.15, D * PD), 0.01, 0.99)
    bw = np.clip(w0[d] + W1 * x + 0.03 * rng.normal(size=x.size), 0, 1)
    bb = np.clip(b0[d] + B1 * x + 0.03 * rng.normal(size=x.size), 0, 1)
    t = x * bb + (1 - x) * bw
    return d, x, t, bb, w0, b0


def pooled_bounds(x, t, n):
    """JKST Prop. 2 (what eiPartialID computes), A = [min x, max x]."""
    w0, c1, d1 = np.linalg.lstsq(np.c_[np.ones_like(x), x, x**2], t, rcond=None)[0]
    lo, hi = -np.inf, np.inf
    for v in (x.min(), x.max()):
        lo = max(lo, -w0 / v, (w0 + c1 + d1 * v - 1) / (1 - v))
        hi = min(hi, (1 - w0) / v, (w0 + c1 + d1 * v) / (1 - v))
    B = lambda w1: np.sum(n * x * (w0 + c1 + d1 * x - w1 * (1 - x))) / np.sum(n * x)
    return B(hi), B(lo)


def riding_bounds(d, x, t, n, coef=None):
    """Model M1: riding-specific intercept & X-slope, common X^2 term.
    Constraint set for w1 is the intersection over ridings of the
    per-riding constraints at that riding's own x-range [l_d, u_d]."""
    if coef is None:
        Z = np.zeros((x.size, 2 * D + 1))
        Z[np.arange(x.size), d] = 1                # riding intercepts w0_d
        Z[np.arange(x.size), D + d] = x            # riding X-slopes  c1_d
        Z[:, -1] = x**2                            # common d1
        beta = np.linalg.lstsq(Z, t, rcond=None)[0]
        w0, c1, d1 = beta[:D], beta[D:2 * D], beta[-1]
    else:
        w0, c1, d1 = coef
    lo, hi = -np.inf, np.inf
    for k in range(D):
        xs = x[d == k]
        for v in (xs.min(), xs.max()):
            lo = max(lo, -w0[k] / v, (w0[k] + c1[k] + d1 * v - 1) / (1 - v))
            hi = min(hi, (1 - w0[k]) / v, (w0[k] + c1[k] + d1 * v) / (1 - v))
    B = lambda w1: np.sum(n * x * (w0[d] + c1[d] + d1 * x - w1 * (1 - x))) / np.sum(n * x)
    return B(hi), B(lo), lo <= hi


def dd_bounds(x, t, n):
    return (np.sum(n * np.maximum(0, t - (1 - x))) / np.sum(n * x),
            np.sum(n * np.minimum(t, x)) / np.sum(n * x))


rng = np.random.default_rng(1)
d, x, t, bb, w0, b0 = simulate(rng)
n = np.full(x.size, N)
truth = np.sum(n * x * bb) / np.sum(n * x)
# population (noise-free) M1 coefficients: c1_d = b0_d - w0_d + w1, d1 = b1 - w1
pop = (w0, b0 - w0 + W1, B1 - W1)
print(f"true B                  : {truth:.4f}")
print("DD (model-free)         : [%.4f, %.4f]" % dd_bounds(x, t, n))
print("pooled M0 (package)     : [%.4f, %.4f]" % pooled_bounds(x, t, n))
print("riding M1, true coefs   : [%.4f, %.4f]  nonempty=%s" % riding_bounds(d, x, t, n, pop))
print("riding M1, estimated    : [%.4f, %.4f]  nonempty=%s" % riding_bounds(d, x, t, n))

# repeat to see coverage
cov = {"pooled": 0, "M1 true coefs": 0, "M1 estimated": 0}
wid = {k: [] for k in cov}
reps = 300
for s in range(reps):
    d, x, t, bb, w0, b0 = simulate(np.random.default_rng(100 + s))
    truth = np.sum(n * x * bb) / np.sum(n * x)
    for k, (L, U, *_) in {"pooled": pooled_bounds(x, t, n),
                          "M1 true coefs": riding_bounds(d, x, t, n, (w0, b0 - w0 + W1, B1 - W1)),
                          "M1 estimated": riding_bounds(d, x, t, n)}.items():
        cov[k] += (L - 1e-9 <= truth <= U + 1e-9)
        wid[k].append(U - L)
for k in cov:
    print(f"{k:15s} covers truth {cov[k] / reps:.2f}, median width {np.median(wid[k]):.4f}")
