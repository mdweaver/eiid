"""Checks of referee report B's claims that the response plan relies on.
Run from anywhere: python3 -I check_B_pi_claims.py   (NumPy and SciPy)."""
from itertools import product
import numpy as np
from scipy.stats import norm

rng = np.random.default_rng(20261118)

def cov(w, a, b):
    w = w / w.sum()
    return (w * (a - w @ a)) @ (b - w @ b)

# ---------------------------------------------------------------------------
# I1. Integer N1 in 2x2 with margins only: the integer D-set is the continuous
#     DD interval intersected with a lattice of spacing N/(N1*N0) = 1/(N mu (1-mu)).
#     So CR / sign sets differ from the continuous ones only within one lattice step.
# ---------------------------------------------------------------------------
def integer_check(N, r, y):
    N, r, y = map(np.asarray, (N, r, y))
    lo = np.maximum(0, y - (N - r)); hi = np.minimum(r, y)
    R1, R0, Y = r.sum(), (N - r).sum(), y.sum()
    Dof = lambda S: S / R1 - (Y - S) / R0
    x, t, p = r / N, y / N, N / N.sum()
    mu = p @ x; var = cov(p, x, x); ER = cov(p, x, t) / var; gam = var / (mu * (1 - mu))
    NM = gam * ER
    # brute force over unit allocations (small cases only)
    if np.prod(hi - lo + 1) <= 20000:
        Ss = {sum(z) for z in product(*[range(a, b + 1) for a, b in zip(lo, hi)])}
        assert Ss == set(range(lo.sum(), hi.sum() + 1))      # every integer S is attained
    Dmin, Dmax = Dof(lo.sum()), Dof(hi.sum())
    a, b = sorted((NM, ER))
    cont_cr = (max(a, Dmin), min(b, Dmax))                   # continuous CR interval
    grid = [Dof(S) for S in range(lo.sum(), hi.sum() + 1) if a - 1e-12 <= Dof(S) <= b + 1e-12]
    step = 1 / R1 + 1 / R0
    return cont_cr, grid, step, (1 - gam) * abs(ER)

# referee's example
cc, grid, step, w = integer_check([4, 4, 4], [1, 2, 3], [0, 1, 0])
assert grid == [] and abs(step - 1/3) < 1e-12 and np.allclose(cc, 0)
print("I1 referee example: continuous CR {0}, integer CR empty, lattice step", step)

# random small populations: integer CR empty iff the continuous CR width < one step
# (up to the position of the lattice), and endpoints differ by < one step
flips = maxdiff = 0; n = 0
for _ in range(3000):
    p_ = rng.integers(3, 7); N = rng.integers(2, 9, p_)
    r = np.array([rng.integers(0, k + 1) for k in N]); y = np.array([rng.integers(0, k + 1) for k in N])
    if len(set(r / N)) < 2: continue
    cc, grid, step, w = integer_check(N, r, y); n += 1
    if not grid: flips += 1; assert w < step + 1e-12
    else: maxdiff = max(maxdiff, (grid[0] - cc[0]) / step, (cc[1] - grid[-1]) / step)
assert maxdiff < 1 + 1e-9
print(f"I1 {n} small populations: integer CR empty in {flips} (all with width < one step); "
      f"endpoint gap <= {maxdiff:.3f} steps")
# realistic size: Wisconsin-like 370 units of ~2000 people
N = rng.integers(500, 4000, 370); x = rng.beta(2, 5, 370)
print("I1 lattice step at 370 units of 500-4000 people:",
      1 / (N.sum() * (x @ N / N.sum()) * (1 - x @ N / N.sum())))

# ---------------------------------------------------------------------------
# I2. District rungs (M5(c), P-N). Referee's design: common positive slopes in both
#     rows, district intercepts. Pooled delta*_W < 0; within each district >= 0.
#     The 2x2 local rows x m1' + (1-x) m0' >= 0 hold everywhere, so P-N's chain
#     "signed slopes => weighted-local => delta_W >= 0" fails for pooled M1.
# ---------------------------------------------------------------------------
x = np.array([.1, .2, .3, .7, .8, .9]); d = (x > .5).astype(int)
m1 = np.where(d == 0, .8, .2) + .1 * x; m0 = m1.copy()
p = np.ones(6) / 6; mu = p @ x
def dW(xs, a1, a0, ps):
    m = ps @ xs
    return m * cov(ps * xs, xs, a1) + (1 - m) * cov(ps * (1 - xs), xs, a0)
pooled = dW(x, m1, m0, p)
within = [dW(x[d == k], m1[d == k], m0[d == k], p[d == k] / p[d == k].sum()) for k in (0, 1)]
local = x * .1 + (1 - x) * .1
assert pooled < 0 and min(within) > 0 and np.all(local > 0)
print(f"I2 pooled delta*_W = {pooled:.5f} < 0; per-district {within[0]:.5f}, {within[1]:.5f} > 0; "
      "local rows all positive")
# The W1 slab identity still holds with the district-fitted tau (what the code uses):
tau = x * m1 + (1 - x) * m0
ER = cov(p, x, tau) / cov(p, x, x)
D = (p @ (x * m1)) / mu - (p @ ((1 - x) * m0)) / (1 - mu)
assert np.isclose(pooled, cov(p, x, x) * (ER - D))
print("I2 W1 identity delta*_W = Var(X)(ER* - D*) holds with district tau: the CR slab rows stay valid under M1")

# ---------------------------------------------------------------------------
# I3. Selection after LOO-CV (2c §6): conditional coverage 1 - .05 / P(|Z| > sqrt 2)
# ---------------------------------------------------------------------------
print("I3 conditional coverage after CV selection:", 1 - .05 / (2 * norm.sf(np.sqrt(2))))

# ---------------------------------------------------------------------------
# I4. Singular-covariance projection (P3 §5): error (Z/sqrt n, -1/n); se of the second
#     coordinate is 0 in the limit, so a zero-width upper bound misses surely.
# ---------------------------------------------------------------------------
n = 10_000; Z = rng.standard_normal(5000)
err2 = -np.ones_like(Z) / n
assert np.all(0 + err2 < 0)          # estimate + 0 * k < truth always
print("I4 null-direction bound with zero s.e. misses in 100% of draws")
print("ALL CHECKS PASS")
