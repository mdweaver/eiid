# Numerical check for Phase 1, Remark 1.4: the homogeneous-form / vertex
# formulation of the outer set (Theta_O) reproduces Proposition 2 of
# Jiang et al. (2019) as implemented in eiPartialID::generateBounds
# (lambda = 0, useXRangeOffset = TRUE).  Data mimic Example 4 of the paper.
# Run: python3 -I check_2x2_equivalence.py
import numpy as np

rng = np.random.default_rng(7)
p = 1000
n = np.full(p, 150.0)
x = rng.uniform(0, 0.95, p)
lg = lambda z: 1 / (1 + np.exp(-z))
bb = lg(2.197225 - 1.791759 * x + (1 - x) * rng.normal(0, .5, p))
bw = lg(2.197225 + (1 - x) * rng.normal(0, .5, p))
t = bb * x + bw * (1 - x)
l, u = x.min() + 1e-5, x.max() - 1e-5

# --- Package formulas (R/main.R, calcSummaryOutputValues_) ---
w0, c1, d1 = np.linalg.lstsq(np.c_[np.ones(p), x, x ** 2], t, rcond=None)[0]
th = np.array([w0, c1, d1])
r = np.sum(n * x * (1 - x)) / np.sum(n * x)
h = np.array([1.0, 1.0, np.sum(n * x * x) / np.sum(n * x)])
g2 = lambda v: np.array([1 / (1 - v), 1 / (1 - v), 1 / (1 - v) - 1])
gl = [(0, np.array([-1 / l, 0, 0])), (-1 / (1 - l), g2(l)),
      (0, np.array([-1 / u, 0, 0])), (-1 / (1 - u), g2(u))]
gu = [(1 / l, np.array([-1 / l, 0, 0])), (0, g2(l)),
      (1 / u, np.array([-1 / u, 0, 0])), (0, g2(u))]
BL = max(-r * g0 + (h - r * g) @ th for g0, g in gu)
BU = min(-r * g0 + (h - r * g) @ th for g0, g in gl)

# --- Homogeneous form: E[T|X] = X'QX with X = (x, 1-x); H = Q + K ---
q22 = w0
q12 = c1 + 2 * w0          # = 2 * Q_12
q11 = d1 + c1 + w0
lo, hi = -np.inf, np.inf
for v in (l, u):            # vertices of A = [l, u]
    X1, X2 = v, 1 - v
    a = q11 * X1 + q12 / 2 * X2          # row 1 (black): a + k*X2 in [0,1]
    lo, hi = max(lo, -a / X2), min(hi, (1 - a) / X2)
    a = q12 / 2 * X1 + q22 * X2          # row 2 (white): a - k*X1 in [0,1]
    lo, hi = max(lo, (a - 1) / X1), min(hi, a / X1)
mu11 = np.sum(n * x * x) / np.sum(n * x)
mu12 = np.sum(n * x * (1 - x)) / np.sum(n * x)
B = lambda k: q11 * mu11 + (q12 / 2 + k) * mu12

print("package Prop. 2 :", BL, BU)
print("vertex form     :", B(lo), B(hi))
print("true B          :", np.sum(n * x * bb) / np.sum(n * x))
assert abs(BL - B(lo)) < 1e-10 and abs(BU - B(hi)) < 1e-10
