"""
Does nesting (suffrage-yes => Republican for governor) tighten the bound on a SINGLE outcome
(veterans' suffrage-yes rate) once linear contextual effects (M0) are assumed for both contests?
2x2 per contest, population level (identification only).  Rows: veterans (share x), others.
For contest j in {Rep, yes}: m_j,vet(x) = Q11 x + (Q12 + k_j)(1-x), m_j,oth(x) = (Q12 - k_j) x + Q22 (1-x).
Separate: k_yes over its own outer set.  Joint: (k_Rep, k_yes) over both outer sets plus the
nesting rows m_yes,r(v) <= m_Rep,r(v) at the hull ends v in {l, u} (both sides are linear in x).
Truth generated so that nesting holds: m_yes,r = s_r * m_Rep,r evaluated linearly
(m_yes,r(x) = m_Rep,r(x) - d_r(x), d_r linear and >= 0 on [l, u]).
"""
import numpy as np
from scipy.optimize import linprog
rng = np.random.default_rng(8)

def rows(H, l, u):
    """Outer-set rows for one contest in terms of k: returns list of (coef_k, const) with const + coef*k in [0,1]."""
    Q11, Q22, Q12 = H[0, 0], H[1, 1], 0.5 * (H[0, 1] + H[1, 0])
    out = []
    for v in (l, u):
        out.append(((1 - v), Q11 * v + Q12 * (1 - v)))      # vet rate
        out.append((-v, Q12 * v + Q22 * (1 - v)))           # other rate
    return out, (Q11, Q22, Q12)

ratios = []
for _ in range(400):
    l, u = sorted(rng.uniform(0.02, 0.6, 2)); u = max(u, l + 0.15)
    # Rep: H[r, s] = rate of group r in a pure-s unit (r,s in {vet, oth})
    HR = rng.uniform(0.2, 0.9, (2, 2))
    # yes = Rep minus a non-negative linear "Rep-but-no" share
    D = rng.uniform(0.0, 1.0, (2, 2)) * HR
    HY = HR - D
    rR, (_, _, _) = rows(HR, l, u); rY, (Q11, Q22, Q12) = rows(HY, l, u)
    xs = rng.uniform(l, u, 4000); n = rng.integers(200, 2000, 4000)
    mu1, mu11, mu12 = np.mean(n * xs), np.mean(n * xs * xs), np.mean(n * xs * (1 - xs))
    c_y = mu12 / mu1                                       # d B*_yes,vet / d k_yes
    def lp(joint):
        A, b = [], []
        for (a, c) in rY:  # 0 <= c + a k_y <= 1, variables (k_R, k_y)
            A += [[0, -a], [0, a]]; b += [c, 1 - c]
        if joint:
            for (a, c) in rR:
                A += [[-a, 0], [a, 0]]; b += [c, 1 - c]
            for (aR, cR), (aY, cY) in zip(rR, rY):        # cY + aY k_y <= cR + aR k_R
                A += [[-aR, aY]]; b += [cR - cY]
        lo = linprog([0, c_y], A_ub=A, b_ub=b, bounds=[(None, None)] * 2).fun
        hi = -linprog([0, -c_y], A_ub=A, b_ub=b, bounds=[(None, None)] * 2).fun
        return hi - lo
    ratios.append(lp(True) / lp(False))
r = np.array(ratios)
txt = (f"400 random designs: width of veterans' yes-rate bound, joint-with-nesting / separate\n"
       f"median {np.median(r):.2f}; share with any gain (>1%): {np.mean(r < 0.99):.2f}; "
       f"10th percentile {np.quantile(r, .1):.2f}; min {r.min():.2f}\n")
print(txt, end=""); open(__file__.replace(".py", "_output.md"), "w").write(txt)
