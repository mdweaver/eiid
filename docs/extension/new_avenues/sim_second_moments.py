"""
Verification simulation for NEW_AVENUES.md, Avenue A (second moments of T given X
and unit size).  2x2, rows: P-voters (share x) and others; column: votes P.

Constraints on the unidentified scalar kappa (= K_{1,12}) compared:
  O    : outer set Theta_O (mean rates in [0,1] at the ends of the composition range)
  O+V0 : + assumption-free variance constraint at the smallest observed unit size
          x(1-x)(m1-m0)^2 <= tau(1-tau) - Var(T | x, n_min)
  O+V1 : + binomial-propensity constraint (extrapolate Var(T|x,n)=A(x)+B(x)/n to n=1)
          x(1-x)(m1-m0)^2 <= tau(1-tau) - A(x) - B(x)          [upper bound on row spread]
  O+V1+R : + non-negative correlation of unit shocks across rows
          (m1-m0)^2 >= D(x) - A(x)/min(x,1-x)^2,  D := (tau(1-tau)-A-B)/(x(1-x))
  "+ord" variants add the ordering m1(x) >= m0(x) (loyalty >= crossover), which turns
  every constraint above into a linear (LP) row in kappa.

Identification ("oracle") uses the population A, B, tau; estimation uses OLS of the
squared residuals on [poly4(x), poly3(x)/n] with p units.  Plug-in only (no CI).
Everything is reported for the model-implied loyalty B* = E[N X m1(X)] / E[N X].
"""
import numpy as np
from scipy import stats

rng = np.random.default_rng(20261102)

# truth: loyalty m1(x)=0.75+0.15x, crossover m0(x)=0.05+0.20x
H11, H12, H21, H22 = 0.90, 0.75, 0.25, 0.05
Q11, Q22, Q12 = H11, H22, 0.5 * (H12 + H21)
KAPPA0 = 0.5 * (H12 - H21)

def m1(x, k): return Q11 * x + (Q12 + k) * (1 - x)
def m0(x, k): return (Q12 - k) * x + Q22 * (1 - x)
def g(x): return m1(x, 0.0) - m0(x, 0.0)          # m1-m0 = g(x)+kappa

def draw_pi(x, phi, rho, size):
    """Beta(m*phi,(1-m)*phi) marginals, Gaussian copula with correlation rho."""
    z = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=size)
    u = stats.norm.cdf(z)
    a1, a0 = m1(x, KAPPA0), m0(x, KAPPA0)
    p1 = stats.beta.ppf(u[:, 0], a1 * phi, (1 - a1) * phi)
    p0 = stats.beta.ppf(u[:, 1], a0 * phi, (1 - a0) * phi)
    return p1, p0

def simulate(p, phi, rho, med_n, sdlog):
    x = np.clip(rng.beta(2, 2, p), 0.05, 0.95)
    n = np.maximum(20, np.round(np.exp(np.log(med_n) + sdlog * rng.standard_normal(p)))).astype(int)
    n1 = np.round(n * x).astype(int); x = n1 / n
    z = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=p)
    u = stats.norm.cdf(z)
    a1, a0 = m1(x, KAPPA0), m0(x, KAPPA0)
    p1 = stats.beta.ppf(u[:, 0], a1 * phi, (1 - a1) * phi)
    p0 = stats.beta.ppf(u[:, 1], a0 * phi, (1 - a0) * phi)
    y = rng.binomial(n1, p1) + rng.binomial(n - n1, p0)
    return x, n, y / n

def kappa_interval(grid, D=None, lowD=None, ordering=False, l=None, u=None):
    """Intersect linear/abs constraints on kappa over the grid. Returns (lo, hi)."""
    lo, hi = -np.inf, np.inf
    for v in (l, u):                         # Theta_O: m_r(v) in [0,1]
        # m1 = Q11 v + (Q12+k)(1-v)
        a, b = Q11 * v + Q12 * (1 - v), (1 - v)
        lo = max(lo, (0 - a) / b); hi = min(hi, (1 - a) / b)
        # m0 = (Q12-k) v + Q22 (1-v)
        a, b = Q12 * v + Q22 * (1 - v), -v
        lo = max(lo, (1 - a) / b); hi = min(hi, (0 - a) / b)
    gg = g(grid)
    if ordering:
        lo = max(lo, np.max(-gg))
    if D is not None:
        s = np.sqrt(np.maximum(D, 0))
        hi = min(hi, np.min(s - gg))
        if not ordering:
            lo = max(lo, np.max(-s - gg))
    if lowD is not None and ordering:
        lo = max(lo, np.max(np.sqrt(np.maximum(lowD, 0)) - gg))
    return lo, hi

def bstar(k, mu1, mu11, mu12):
    return (Q11 * mu11 + (Q12 + k) * mu12) / mu1

def design_pop_moments(grid, phi, rho, M=40000):
    """Population A(x), B(x), tau(x) at grid points (Monte Carlo over pi)."""
    A, B, tau = [], [], []
    for x in grid:
        p1, p0 = draw_pi(x, phi, rho, M)
        ti = x * p1 + (1 - x) * p0
        A.append(ti.var()); tau.append(ti.mean())
        B.append(np.mean(x * p1 * (1 - p1) + (1 - x) * p0 * (1 - p0)))
    return np.array(A), np.array(B), np.array(tau)

def fit_moments(x, n, t):
    Xq = np.column_stack([x ** k for k in range(3)])
    beta = np.linalg.lstsq(Xq, t, rcond=None)[0]
    r2 = (t - Xq @ beta) ** 2
    Z = np.column_stack([x ** k for k in range(5)] + [x ** k / n for k in range(4)])
    c = np.linalg.lstsq(Z, r2, rcond=None)[0]
    tau = lambda v: np.column_stack([v ** k for k in range(3)]) @ beta
    A = lambda v: np.column_stack([v ** k for k in range(5)]) @ c[:5]
    Bf = lambda v: np.column_stack([v ** k for k in range(4)]) @ c[5:]
    return tau, A, Bf

def run_design(name, phi, rho, med_n, sdlog, ps=(64, 500, 5000), reps=200):
    out = []
    # population weights from a big sample
    x, n, _ = simulate(200000, phi, rho, med_n, sdlog)
    mu1, mu11, mu12 = np.mean(n * x), np.mean(n * x * x), np.mean(n * x * (1 - x))
    truth = bstar(KAPPA0, mu1, mu11, mu12)
    l, u = 0.05, 0.95
    grid = np.linspace(0.08, 0.92, 29)
    A, B, tau = design_pop_moments(grid, phi, rho, M=200000)
    om = tau * (1 - tau) - A - B
    D = om / (grid * (1 - grid))
    nmin = 20
    D0 = (tau * (1 - tau) - A - B / nmin) / (grid * (1 - grid))
    lowD = D - A / np.minimum(grid, 1 - grid) ** 2
    dd_lo = np.mean(n * np.maximum(0, x - (1 - x) * 0 - 0)) # placeholder (not used)
    rows = {}
    for lab, kw in [("O", {}), ("O+V0", dict(D=D0)), ("O+V1", dict(D=D)),
                    ("O+ord", dict(ordering=True)), ("O+V1+ord", dict(D=D, ordering=True)),
                    ("O+V1+R+ord", dict(D=D, lowD=lowD, ordering=True))]:
        k = kappa_interval(grid, l=l, u=u, **kw)
        rows[lab] = (bstar(k[0], mu1, mu11, mu12), bstar(k[1], mu1, mu11, mu12))
    out.append(f"\n### Design {name}: phi={phi} (Beta dispersion), rho={rho}, unit size median {med_n}, sdlog {sdlog}")
    out.append(f"truth B* = {truth:.4f};  mean Var(pi1-pi0|x) over grid = {np.mean(D - g(grid)-0 if False else D - (g(grid)+KAPPA0)**2):.4f}")
    out.append("Identification (population moments):")
    out.append("| constraint set | bound on loyalty B* | width | contains truth |")
    out.append("|---|---|---|---|")
    w0 = rows["O"][1] - rows["O"][0]
    for lab, (a, b) in rows.items():
        out.append(f"| {lab} | [{a:.3f}, {b:.3f}] | {b-a:.3f} ({(b-a)/w0:.2f} of O) | {a-1e-9 <= truth <= b+1e-9} |")
    # integrated (single-constraint) version: identification
    wgt = grid * (1 - grid)
    def kap_int(Omean, Lmean, xs, gl=l, gu=u, ordering=True):
        lo, hi = kappa_interval(xs, l=gl, u=gu, ordering=ordering)
        w = xs * (1 - xs); gg = g(xs)
        # mean w (g+k)^2 = a k^2 + b k + c
        a, b, c = w.mean(), 2 * (w * gg).mean(), (w * gg * gg).mean()
        def root(rhs):
            disc = b * b - 4 * a * (c - rhs)
            return (-b + np.sqrt(disc)) / (2 * a) if disc >= 0 else -np.inf
        hi = min(hi, root(Omean))
        if Lmean is not None:
            lo = max(lo, root(max(Lmean, 0)) if Lmean > 0 else lo)
        return lo, hi
    xs_pop = x[:20000]; xs_pop = xs_pop[(xs_pop > 0.08) & (xs_pop < 0.92)]
    # population omega and lower bound at the sample x (interpolate grid)
    omp = np.interp(xs_pop, grid, om); Ap = np.interp(xs_pop, grid, A)
    lowp = omp - xs_pop * (1 - xs_pop) * Ap / np.minimum(xs_pop, 1 - xs_pop) ** 2
    kI = kap_int(omp.mean(), None, xs_pop)
    kIR = kap_int(omp.mean(), lowp.mean(), xs_pop)
    out.append(f"| O+V1+ord, integrated (one constraint) | [{bstar(kI[0],mu1,mu11,mu12):.3f}, {bstar(kI[1],mu1,mu11,mu12):.3f}] | {bstar(kI[1],mu1,mu11,mu12)-bstar(kI[0],mu1,mu11,mu12):.3f} | {bstar(kI[0],mu1,mu11,mu12)-1e-9<=truth<=bstar(kI[1],mu1,mu11,mu12)+1e-9} |")
    out.append(f"| O+V1+R+ord, integrated | [{bstar(kIR[0],mu1,mu11,mu12):.3f}, {bstar(kIR[1],mu1,mu11,mu12):.3f}] | {bstar(kIR[1],mu1,mu11,mu12)-bstar(kIR[0],mu1,mu11,mu12):.3f} | {bstar(kIR[0],mu1,mu11,mu12)-1e-9<=truth<=bstar(kIR[1],mu1,mu11,mu12)+1e-9} |")
    # estimation
    out.append("")
    out.append("Estimation, reps=%d.  'pointwise plug-in': grid constraints with estimated A,B (min over grid, no correction)." % reps)
    out.append("'integrated, 95%': one averaged constraint; its right-hand side moved by 1.645 HC0 standard errors in the conservative direction.")
    out.append("| p | pointwise plug-in O+V1+ord: width / cover / empty | integrated 95% O+V1+ord: width / cover / empty | integrated 95% O+V1+R+ord: width / cover / empty |")
    out.append("|---|---|---|---|")
    for p in ps:
        res = {"a": [], "b": [], "c": []}
        for _ in range(reps):
            xs, ns, ts = simulate(p, phi, rho, med_n, sdlog)
            tauf, Af, Bf = fit_moments(xs, ns, ts)
            gl, gu = np.quantile(xs, [0.02, 0.98])
            gr = np.linspace(max(gl, 0.06), min(gu, 0.94), 25)
            om_h = tauf(gr) * (1 - tauf(gr)) - Af(gr) - Bf(gr)
            Dh = om_h / (gr * (1 - gr))
            res["a"].append(kappa_interval(gr, l=l, u=u, D=Dh, ordering=True))
            # integrated: omega_bar = mean[tau(1-tau)] - mean[A+B]; A+B linear in coefs
            keep = (xs > 0.06) & (xs < 0.94); xk = xs[keep]
            Xq = np.column_stack([xs ** k for k in range(3)])
            bq = np.linalg.lstsq(Xq, ts, rcond=None)[0]
            r2 = (ts - Xq @ bq) ** 2
            Z = np.column_stack([xs ** k for k in range(5)] + [xs ** k / ns for k in range(4)])
            ZtZi = np.linalg.pinv(Z.T @ Z); cf = ZtZi @ Z.T @ r2
            e = r2 - Z @ cf
            V = ZtZi @ (Z.T * e ** 2) @ Z @ ZtZi
            L = np.concatenate([[np.mean(xk ** k) for k in range(5)], [np.mean(xk ** k) for k in range(4)]])
            ABbar = L @ cf; se = np.sqrt(L @ V @ L)
            tk = np.column_stack([xk ** k for k in range(3)]) @ bq
            Obar_up = np.mean(tk * (1 - tk)) - ABbar + 1.645 * se
            LA = np.concatenate([[np.mean(xk ** k * xk * (1 - xk) / np.minimum(xk, 1 - xk) ** 2) for k in range(5)], np.zeros(4)])
            low = np.mean(tk * (1 - tk)) - ABbar - LA @ cf
            se_low = np.sqrt((L + LA) @ V @ (L + LA))
            res["b"].append(kap_int(Obar_up, None, xk))
            res["c"].append(kap_int(Obar_up, low - 1.645 * se_low, xk))
        line = f"| {p} |"
        for key in ("a", "b", "c"):
            K = np.array(res[key]); emp = K[:, 0] > K[:, 1]
            lo = bstar(K[:, 0], mu1, mu11, mu12); hi = bstar(K[:, 1], mu1, mu11, mu12)
            w = np.median((hi - lo)[~emp]) if (~emp).any() else np.nan
            cov = np.mean((lo <= truth + 1e-9) & (truth <= hi + 1e-9))
            line += f" {w:.3f} / {cov:.2f} / {emp.mean():.2f} |"
        out.append(line)
    return out

if __name__ == "__main__":
    lines = ["# Avenue A verification: second moments and unit size (2x2)",
             "", "Truth: loyalty m1(x)=0.75+0.15x, crossover m0(x)=0.05+0.20x; kappa0=%.3f." % KAPPA0,
             "In the estimation rows Q is treated as known so that only the variance step is noisy."]
    for args in [("low heterogeneity, NZ-like sizes", 60, 0.5, 300, 0.9),
                 ("high heterogeneity, NZ-like sizes", 8, 0.5, 300, 0.9),
                 ("low heterogeneity, large units", 60, 0.5, 3000, 0.5),
                 ("low heterogeneity, independent shocks", 60, 0.0, 300, 0.9)]:
        lines += run_design(*args)
    txt = "\n".join(lines)
    print(txt)
    open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
