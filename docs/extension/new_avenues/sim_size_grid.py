"""Avenue A: how many units, and how small, are needed? Valid (integrated, 1.645-SE) A1 and A1+A2
bounds with ordering, low-heterogeneity design of sim_second_moments.py, Q known; 60 reps per cell."""
import sys; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sim_second_moments as S
S.rng = np.random.default_rng(5)
phi, rho = 60, 0.5
x, n, _ = S.simulate(100000, phi, rho, 300, 0.9)
mu1, mu11, mu12 = np.mean(n*x), np.mean(n*x*x), np.mean(n*x*(1-x))
truth = S.bstar(S.KAPPA0, mu1, mu11, mu12)
def one(p, med):
    xs, ns, ts = S.simulate(p, phi, rho, med, 0.9)
    keep = (xs > 0.06) & (xs < 0.94); xk = xs[keep]
    Xq = np.column_stack([xs**k for k in range(3)]); bq = np.linalg.lstsq(Xq, ts, rcond=None)[0]
    r2 = (ts - Xq @ bq)**2
    Z = np.column_stack([xs**k for k in range(5)] + [xs**k/ns for k in range(4)])
    G = np.linalg.pinv(Z.T @ Z); cf = G @ Z.T @ r2; e = r2 - Z @ cf
    V = G @ (Z.T * e**2) @ Z @ G
    L = np.concatenate([[np.mean(xk**k) for k in range(5)], [np.mean(xk**k) for k in range(4)]])
    tk = np.column_stack([xk**k for k in range(3)]) @ bq
    up = np.mean(tk*(1-tk)) - L @ cf + 1.645*np.sqrt(L @ V @ L)
    LA = np.concatenate([[np.mean(xk**k * xk*(1-xk)/np.minimum(xk,1-xk)**2) for k in range(5)], np.zeros(4)])
    low = np.mean(tk*(1-tk)) - L @ cf - LA @ cf - 1.645*np.sqrt((L+LA) @ V @ (L+LA))
    w = xk*(1-xk); g = S.g(xk); a, b, c = w.mean(), 2*(w*g).mean(), (w*g*g).mean()
    def root(rhs):
        d = b*b - 4*a*(c-rhs); return (-b+np.sqrt(d))/(2*a) if d >= 0 else -np.inf
    lo, hi = S.kappa_interval(xk, l=0.05, u=0.95, ordering=True)
    hi1 = min(hi, root(up)); lo2 = max(lo, root(low)) if low > 0 else lo
    f = lambda k: S.bstar(k, mu1, mu11, mu12)
    return f(hi1)-f(lo), f(hi1)-f(lo2), f(lo) <= truth <= f(hi1)+1e-9, f(lo2) <= truth <= f(hi1)+1e-9
print("O+ordering width (no variance info): 0.363; identification limits: A1 0.262, A1+A2 ~0.00")
print("| median poll size | units | A1 width | A1 cover | A1+A2 width | A1+A2 cover |")
for med in (30, 100, 300, 1000):
    for p in (1000, 5000, 25000):
        R = np.array([one(p, med) for _ in range(60)])
        print(f"| {med} | {p} | {np.median(R[:,0]):.3f} | {R[:,2].mean():.2f} | {np.median(np.where(R[:,1]>=0,R[:,1],np.nan)[~np.isnan(np.where(R[:,1]>=0,R[:,1],np.nan))]) if (R[:,1]>=0).any() else float('nan'):.3f} | {R[:,3].mean():.2f} |", flush=True)
