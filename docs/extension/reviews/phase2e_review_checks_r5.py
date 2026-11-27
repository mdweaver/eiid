"""PHASE2e proof review, round 5 (Prop. E36). Independent of verify_phase2e.py: own LP for CI_LCE(eps) in
(m_g0,i, a_g, b_g, c_g, e_g, gamma, delta), own E35 LP (M = 0, Lam = 0, margins only), eps* by (i) bisection and
(ii) one LP with a SINGLE model witness (|D_W(n) - D~| <= t_E), the simplification noted in R5-1.
Run: python3 phase2e_review_checks_r5.py"""
import numpy as np
from scipy.optimize import linprog

res = []
def check(name, ok, info=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{info}]" if info else ""))

x = np.array([.25, .5, .75]); K = 3; aV = 200 * x; aN = 200 - aV; wV = aV / aV.sum()
y0 = np.array([1, 41, 81.]); y1 = np.array([51, 31, 1.])
T0, T1 = y0 / 200, y1 / 200
Z = np.column_stack([np.ones(K), x, x ** 2]); P = np.linalg.pinv(Z)
bh = np.concatenate([P @ T0, P @ (T1 - T0)])

# model block variables: m (6: V then N), a,b (4), c,e (4), gamma (3), delta (3) = 20 ; optional eps (1)
NM = 20
def model_rows(off, nv, eps, h, eps_idx=None):
    Aeq, beq, Aub, bub, bnd = [], [], [], [], {}
    mi = lambda g, i: off + 3 * g + i
    ag, bg = (lambda g: off + 6 + 2 * g), (lambda g: off + 7 + 2 * g)
    cg, eg = (lambda g: off + 10 + 2 * g), (lambda g: off + 11 + 2 * g)
    gm, dl = (lambda j: off + 14 + j), (lambda j: off + 17 + j)
    for g in range(2):
        for i in range(K):
            bnd[mi(g, i)] = (0, 1)
            r = np.zeros(nv); r[mi(g, i)] = 1; r[ag(g)] = -1; r[bg(g)] = -x[i]
            if eps_idx is None:
                Aub += [r, -r]; bub += [eps, eps]
            else:
                r1 = r.copy(); r1[eps_idx] = -1; r2 = -r; r2[eps_idx] = -1; Aub += [r1, r2]; bub += [0, 0]
            p = np.zeros(nv); p[mi(g, i)] = 1; p[cg(g)] = 1; p[eg(g)] = x[i]
            Aub += [p, -p]; bub += [1, 0]
    for j in range(3):
        r = np.zeros(nv); r[gm(j)] = -1
        for i in range(K):
            r[mi(0, i)] += P[j, i] * x[i]; r[mi(1, i)] += P[j, i] * (1 - x[i])
        Aeq.append(r); beq.append(0)
        bnd[gm(j)] = (bh[j] - h, bh[j] + h); bnd[dl(j)] = (bh[3 + j] - h, bh[3 + j] + h)
    for coef in ({dl(0): 1, cg(1): -1}, {dl(1): 1, cg(0): -1, cg(1): 1, eg(1): -1}, {dl(2): 1, eg(0): -1, eg(1): 1}):
        r = np.zeros(nv)
        for k, v in coef.items(): r[k] = v
        Aeq.append(r); beq.append(0)
    tgt = np.zeros(nv); tgt[cg(0)] = wV.sum(); tgt[eg(0)] = wV @ x; tgt[cg(1)] = -wV.sum(); tgt[eg(1)] = -wV @ x
    return Aeq, beq, Aub, bub, bnd, tgt

def lp(c, Aeq, beq, Aub, bub, bnd, nv):
    B = [bnd.get(k, (None, None)) for k in range(nv)]
    r = linprog(c, A_ub=np.array(Aub) if Aub else None, b_ub=np.array(bub) if Aub else None,
                A_eq=np.array(Aeq) if Aeq else None, b_eq=np.array(beq) if Aeq else None, bounds=B, method="highs")
    return r

def ci(eps, h=0.0):
    Aeq, beq, Aub, bub, bnd, t = model_rows(0, NM, eps, h)
    return [s * lp(s * t, Aeq, beq, Aub, bub, bnd, NM).fun for s in (1, -1)]

# E35 block (M = 0, Lam = 0): veteran yes x1_i (base), x2_i (post); theta_i = x1/aV + (dN_i) in [0,1] where
# dN_i = ((y1-x2) - (y0-x1))/aN ; ATT = sum w (x2/aV - theta) = D_W(n).
def e35_block(off, nv):
    Aub, bub, bnd = [], [], {}
    att = np.zeros(nv); c0 = 0.0
    for i in range(K):
        i1, i2 = off + i, off + K + i
        bnd[i1] = (max(0, y0[i] - aN[i]), min(aV[i], y0[i])); bnd[i2] = (max(0, y1[i] - aN[i]), min(aV[i], y1[i]))
        th = np.zeros(nv); th[i1] = 1 / aV[i] + 1 / aN[i]; th[i2] = -1 / aN[i]; th0 = (y1[i] - y0[i]) / aN[i]
        Aub += [th, -th]; bub += [1 - th0, th0]
        att[i2] += wV[i] / aV[i]; att -= wV[i] * th; c0 -= wV[i] * th0
    return Aub, bub, bnd, att, c0

def causal(eps, h=0.0):
    lo, hi = ci(eps, h); nv = 2 * K
    Aub, bub, bnd, att, c0 = e35_block(0, nv)
    Aub += [att, -att]; bub += [hi - c0, -(lo - c0)]
    return [s * lp(s * att, [], [], Aub, bub, bnd, nv).fun + c0 for s in (1, -1)]

def estar_bisect(h=0.0):
    lo, hi = 0.0, 0.5
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if causal(m, h)[0] > 1e-12 else (lo, m)
    return hi

def estar_one_witness(h=0.0, tE=0.0):
    nv = 2 * K + NM + 1; ei = nv - 1
    Aub, bub, bnd, att, c0 = e35_block(0, nv)
    Aeq, beq, A2, b2, bnd2, tgt = model_rows(2 * K, nv, None, h, eps_idx=ei)
    bnd.update(bnd2); bnd[ei] = (0, None)
    Aub += A2; bub += b2
    Aub += [att - tgt, tgt - att]; bub += [tE - c0, tE + c0]     # |D_W(n) - D~| <= tE  (D_W = ATT here)
    Aub += [att]; bub += [-c0]                                    # ATT <= 0
    c = np.zeros(nv); c[ei] = 1
    return lp(c, Aeq, beq, Aub, bub, bnd, nv).fun

R = {e: (ci(e), causal(e)) for e in (0, 0.001, 0.01, 0.5)}
for e, (a, b) in R.items(): print(f"  eps={e}: CI_LCE {np.round(a,4)}  causal {np.round(b,4)}")
check("R5.1 example reproduced: eps=0 [0.0267,0.0800]; .001 [0.0107,0.0960]; .01 [-0.1333,0.2400]; 1/2 CI [-0.4933,0.6867], causal [-0.4933,0.2733]",
      np.allclose(R[0][1], [.0267, .08], atol=6e-5) and np.allclose(R[0.001][1], [.0107, .096], atol=6e-5)
      and np.allclose(R[0.01][1], [-.1333, .24], atol=6e-5) and np.allclose(R[0.5][0], [-.4933, .6867], atol=6e-5)
      and np.allclose(R[0.5][1], [-.4933, .2733], atol=6e-5))
eb, ew = estar_bisect(), estar_one_witness()
ebh, ewh = estar_bisect(0.0002), estar_one_witness(0.0002)
check("R5.2 eps* = 1/600 by bisection and by a one-witness LP; box h = 0.0002 gives 0.00136 by both",
      abs(eb - 1 / 600) < 1e-7 and abs(ew - 1 / 600) < 1e-7 and abs(ebh - ewh) < 1e-7 and abs(ebh - 0.00136) < 6e-6,
      f"{eb:.7f} {ew:.7f}; box {ebh:.7f} {ewh:.7f}")

# one witness == two witnesses for t_E > 0: CI(eps) widened by t_E vs single-witness LP, at a grid of eps
tE = 0.02
def estar_bisect_t(tE):
    lo, hi = 0.0, 0.5
    for _ in range(60):
        m = (lo + hi) / 2
        a, b = ci(m); nv = 2 * K
        Aub, bub, bnd, att, c0 = e35_block(0, nv); Aub += [att, -att]; bub += [b + tE - c0, -(a - tE - c0)]
        l = lp(att, [], [], Aub, bub, bnd, nv).fun + c0
        lo, hi = (m, hi) if l > 1e-12 else (lo, m)
    return hi
et_b, et_w = estar_bisect_t(tE), estar_one_witness(0.0, tE)
check("R5.3 (R5-1) with t_E = 0.02 a single witness |D_W - D~| <= t_E gives the same eps* as CI widened by t_E (two witnesses redundant)",
      abs(et_b - et_w) < 1e-7, f"{et_b:.7f} vs {et_w:.7f}")

# (e) eps*(M) convex nonincreasing: realised route at M>0 needs tau rows; use E35 with tolerance via theta slack
def estar_M(M):
    nv = 2 * K + NM + 1 + K; ei = 2 * K + NM; sl = lambda i: 2 * K + NM + 1 + i
    Aub, bub, bnd, att, c0 = e35_block(0, nv)   # theta here is the M=0 imputation; real theta = that + s_i
    Aub, bub = [], []
    att = np.zeros(nv); c0 = 0.0
    dw = np.zeros(nv); d0 = 0.0
    for i in range(K):
        i1, i2 = i, K + i
        th = np.zeros(nv); th[i1] = 1 / aV[i] + 1 / aN[i]; th[i2] = -1 / aN[i]; th[sl(i)] = 1; th0 = (y1[i] - y0[i]) / aN[i]
        Aub += [th, -th]; bub += [1 - th0, th0]
        tau = M * np.sqrt(1 / aV[i] + 1 / aN[i]); bnd[sl(i)] = (-tau, tau)
        att[i2] += wV[i] / aV[i]; att -= wV[i] * th; c0 -= wV[i] * th0
        dw[i2] += wV[i] * (1 / aV[i] + 1 / aN[i]); dw[i1] -= wV[i] * (1 / aV[i] + 1 / aN[i]); d0 -= wV[i] * th0
    Aeq, beq, A2, b2, bnd2, tgt = model_rows(2 * K, nv, None, 0.0, eps_idx=ei)
    bnd.update(bnd2); bnd[ei] = (0, None); Aub += A2; bub += b2
    Aub += [dw - tgt, tgt - dw]; bub += [-d0, d0]; Aub += [att]; bub += [-c0]
    c = np.zeros(nv); c[ei] = 1
    r = lp(c, Aeq, beq, Aub, bub, bnd, nv)
    return r.fun if r.status == 0 else np.inf
Ms = np.linspace(0, 0.3, 13); es = np.array([estar_M(m) for m in Ms])
d1 = np.diff(es); d2 = np.diff(es, 2)
print("  eps*(M):", np.round(es, 6))
check("R5.4 eps*(M) nonincreasing and convex on a grid (E36(e))", np.all(d1 <= 1e-9) and np.all(d2 >= -1e-9) and abs(es[0] - 1 / 600) < 1e-7)

print(f"\n{sum(res)}/{len(res)} PASS")
