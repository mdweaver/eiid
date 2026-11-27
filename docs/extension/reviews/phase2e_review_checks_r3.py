"""
Round-3 checks for the proof review of PHASE2e (independent-review revisions, sections 5.8-5.10).
Reviewer's own code; nothing imported from verify_phase2e.py or the independent reviewer's scripts.
Only new claims are checked here; claims already covered by passing checks are not rebuilt.
Run: python3 phase2e_review_checks_r3.py > phase2e_review_checks_r3_output.txt
"""
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(20261126)
RES = []


def check(name, ok, detail=""):
    RES.append(bool(ok))
    print(("PASS  " if ok else "FAIL  ") + name + (("\n      " + detail) if detail else ""), flush=True)


# ------------------------------------------------------------------------------------------
# 1. E30(c): three-precinct LCE combinations and epsilon* = 1/600, by an independent LP.
#    Variables per precinct i: mV0, mN0, dV, dN (12), plus a_g, b_g for curvature (4), plus eps (1).
# ------------------------------------------------------------------------------------------
x = np.array([.25, .5, .75]); T0 = np.array([.005, .205, .405]); T1 = np.array([.255, .155, .005])
wV = x / x.sum()
iv = lambda k, i: 4 * i + k      # k: 0 mV0, 1 mN0, 2 dV, 3 dN
NV = 12 + 4 + 1


def e30(change=False, base=False, post=False, eps=None, target_side=None):
    Aeq, beq, Aub, bub = [], [], [], []
    def row(d):
        v = np.zeros(NV)
        for k, c in d.items():
            v[k] += c
        return v
    for i in range(3):
        Aeq.append(row({iv(0, i): x[i], iv(1, i): 1 - x[i]})); beq.append(T0[i])
        Aeq.append(row({iv(0, i): x[i], iv(1, i): 1 - x[i], iv(2, i): x[i], iv(3, i): 1 - x[i]})); beq.append(T1[i])
        for g in (0, 1):   # factual post rates in [0,1]
            Aub.append(row({iv(g, i): 1, iv(2 + g, i): 1})); bub.append(1)
            Aub.append(row({iv(g, i): -1, iv(2 + g, i): -1})); bub.append(0)
        Aub.append(row({iv(0, i): 1, iv(3, i): 1})); bub.append(1)      # WPT rate limit
        Aub.append(row({iv(0, i): -1, iv(3, i): -1})); bub.append(0)
    curv = lambda k: row({iv(k, 0): 1, iv(k, 1): -2, iv(k, 2): 1})   # equally spaced x: affine iff second difference 0
    for g in (0, 1):
        if change:
            Aeq.append(curv(2 + g)); beq.append(0)
        if base:
            Aeq.append(curv(g)); beq.append(0)
        if post:
            Aeq.append(curv(g) + curv(2 + g)); beq.append(0)
    if eps is not None:   # |m_g0(x_i) - a_g - b_g x_i| <= eps
        for g in (0, 1):
            for i in range(3):
                for s in (1, -1):
                    Aub.append(row({iv(g, i): s, 12 + 2 * g: -s, 13 + 2 * g: -s * x[i], 16: -1})); bub.append(0)
    tgt = row({**{iv(2, i): wV[i] for i in range(3)}})
    for i in range(3):
        tgt[iv(3, i)] -= wV[i]
    bnds = [(0, 1), (0, 1), (-1, 1), (-1, 1)] * 3 + [(None, None)] * 4 + [(0, None) if eps == "free" else (eps or 0, eps or 0) if eps is not None else (0, 0)]
    if eps == "free":
        Aub.append(tgt); bub.append(0)        # target <= 0
        c = np.zeros(NV); c[16] = 1
        r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bnds, method="highs")
        return r.fun
    out = []
    for s in (1, -1):
        r = linprog(s * tgt, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bnds, method="highs")
        out.append(s * r.fun)
    return tuple(out)


res = {"none": e30(), "change": e30(change=True), "base+change": e30(change=True, base=True),
       "base+post": e30(base=True, post=True), "post+change": e30(change=True, post=True),
       "all three": e30(True, True, True)}
exp = {"none": (-29 / 50, 41 / 150), "change": (-37 / 75, 11 / 150)}
for k in ("base+change", "base+post", "post+change", "all three"):
    exp[k] = (2 / 75, 11 / 150)
ok = all(np.allclose(res[k], exp[k], atol=1e-7) for k in exp)
eps_star = e30(change=True, eps="free")
lo_below = e30(change=True, eps=1 / 600 - 1e-4)[0]
lo_half = e30(change=True, eps=0.5)
check("E30(c) [V] six LCE combinations by an independent LP; eps* = 1/600 by one LP; eps = 1/2 = change-only",
      ok and abs(eps_star - 1 / 600) < 1e-9 and lo_below > 0 and np.allclose(lo_half, res["change"]),
      "; ".join(f"{k} [{v[0]*100:.2f},{v[1]*100:.2f}]" for k, v in res.items()) + f"; eps* = {eps_star:.6f}")


# ------------------------------------------------------------------------------------------
# 2. E31 (empirical Bernstein, Maurer-Pontil Thm 11 rescaled, two-sided) and E29(b) (Hoeffding):
#    coverage of E[D_W] by D_W(n0) +- r under precinct-wide common shocks, unequal weights,
#    non-identical precinct laws. (The projection over F only enlarges the set.)
# ------------------------------------------------------------------------------------------
def sim_cov(K=40, reps=20000, alpha=.05):
    w = rng.dirichlet(np.ones(K) * 2)                     # unequal veteran weights
    pV0 = rng.uniform(.1, .9, K); pN0 = rng.uniform(.1, .9, K)
    # each precinct: with prob q_i a common shock flips everyone (all-yes or all-no), else rates stay near base
    q = rng.uniform(0, .6, K)
    mu_i = np.zeros(K)
    Z = np.zeros((reps, K))
    for i in range(K):
        shock = rng.random(reps) < q[i]
        dir_ = rng.random(reps) < .5
        dV = np.where(shock, np.where(dir_, 1 - pV0[i], -pV0[i]), 0.0)          # all veterans to yes / no
        dN = np.where(shock, np.where(rng.random(reps) < .3, 1 - pN0[i], -pN0[i]), 0.0)
        Z[:, i] = w[i] * (dV - dN)
        mu_i[i] = w[i] * (q[i] * (.5 * (1 - pV0[i]) - .5 * pV0[i]) - q[i] * (.3 * (1 - pN0[i]) - .7 * pN0[i]))
    DW, mu = Z.sum(1), mu_i.sum()
    R = 4 * K * w.max()
    Y = K * Z
    s2 = Y.var(1, ddof=1)
    lg = np.log(4 / alpha)
    r_eb = np.sqrt(2 * s2 * lg / K) + 7 * R * lg / (3 * (K - 1))
    r_h = np.sqrt(.5 * np.log(2 / alpha) * np.sum((4 * w) ** 2))
    return np.mean(np.abs(DW - mu) <= r_eb), np.mean(np.abs(DW - mu) <= r_h), r_eb.mean(), r_h, 1 / np.sum(w ** 2)


ceb, ch, reb, rh, neff = sim_cov()
check("E31/E29(b) [V] coverage >= .95 under precinct-wide shocks, unequal weights, non-identical laws",
      ceb >= .95 and ch >= .95, f"empirical Bernstein {ceb:.4f} (mean radius {reb:.3f}); Hoeffding {ch:.4f} (radius {rh:.3f}); n_eff {neff:.1f}")

# E31(c) arithmetic
check("E31(c) [V] additive term 7*4*log(80)/(3*29) = 1.41; Hoeffding radius for 30 equal weights 0.992",
      abs(7 * 4 * np.log(80) / 87 - 1.41) < .005 and abs(np.sqrt(8 * np.log(40) / 30) - .992) < .001)


# ------------------------------------------------------------------------------------------
# 3. E32: simultaneous Hoeffding event over J preselected rows (two-sided per row) for averages of
#    L_it, U_it, T_it in [0,1]; probability of the event >= 1 - alpha under common within-precinct
#    shocks across elections (2x2 precincts with binomial-plus-shock totals).
# ------------------------------------------------------------------------------------------
def e32_event(K=30, reps=20000, alpha=.05):
    xs = rng.uniform(.1, .9, K)
    Lf = lambda t, x: np.maximum(0, (t - 1 + x) / x)
    Uf = lambda t, x: np.minimum(1, t / x)
    # law: precinct-wide shock s ~ U(-.3,.3) shared by both elections
    base = rng.uniform(.2, .8, K)
    def draw(n):
        s = rng.uniform(-.3, .3, (n, K))
        T0 = np.clip(base + s, 0, 1); T1 = np.clip(base + s + .1 * (xs - .5), 0, 1)
        return T0, T1
    big0, big1 = draw(200000)
    H = {"L0": Lf(big0, xs), "U0": Uf(big0, xs), "T0": big0, "L1": Lf(big1, xs), "U1": Uf(big1, xs), "T1": big1}
    EH = {k: v.mean(0) for k, v in H.items()}
    rows = []
    for k in H:                                   # two weightings per (t,H): uniform and by x
        for h in (np.ones(K) / K, xs / xs.sum()):
            rows.append((k, h))
    J = len(rows)
    T0, T1 = draw(reps)
    Hs = {"L0": Lf(T0, xs), "U0": Uf(T0, xs), "T0": T0, "L1": Lf(T1, xs), "U1": Uf(T1, xs), "T1": T1}
    ok = np.ones(reps, bool)
    for k, h in rows:
        c = np.sqrt(.5 * np.log(2 * J / alpha) * np.sum(h ** 2))
        ok &= np.abs(Hs[k] @ h - EH[k] @ h) <= c
    return ok.mean(), J


pe, J = e32_event()
check("E32 [V] the simultaneous moment event (J rows, two-sided Hoeffding + union) holds with probability >= .95 "
      "under shocks shared across elections", pe >= .95, f"J = {J}, P(event) = {pe:.4f}")

print(f"\n{sum(RES)}/{len(RES)} PASS")
